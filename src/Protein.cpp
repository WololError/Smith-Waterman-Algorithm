#include "../headers/Protein.h"
#include "../headers/blast.h"
#include "../headers/SmithWaterman.h"

//initialise la liste des protéines à partir des fichiers .phr et .psq
vector<Protein> Protein::initProtlist(const string& phrfile, const string& psqfile, const dataPin pin){
    
    ifstream phr(phrfile, ios::binary);
    if (!phr) throw runtime_error("Impossible d'ouvrir le fichier .phr");

    ifstream psq(psqfile, ios::binary);
    if (!psq) throw runtime_error("Impossible d'ouvrir le fichier .phr");

    vector<Protein> proteins(pin.get_nop());

    vector<uint32_t> sequence_offsets = pin.get_so();
    vector<uint32_t> header_offsets = pin.get_ho();

	//lecture des protéines via les offsets
    int i = 0;
    while (i < sequence_offsets.size() - 1){
        string p  = read_header(phr, header_offsets[i],header_offsets[i + 1]);
        string seq = read_sequence(psq, sequence_offsets[i], sequence_offsets[i + 1]);
        proteins[i].id = p;
        proteins[i].sequence = seq;
        i++;
    }

    phr.close();
    psq.close();
    return proteins;
}

const string& Protein::getseq() const{
    return this->sequence;
}

const string& Protein::getid() const{
    return this->id;
}

int Protein::getscore() const{
    return this->sw_score;
}

//initialise une priority_queue avec les protéines et leurs scores SW par rapport a une query
priority_queue<Protein> Protein::initProtqueue(const string& phrfile, const string& psqfile, const dataPin& pin, const query& query, Blosum& blosum, int GEP, int GOP){
    
    ifstream phr(phrfile, ios::binary);
    if (!phr) throw runtime_error("Impossible d'ouvrir le fichier .phr");

    ifstream psq(psqfile, ios::binary);
    if (!psq) throw runtime_error("Impossible d'ouvrir le fichier .psq");
    
    priority_queue<Protein> pq;
	
    vector<uint32_t> sequence_offsets = pin.get_so();
    vector<uint32_t> header_offsets = pin.get_ho();

	//lecture de chaque protéine et calcul du score SW
    int i = 0;
    while (i < sequence_offsets.size() - 1){    
        Protein P;
        P.id = read_header(phr, header_offsets[i],header_offsets[i + 1]);
        P.sequence = read_sequence(psq, sequence_offsets[i], sequence_offsets[i + 1]);
        P.sw_score = SWmatrix(query.get_seq(), P, blosum, GOP, GEP);
        pq.push(P);
        i++;
    }

    phr.close();
    psq.close();
    return pq;
}

//affiche les 20 meilleures protéines
void Protein::print20best(priority_queue<Protein>& pq){
    int i = 0;
    while (i < 20){
        const Protein& p = pq.top();
        cout << p.id << " " << p.sw_score << endl;
        pq.pop();
        i++;
    }
}

//calcule les scores SW pour un intervalle de protéines pour le multithreading
void Protein::computeSW(int start, int end, const query& query, const Blosum& blosum, const string& phrfile, const string& psqfile, const dataPin& pin, int GEP, int GOP, priority_queue<Protein>& thread_results) {
    const int TOP_K = 20;
    ifstream phr(phrfile, ios::binary);
    if (!phr) throw runtime_error("Impossible d'ouvrir le fichier .phr");
    ifstream psq(psqfile, ios::binary);
    if (!psq) throw runtime_error("Impossible d'ouvrir le fichier .psq");

    vector<uint32_t> sequence_offsets = pin.get_so();
    vector<uint32_t> header_offsets = pin.get_ho();

	//parcours des protéines assignées a ce thread
    int i = start;
    while (i < end) {
        Protein P;
        P.id = read_header(phr, header_offsets[i], header_offsets[i + 1]);
        P.sequence = read_sequence(psq, sequence_offsets[i], sequence_offsets[i + 1]);
        P.sw_score = SWmatrix(query.get_seq(), P, blosum, GOP, GEP);
        //on ne garde que les TOP_K protéines
        if (thread_results.size() < TOP_K) {
            thread_results.push(move(P));
        } else if (P.sw_score > thread_results.top().sw_score) {
            thread_results.pop(); 
            thread_results.push(move(P));
        }
        ++i;
    }
}

//initialise une priority_queue avec les meilleurs score de SW en multithreading
priority_queue<Protein> Protein::initProtqueueMT(const string& phrfile, const string& psqfile, const dataPin& pin, const query& query, Blosum& blosum, int GEP, int GOP) {
    unsigned int num_threads = thread::hardware_concurrency();
    int total_proteins = pin.get_nop();
    int chunk_size = total_proteins / num_threads;
    
    vector<priority_queue<Protein>> all_thread_results(num_threads);
    vector<thread> workers;

	//lancement d'un thread par portion de protéines
    unsigned int i = 0;
    while (i < num_threads) {
        int start_index = i * chunk_size;
        int end_index;

        if (i == num_threads - 1) { end_index = total_proteins; } else { end_index = start_index + chunk_size; }   

        if (start_index >= end_index){
            ++i;
            continue;
        } 
        //création du thread pour cette portion   
        workers.emplace_back(&Protein::computeSW, start_index, end_index, cref(query), cref(blosum), phrfile, psqfile, cref(pin), GEP, GOP, ref(all_thread_results[i]));
        ++i;
    }

	//attente de la fin des threads
    for (auto& worker : workers) {
        if (worker.joinable()) worker.join();
    }

    //fusion de toutes les queues 
    priority_queue<Protein> final_pq = mergeQueues(all_thread_results);
    return final_pq;
}

//fusionne plusieurs priority_queue en une seule
priority_queue<Protein> Protein::mergeQueues(vector<priority_queue<Protein>>& all_queues) {
    
    priority_queue<Protein> final_pq;
    
    for (auto& thread_pq : all_queues) {
        while (!thread_pq.empty()) {
            final_pq.push(thread_pq.top());
            thread_pq.pop();
        }
    }
    
    return final_pq;
}
