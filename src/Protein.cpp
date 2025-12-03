#include "../headers/Protein.h"
#include "../headers/blast.h"
#include "../headers/SmithWaterman.h"

int m = 573661;
int c = 0;
vector<Protein> Protein::initProtlist(const string& phrfile, const string& psqfile, const dataPin pin){
    
    ifstream phr(phrfile, ios::binary);
    if (!phr) throw runtime_error("Impossible d'ouvrir le fichier .phr");

    ifstream psq(psqfile, ios::binary);
    if (!psq) throw runtime_error("Impossible d'ouvrir le fichier .phr");

    vector<Protein> proteins(pin.numberOfprot);

    for(int i = 0; i < pin.sequence_offsets.size() - 1 ; i++){
        string p  = read_header(phr, pin.header_offsets[i],pin.header_offsets[i + 1]);
        string seq = read_sequence(psq, pin.sequence_offsets[i], pin.sequence_offsets[i + 1]);
        proteins[i].id = p;
        proteins[i].sequence = seq;
    }

    phr.close();
    psq.close();
    return proteins;
}

string Protein::getseq() const{
    return this->sequence;
}

string Protein::getid() const{
    return this->id;
}

int Protein::getscore() const{
    return this->sw_score;
}

priority_queue<Protein> Protein::initProtqueue(const query& q, Blosum& blosum, string& phrfile, const string& psqfile, const dataPin& pin, int GEP, int GOP){
    
    ifstream phr(phrfile, ios::binary);
    if (!phr) throw runtime_error("Impossible d'ouvrir le fichier .phr");

    ifstream psq(psqfile, ios::binary);
    if (!psq) throw runtime_error("Impossible d'ouvrir le fichier .phr");
    
    priority_queue<Protein> pq;

    for(int i = 0; i < pin.sequence_offsets.size() - 1 ; i++){
        
        Protein P;
        P.id = read_header(phr, pin.header_offsets[i],pin.header_offsets[i + 1]);
        P.sequence = read_sequence(psq, pin.sequence_offsets[i], pin.sequence_offsets[i + 1]);
        P.sw_score = SWmatrix(q, P, blosum, GEP, GOP);
        
        cout << P.id + " " << P.sw_score << " || il en reste : " << pin.sequence_offsets.size() - i << endl;
        pq.push(P);
    }

    phr.close();
    psq.close();
    return pq;
}

void Protein::print20best(priority_queue<Protein>& pq){
    for(int i = 0; i < 20 ; i++){
        const Protein& p = pq.top();
        cout << p.id << " " << p.sw_score << endl;
        pq.pop();
    }
}

void Protein::computeSW(int start, int end, const query& query, const Blosum& blosum, const string& phrfile, const string& psqfile, const dataPin& pin, int GEP, int GOP, priority_queue<Protein>& thread_results) {
    const int TOP_K = 20;
    ifstream phr(phrfile, ios::binary);
    if (!phr) throw runtime_error("Impossible d'ouvrir le fichier .phr");
    ifstream psq(psqfile, ios::binary);
    if (!psq) throw runtime_error("Impossible d'ouvrir le fichier .psq");

    for (int i = start; i < end; ++i) {
        Protein P;
        P.id = read_header(phr, pin.header_offsets[i], pin.header_offsets[i + 1]);
        P.sequence = read_sequence(psq, pin.sequence_offsets[i], pin.sequence_offsets[i + 1]);
        P.sw_score = SWmatrix(query, P, blosum, GEP, GOP);
        if (thread_results.size() < TOP_K) {
            thread_results.push(move(P));
        } else if (P.sw_score > thread_results.top().sw_score) {
            thread_results.pop(); 
            thread_results.push(move(P));
        }
    }
}

vector<Protein> Protein::createVector(const string& phrfile, const string& psqfile, const dataPin& pin, const query& query, Blosum& blosum, int GEP, int GOP) {
    unsigned int num_threads = thread::hardware_concurrency();
    int total_proteins = pin.numberOfprot;
    int chunk_size = total_proteins / num_threads;
    const int TOP_K = 20;
    
    vector<priority_queue<Protein>> all_thread_results(num_threads);
    vector<thread> workers;

    for (unsigned int i = 0; i < num_threads; ++i) {
        int start_index = i * chunk_size;
        int end_index;

        if (i == num_threads - 1) { end_index = total_proteins; } else { end_index = start_index + chunk_size; }   

        if (start_index >= end_index) continue;
            
        workers.emplace_back(&Protein::computeSW, start_index, end_index, cref(query), cref(blosum), phrfile, psqfile, cref(pin), GEP, GOP, ref(all_thread_results[i]));
    }

    for (auto& worker : workers) {
        if (worker.joinable()) worker.join();
    }

    
    priority_queue<Protein> final_pq = mergeQueues(all_thread_results, TOP_K);
    
    vector<Protein> final_results;
    while (!final_pq.empty()) {
        final_results.push_back(final_pq.top());
        final_pq.pop();
    }
    
    return final_results;
}


void Protein::printbetter(vector<Protein>& v){
    for(int i = 0; i < 20 ; i++){
        cout << v[i].getid() << " " << v[i].getscore() << endl;
    }
}

priority_queue<Protein> Protein::mergeQueues(vector<priority_queue<Protein>>& all_queues, int top_k) {
    priority_queue<Protein> final_pq;
    
    for (auto& thread_pq : all_queues) {
        while (!thread_pq.empty()) {
            final_pq.push(thread_pq.top());
            thread_pq.pop();
        }
    }
    
    priority_queue<Protein> result;
    for (int i = 0; i < top_k && !final_pq.empty(); ++i) {
        result.push(final_pq.top());
        final_pq.pop();
    }
    
    return result;
}