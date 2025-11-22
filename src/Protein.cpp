#include "../headers/Protein.h"
#include "../headers/blast.h"
#include "../headers/SmithWaterman.h"

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

priority_queue<Protein> Protein::initProtqueue(const query& q, Blosum& blosum, string& phrfile, const string& psqfile, const dataPin& pin, int GEP, int GOP){
    
    ifstream phr(phrfile, ios::binary);
    if (!phr) throw runtime_error("Impossible d'ouvrir le fichier .phr");

    ifstream psq(psqfile, ios::binary);
    if (!psq) throw runtime_error("Impossible d'ouvrir le fichier .phr");
    
    priority_queue<Protein> pq;

    for(int i = 0; i < 20 ; i++){
        
        Protein P;
        P.id = read_header(phr, pin.header_offsets[i],pin.header_offsets[i + 1]);
        P.sequence = read_sequence(psq, pin.sequence_offsets[i], pin.sequence_offsets[i + 1]);
        P.sw_score = SWmatrix(q, P, blosum, GEP, GOP);
        
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