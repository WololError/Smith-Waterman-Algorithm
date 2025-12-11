#ifndef PROTEIN_H
#define PROTEIN_H

#include <string>
#include <vector>
#include <queue>
#include <thread>
#include <mutex>

#include "blast.h"
#include "blosum.h"

struct query;


using namespace std;

class Protein{

private:
    string id;
    string sequence;
    int sw_score = 0;
    
public:
    static vector<Protein> initProtlist(const string& phrfile, const string& psqfile, const dataPin pin);
    const string& getseq() const;
    const string& getid() const;
    int getscore() const;
    
    bool operator< (const Protein &other) const {
            if (sw_score != other.sw_score) {
                return sw_score < other.sw_score;
            }
            return id > other.id;
        }

    bool operator>(const Protein &other) const {
            if (sw_score != other.sw_score){
                return sw_score > other.sw_score;
            }
        return id < other.id; 
    }

    bool operator==(const Protein &other) const {
        return sw_score == other.sw_score && id == other.id;
    }
    
    static priority_queue<Protein> initProtqueue(const string& phrfile, const string& psqfile, const dataPin& pin, const query& query, Blosum& blosum, int GEP, int GOP);
    

    static void computeSW(int start, int end, const query& query, const Blosum& blosum, const string& phrfile, const string& psqfile, const dataPin& pin, int GEP, int GOP, priority_queue<Protein>& thread_results);
    static priority_queue<Protein> initProtqueueMT(const string& phrfile, const string& psqfile, const dataPin& pin, const query& query, Blosum& blosum, int GEP, int GOP);
    static priority_queue<Protein> mergeQueues(vector<priority_queue<Protein>>& all_queues);
    
    static void print20best(priority_queue<Protein>& pq);
};

#endif