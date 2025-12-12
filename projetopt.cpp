#include <iostream>
#include <fstream>
#include <chrono>

#include "headers/fasta.h"
#include "headers/blast.h"
#include "headers/blosum.h"
#include "headers/Protein.h"
#include "headers/SmithWaterman.h"

int main(int argc, char** argv){
    auto start = chrono::high_resolution_clock::now();
    string fastafile = argv[1];
    string pinfile = string(argv[2]) + ".pin";
    string psqfile = string(argv[2]) + ".psq";
    string phrfile = string(argv[2]) + ".phr";
    string blosumfile = string(argv[3]);

    int GOP = atoi(argv[4]); // GOP = Gap Open Penalty
    int GEP = atoi(argv[5]); // GEP = Gap Extension Penalty

    query query;
    query.getIdandsequence(fastafile);

    dataPin pindata;
    pindata.read_pin(pinfile);
    
    Blosum scoring(blosumfile);
    

// YOUR CODE HERE

    priority_queue<Protein> best20Prot = Protein::initProtqueueMT(phrfile, psqfile, pindata, query,scoring,GEP,GOP);

    Protein::print20best(best20Prot);
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
    // std::cout << "Elapsed time: " << duration.count() << " ms\n";
    return 0;
}