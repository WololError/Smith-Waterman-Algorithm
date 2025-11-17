#include <iostream>
#include <fstream>
#include <chrono>
#include "headers/fasta.h"
#include "headers/blast.h"


int main(int argc, char** argv){

    string fastafile = argv[1];
    string pinfile = string(argv[2]) + ".pin";
    string psqfile = string(argv[2]) + ".psq";
    string phrfile = string(argv[2]) + ".phr";

    Prot query = getIdandsequence(fastafile);
    dataPin pindata = read_pin(pinfile);
    // cout << pindata.numberOfprot << endl;
    string seq;
    auto start = std::chrono::high_resolution_clock::now();
    ifstream file(psqfile, ios::binary);
    if (!file) throw runtime_error("Impossible d'ouvrir le fichier psq");
    
    for(int i = 0; i < pindata.sequence_offsets.size() - 1 ; i++){
        seq = read_sequence(file, pindata.sequence_offsets[i], pindata.sequence_offsets[i + 1]);
        if (seq == query.sequence){
            //cout << "protein Trouvee !" << " i = " << i << endl;
            cout << read_header(phrfile, pindata.header_offsets[i],pindata.header_offsets[i + 1]) << endl;
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end - start;
            // std::cout << "Temps ecoule : " << elapsed.count() << " secondes\n";
            break;
        }
    }
    file.close();
    return 0;
}