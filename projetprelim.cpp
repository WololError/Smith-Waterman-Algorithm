#include <iostream>
#include <fstream>
#include <chrono>

#include "headers/fasta.h"
#include "headers/blast.h"
#include "headers/Protein.h"

int main(int argc, char** argv){

    string fastafile = argv[1];
    string pinfile = string(argv[2]) + ".pin";
    string psqfile = string(argv[2]) + ".psq";
    string phrfile = string(argv[2]) + ".phr";

    query query;
    query.getIdandsequence(fastafile);

    dataPin pindata;
    pindata.read_pin(pinfile);

    vector<Protein> proteins = Protein::initProtlist(phrfile, psqfile, pindata); 

    query.findquery(proteins);

    return 0;
}