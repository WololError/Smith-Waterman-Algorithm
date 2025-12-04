#ifndef FASTA_H
#define FASTA_H

#include <iostream>
#include <fstream>
#include <vector>
#include "Protein.h"
using namespace std;

//structure qui représente une protéine
struct query {
    string id;
    string sequence;
};

query getIdandsequence(const string& files);
void findquery(const query q, const vector<Protein>& v);
#endif