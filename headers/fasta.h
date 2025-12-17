#ifndef FASTA_H
#define FASTA_H

#include <iostream>
#include <fstream>
#include <vector>
#include "Protein.h"

using namespace std;

class query {
    
private:
    string id;
    string sequence;

public:
    const string& get_id() const;
    const string& get_seq() const;
    void getIdandsequence(const string& files);
    void findquery(const vector<Protein>& v);
};

#endif