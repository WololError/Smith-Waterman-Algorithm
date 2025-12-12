#ifndef BLAST_H
#define BLAST_H

#include <iostream>
#include <fstream>
#include <cstdint>
#include <vector>

using namespace std;

// structure pr stocker les informarions intéressante du fichier pin
class dataPin {

private:
    int numberOfprot;
    vector<uint32_t> header_offsets;
    vector<uint32_t> sequence_offsets;

public:
    int get_nop() const;
    const vector<uint32_t>& get_ho() const;
    const vector<uint32_t>& get_so() const;

    void read_pin(const string& pin_path);
};

uint32_t swapbytes(uint32_t val);
string read_sequence(ifstream& file, const int a,const int b);
string read_header(ifstream& file, const int a, const int b);

#endif