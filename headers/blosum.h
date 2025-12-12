#ifndef BLOSUM_H
#define BLOSUM_H

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <string>
#include <unordered_map>
#include <sstream>

using namespace std;

class Blosum {

private:
    vector<int> matrix;
    int size;
    unordered_map<char, int> indexMap;

    int parseBlosumSize(const string& blosumfile) const;
    unordered_map<char, int> parseIndexMap(const string& blosumfile) const;
    static vector<int> linetovector(string& line);
;
public:
    Blosum(const string& blosumfile);
    int Score(char acide1, char acide2) const;
    void printMatrix() const;

};



#endif
