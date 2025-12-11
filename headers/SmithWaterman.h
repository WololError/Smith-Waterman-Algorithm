#ifndef SMITHWATERMAN_H
#define SMITHWATERMAN_H

#include <vector>
#include <algorithm>

#include "fasta.h"
#include "Protein.h"
#include "blosum.h"

using namespace std;

int SWmatrix(const string& query, const Protein& prot,const Blosum& blosum, const int GOP, const int GEP);

#endif