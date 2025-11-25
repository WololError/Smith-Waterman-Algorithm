#ifndef SMITHWATERMAN_H
#define SMITHWATERMAN_H

#include <vector>
#include <algorithm>

#include "fasta.h"
#include "Protein.h"
#include "blosum.h"

using namespace std;

int SWmatrix(const query& query, 
	const Protein& prot, Blosum& blosum, const int gap_open_penalty, const int gap_extension);

#endif