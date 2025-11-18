#include "../headers/SmithWaterman.h"
#include "../headers/fasta.h"
#include "../headers/blast.h"
#include "../headers/Protein.h"
#include "../headers/blosum.h"

vector<Protein> Top20Prot(const vector<Protein> proteins, 
	const query prot_to_compare, Blosum scoring, const int gap_open_penalty, const int gap_extension_penalty)
{
	vector<Protein> bestprots;
	vector<vector<int>> matrix;

	const int query_len = prot_to_compare.sequence.size();

	for (int i = 0; i < proteins.size(); i++)
	{
		matrix = SWmatrix(query_len, prot_to_compare, proteins[i]);
	}

	
	return bestprots;
}

vector<vector<int>> SWMatrix(int query_len, const query base, 
	Protein prot, Blosum blosum, const int gap_open_penalty, const int gap_extension_penalty)
{
	// Init de 0 dans toute la matrice (ajout d'une ligne et d'une colonne supp.)
	int prot_len = prot.sequence.size();
	vector<vector<int>> H(query_len + 1, vector<int>(prot_len + 1, 0));
	
	int gap_len = 0;
	int gap_penalty = gap_open_penalty + gap_extension_penalty(len);
	// Application de la formule de récurrence
	for (int i = 1;  i <= query_len; i++)
	{
		for (int j = 1; j <= prot_len; j++)
		{
			H[i][j] = max({0, 
				H[i-1][j-1] + blosum.Score(base[i], prot[j]), 
				H[i-1][j] - gap_penalty, 
				H[i][j-1] - gap_penalty }
			);

		}
	}


}
