#include "../headers/SmithWaterman.h"
#include <iomanip>

pair<string, string> GetProtGapped(vector<vector<int>>& H, const query& query, const Protein& prot, const Blosum& blosum, int i, int j, int gap_penalty)
{
    /*Cette fonction permet de print les alignements avec les tirets qui représentent les gaps, permet principalement de debugger*/

    // Traceback du meilleur alignement, i et j correspondent au max.
    string prot_sequence = prot.getseq();
    string rev_on_query; // On stocke les caractères en ordre inverse
    string rev_on_prot;
    bool traceback = true;
    
    while (H[i][j] > 0 && i > 0 && j > 0 )
    {
        if (H[i][j] == H[i-1][j-1] + blosum.Score(query.sequence[i-1], prot_sequence[j-1]))
        {
            rev_on_query = rev_on_query + query.sequence[i-1];
            rev_on_prot = rev_on_prot + prot_sequence[j-1];
            i = i - 1;
            j = j - 1;
            continue;
        }
        if (H[i][j] == H[i-1][j] - gap_penalty)
        {
            rev_on_query = rev_on_query + query.sequence[i-1];
            rev_on_prot = rev_on_prot + "-";
            i = i - 1;
            continue;
        }
        if (H[i][j] == H[i][j-1] - gap_penalty)
        {
            rev_on_query = rev_on_query + "-";
            rev_on_prot = rev_on_prot + prot_sequence[j-1];
            j = j - 1;
            continue;
        }
        else
        {
            break;
        }
    }

    reverse(rev_on_query.begin(), rev_on_query.end()); // On inverse prcq c'était dans le mauvais sens
    reverse(rev_on_prot.begin(), rev_on_prot.end());

    return {rev_on_query, rev_on_prot};
}

int SWmatrix(const query& query, 
	const Protein& prot, Blosum& blosum, const int GOP, const int GEP)
{
    /*Renvoie le score maximal d'alignement de la protéine à comparer grâce à l'algorithme de Smith-Waterman*/

	// Init de 0 dans toute la matrice (ajout d'une ligne et d'une colonne supp.)
    string prot_sequence = prot.getseq();
	int prot_len = prot_sequence.size();
    int query_len = query.sequence.size();

	static vector<vector<int>> H;
    static vector<vector<int>> E;
    static vector<vector<int>> F;

    H.assign(query_len + 1, vector<int>(prot_len + 1));
	E = H;
    F = H;

    int max_score = 0;
    int i_max, j_max, i, j;

    // Double boucle créant la matrice sw grâce à l'implémentation GOTOH
    
    for (i = 1;  i <= query_len; i++)
	{
		for (j = 1; j <= prot_len; j++)
		{
			E[i][j] = max(H[i][j-1] - GOP - GEP, E[i][j-1] - GEP);
    
            F[i][j] = max(H[i-1][j] - GOP - GEP, F[i-1][j] - GEP);
     
            H[i][j] = max(0, 
				max(H[i-1][j-1] + blosum.Score(query.sequence[i-1], prot_sequence[j-1]), 
				max(E[i][j], F[i][j])));
    
            if (H[i][j] > max_score)
            {
                max_score = H[i][j];
                i_max = i;
                j_max = j;
            }
		}
	}
    // Traceback pour afficher et debug les tests mais sinon inutile 
    // pair<string, string> two_prots_gapped = GetProtGapped(H, query, prot, blosum, i_max, j_max, gap_penalty);
    
    return max_score;
}
