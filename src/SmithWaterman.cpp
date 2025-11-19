#include "../headers/SmithWaterman.h"
#include "../headers/fasta.h"
#include "../headers/blast.h"
#include "../headers/Protein.h"
#include "../headers/blosum.h"


int SWmatrix(query query, 
	Protein prot, Blosum blosum, const int gap_open_penalty, const int gap_extension_penalty)
{
	// Init de 0 dans toute la matrice (ajout d'une ligne et d'une colonne supp.)
	int prot_len = prot.sequence.size();
    int query_len = query.sequence.size();
	vector<vector<int>> H(query_len + 1, vector<int>(prot_len + 1, 0));
	
    int max_score = 0;
    int i_max, j_max, i, j;;
	int gap_len = 0;
	int gap_penalty = gap_open_penalty + gap_extension_penalty; // JSP s'il faut jouter malus si le gap est plus grand

    // Double boucle créant la matrice sw
    for (i = 1;  i <= query_len; i++)
	{
		for (j = 1; j <= prot_len; j++)
		{
			H[i][j] = max(0, 
				max(H[i-1][j-1] + blosum.Score(query.sequence[i-1], prot.sequence[j-1]), 
				max(H[i-1][j] - gap_penalty, 
				H[i][j-1] - gap_penalty))
			); // Pas de fonction max d'une liste, obligé de faire comme ça

            if (H[i][j] > max_score)
            {
                max_score = H[i][j];
                i_max = i;
                j_max = j;
            }
		}
	}

    // Traceback du meilleur alignement
    // string rev_on_query; // On stocke les caractères en ordre inverse
    // string rev_on_prot;
    // bool traceback = true;
    // i = i_max; // On redéfinit i et j comme la case où y a le max et ils vaudront tjs le max dans l'ordre
    // j = j_max;
    // while (H[i][j] > 0 && i > 0 && j > 0 )
    // {
    //     if (H[i][j] == H[i-1][j-1] + blosum.Score(query.sequence[i-1], prot.sequence[j-1]))
    //     {
    //         rev_on_query = rev_on_query + query.sequence[i-1];
    //         rev_on_prot = rev_on_prot + prot.sequence[j-1];
    //         i = i - 1;
    //         j = j - 1;
    //         continue;
    //     }
    //     if (H[i][j] == H[i-1][j] - gap_penalty)
    //     {
    //         rev_on_query = rev_on_query + query.sequence[i-1];
    //         rev_on_prot = rev_on_prot + "-";
    //         i = i - 1;
    //         continue;
    //     }
    //     if (H[i][j] == H[i][j-1] - gap_penalty)
    //     {
    //         rev_on_query = rev_on_query + "-";
    //         rev_on_prot = rev_on_prot + prot.sequence[j-1];
    //         j = j - 1;
    //         continue;
    //     }
    //     else
    //     {
    //         break;
    //     }
    // }

    return max_score;
}

vector<Protein> Top20Prot(vector<Protein> proteins, 
	const query query, Blosum blosum, const int gap_open_penalty, const int gap_extension_penalty)
{
	vector<Protein> bestprots;
	vector<vector<int>> matrix;

	for (int i = 0; i < proteins.size(); i++)
	{
        proteins[i].sw_score = SWmatrix(query, proteins[i], blosum, gap_open_penalty, gap_extension_penalty);
		if (proteins[i].sw_score > bestprots.back().sw_score)
        {
            if (i > 20)
            {
                bestprots.pop_back();
            }
            bestprots.push_back(proteins[i])
        }
	}

	
	return bestprots;
}
