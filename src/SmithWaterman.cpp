#include "../headers/SmithWaterman.h"

pair<string, string> GetProtGapped(vector<vector<int>>& H, const query& query, const Protein& prot, const Blosum& blosum, int i, int j, int gap_penalty)
{
    /*Cette fonction permet de print les alignements avec les tirets qui représentent les gaps, permet principalement de debugger*/

    // Traceback du meilleur alignement, i et j correspondent au max.
    string rev_on_query; // On stocke les caractères en ordre inverse
    string rev_on_prot;
    bool traceback = true;
    
    while (H[i][j] > 0 && i > 0 && j > 0 )
    {
        if (H[i][j] == H[i-1][j-1] + blosum.Score(query.sequence[i-1], prot.sequence[j-1]))
        {
            rev_on_query = rev_on_query + query.sequence[i-1];
            rev_on_prot = rev_on_prot + prot.sequence[j-1];
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
            rev_on_prot = rev_on_prot + prot.sequence[j-1];
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
	const Protein& prot, Blosum& blosum, const int gap_open_penalty, const int gap_extension_penalty)
{
    /*Renvoie le score maximal d'alignement de la protéine à comparer grâce à l'algorithme de Smith-Waterman*/

	// Init de 0 dans toute la matrice (ajout d'une ligne et d'une colonne supp.)
	int prot_len = prot.sequence.size();
    int query_len = query.sequence.size();

	static vector<vector<int>> H;

    H.assign(query_len + 1, vector<int>(prot_len + 1));
	
    int max_score = 0;
    int i_max, j_max, i, j;;
	int gap_len = 0;
	int gap_penalty = gap_open_penalty + gap_extension_penalty; // JSP s'il faut ajouter malus si le gap est plus grand

    // Double boucle créant la matrice sw
    cout << prot.sequence << endl;
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

    // Traceback pour afficher et debug les tests mais sinon inutile 
    // pair<string, string> two_prots_gapped = GetProtGapped(H, query, prot, blosum, i_max, j_max, gap_penalty);

    return max_score;
}

vector<Protein> Top20Prot(vector<Protein>& proteins, 
	const query& query, Blosum& blosum, const int gap_open_penalty, const int gap_extension_penalty)
{
    /*Retourne les 20 meilleurs protéines dans un vecteur de Protein*/
	vector<Protein> bestprots;
    Protein init_prot;
    bestprots.push_back(init_prot);

	for (int i = 0; i < proteins.size(); i++)
	{
        cout << i;
        proteins[i].sw_score = SWmatrix(query, proteins[i], blosum, gap_open_penalty, gap_extension_penalty);
		if (proteins[i].sw_score > bestprots.back().sw_score)
        {
            if (bestprots.size() > 20)
            {
                sort(bestprots.begin(), bestprots.end());
                bestprots.pop_back();
            }
            bestprots.push_back(proteins[i]);
        }
	}
    cout << endl;

	
	return bestprots;
}
