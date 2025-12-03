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

#include "../headers/SmithWaterman.h"
#include <iomanip>

int SWmatrix(const query& query, const Protein& prot,const Blosum& blosum, const int GOP, const int GEP) {
    
    string prot_sequence = prot.getseq();
    int prot_len = prot_sequence.size();
    int query_len = query.sequence.size();
    
    vector<int> H((query_len + 1)*(prot_len + 1), 0); 
    vector<int> F = H;
    vector<int> E = H;
    
    int max_score = 0;
    int i_curr, i_left, i_up;
    
    for (int i = 1; i <= query_len; i++)
    {
        for (int j = 1; j <= prot_len; j++)
        {
            i_curr = i*(prot_len+1) + j;
            i_left = i*(prot_len+1) + (j-1);
            i_up = (i-1)*(prot_len+1) + j;
            
            E[i_curr] = max(H[i_left] - GOP - GEP, E[i_left] - GEP);
            
            F[i_curr] = max(H[i_up] - GOP - GEP, F[i_up] - GEP);
            
            H[i_curr] = max(0, 
                max(H[(i-1)*(prot_len+1) + (j-1)] + blosum.Score(query.sequence[i-1], prot_sequence[j-1]), 
                max(E[i_curr], F[i_curr])));

            if (H[i_curr] > max_score)
            {
                max_score = H[i_curr];
            }
        }
    }
    
    return max_score;
}