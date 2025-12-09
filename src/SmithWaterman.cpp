#include "../headers/SmithWaterman.h"

//calcule le score de Smith-Walterman entre une query et une protéine
//blosum : matrice du substitution
//GOP : gap opening penalty, GEP : gap extension penalty
int SWmatrix(const query& query, const Protein& prot,const Blosum& blosum, const int GOP, const int GEP) {
    
    string prot_sequence = prot.getseq();
    int prot_len = prot_sequence.size();
    int query_len = query.sequence.size();
    
    //matrices H, F, E stockées en une dimension pour performance
    //H : score principal, E : gap horizontal, F : gap vertical
    vector<int> H((query_len + 1)*(prot_len + 1), 0); 
    vector<int> F = H;
    vector<int> E = H;
    
    int max_score = 0;
    int i_curr, i_left, i_up;
    
    //parcours de toutes les positions de la matrice
    for (int i = 1; i <= query_len; i++)
    {
        for (int j = 1; j <= prot_len; j++)
        {
			//indices dans la matrice
            i_curr = i*(prot_len+1) + j;
            i_left = i*(prot_len+1) + (j-1);
            i_up = (i-1)*(prot_len+1) + j;
            
            //calcul du score pour un gap horizontal
            E[i_curr] = max(H[i_left] - GOP - GEP, E[i_left] - GEP);
            
            //calcul du score pour un gap vertical
            F[i_curr] = max(H[i_up] - GOP - GEP, F[i_up] - GEP);
            
            //score principal de Smith-Walterman
            H[i_curr] = max(0, 
                max(H[(i-1)*(prot_len+1) + (j-1)] + blosum.Score(query.sequence[i-1], prot_sequence[j-1]), 
                max(E[i_curr], F[i_curr])));
			
			//mise à jour du score max
            if (H[i_curr] > max_score)
            {
                max_score = H[i_curr];
            }
        }
    }
    
    return max_score;
}
