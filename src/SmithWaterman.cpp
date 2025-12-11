#include "../headers/SmithWaterman.h"

//calcule le score de Smith-Walterman entre une query et une protéine
//blosum : matrice du substitution
//GOP : gap opening penalty, GEP : gap extension penalty
int SWmatrix(const string& query_string, const Protein& prot,const Blosum& blosum, const int GOP, const int GEP) {
    
    string prot_sequence = prot.getseq();
    int prot_len = prot_sequence.size();
    
    int query_len = query_string.size();
    
    //matrices H, F, E stockées en une dimension pour performance
    //H : score principal, E : gap horizontal, F : gap vertical
    vector<int> H((query_len + 1)*(prot_len + 1), 0); 
    vector<int> F = H;
    vector<int> E = H;
    
    int max_score = 0;
    int i_curr, i_left, i_up;
    
    //parcours de toutes les positions de la matrice
    int i = 1;
    while (i <= query_len)
    {
        int j = 1;
        while (j <= prot_len)
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
                max(H[(i-1)*(prot_len+1) + (j-1)] + blosum.Score(query_string[i-1], prot_sequence[j-1]), 
                max(E[i_curr], F[i_curr])));
			
			//mise à jour du score max
            if (H[i_curr] > max_score)
            {
                max_score = H[i_curr];
            }
            j++;
        }
        i++;
    }
    
    return max_score;
}
