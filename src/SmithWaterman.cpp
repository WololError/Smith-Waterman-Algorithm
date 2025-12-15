#include "../headers/SmithWaterman.h"

/* Calcule le score d’alignement local Smith-Waterman entre une séquence requête
 * et une protéine donnée.
 * L’algorithme utilise une matrice de substitution BLOSUM ainsi que des pénalités
 * de gap avec ouverture et extension. L’implémentation est optimisée en mémoire
 * en ne conservant que deux lignes de la matrice principale.
 *
 * @param query_string Séquence de la protéine requête.
 * @param prot_sequence Séquence de la protéine de la base à comparer.
 * @param blosum Objet Blosum contenant la matrice de substitution.
 * @param GOP Gap Opening Penalty.
 * @param GEP Gap Extension Penalty.
 *
 * @return Score maximal d’alignement local entre la requête et la protéine.
 */
int SWmatrix(const string& query_string, const string& prot_sequence,const Blosum& blosum, const int GOP, const int GEP) {
    
    int prot_len = prot_sequence.size();
    
    int query_len = query_string.size();
    
    //matrices H, F, E stockées en une dimension pour performance
    //H : score principal, E : gap horizontal, F : gap vertical
    vector<int> H_prev(prot_len + 1, 0);
    vector<int> H_curr(prot_len + 1, 0);
    vector<int> E(prot_len + 1, 0);
    vector<int> F(prot_len + 1, 0);

    int max_score = 0;
    int diag_prev, temp, diag_score;

    for (int i = 1; i <= query_len; i++) {
        diag_prev = 0;
        for (int j = 1; j <= prot_len; j++) {
            temp = H_prev[j];

            E[j] = max(H_curr[j-1] - GOP - GEP, E[j-1] - GEP);

            F[j] = max(H_prev[j] - GOP - GEP, F[j] - GEP);

            diag_score = diag_prev + blosum.Score(query_string[i-1], prot_sequence[j-1]);

            H_curr[j] = max({0, diag_score, E[j], F[j]});

            if (H_curr[j] > max_score)
                max_score = H_curr[j];

            diag_prev = temp;
        }

        swap(H_prev, H_curr);
        fill(H_curr.begin(), H_curr.end(), 0);
    }

    
    return max_score;
}
