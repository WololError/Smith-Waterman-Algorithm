#include "../headers/fasta.h"
//fonction qui lit une fichier fasta qui ne contient qu'un seul prot et qui renvoie une objet de la structure Prot
//qui servira à la représenter
query getIdandsequence(const string& filefasta){
    ifstream fichier(filefasta);
    
    if (!fichier.good()) {
        std::cerr << "fichier fasta non ouvert" << std::endl;
        return {"",""};
    }

    string ligne, id, sequence;
    while(getline(fichier, ligne)){
        
        if (ligne[0] == '>'){
            id = ligne.erase(0,1);  // on récupère l'id
        }
        
        else{
            sequence += ligne;     // on actualise le sequence car elle prend plusieurs lignes
        }
    }
    fichier.close();
    query query;
    query.id = id;
    query.sequence = sequence;
    return query;
}

void findquery(const query q, const vector<Protein> v){
    for(int i = 0; i < v.size(); i++){
        if (q.sequence == v[i].getseq()){
            cout << v[i].getid() << endl; 
            break;
        }
    }
}