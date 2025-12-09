#include "../headers/fasta.h"
//fonction qui lit une fichier fasta qui ne contient qu'une seul protéine et qui renvoie une objet de la structure Prot
//qui servira à la représenter
query getIdandsequence(const string& filefasta){
    ifstream fichier(filefasta);
    
    //verifie que le fichier est ouvert
    if (!fichier.good()) {
        std::cerr << "fichier fasta non ouvert" << std::endl;
        return {"",""};
    }

    string ligne, id, sequence;
    while(getline(fichier, ligne)){
        
        if (ligne[0] == '>'){
			//récupère l'id en enlevant l'entête
            id = ligne.erase(0,1);
        }
        else{
			//concatène les lignes de séquence
            sequence += ligne;
        }
    }
    //construit l'objet query
    fichier.close();
    query query;
    query.id = id;
    query.sequence = sequence;
    return query;
}

//fonction qui cherche une séquence donnée dans un vecteur de protéines
void findquery(const query q, const vector<Protein>& v){
    for(long unsigned int i = 0; i < v.size(); i++){
		//si la séquence correspond à celle d'une protéine
        if (q.sequence == v[i].getseq()){
            cout << v[i].getid() << endl; 
            break;
        }
    }
}
