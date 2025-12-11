#include "../headers/fasta.h"
//fonction qui lit une fichier fasta qui ne contient qu'une seul protéine et qui renvoie une objet de la structure Prot
//qui servira à la représenter

const string& query::get_seq() const{
    return this->sequence;
}

const string& query::get_id() const{
    return this->id;
}

void query::getIdandsequence(const string& filefasta){
    ifstream fichier(filefasta);
    
    //verifie que le fichier est ouvert
    if (!fichier.good()) {
        std::cerr << "fichier fasta non ouvert" << std::endl;
        this->id = "";
        this->sequence = "";
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
    this->id = id;
    this->sequence = sequence;
}

//fonction qui cherche une séquence donnée dans un vecteur de protéines
void query::findquery(const vector<Protein>& v){
    for(long unsigned int i = 0; i < v.size(); i++){
		//si la séquence correspond à celle d'une protéine
        if (this->sequence == v[i].getseq()){
            cout << v[i].getid() << endl; 
            break;
        }
    }
}
