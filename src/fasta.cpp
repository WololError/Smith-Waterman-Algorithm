#include "../headers/fasta.h"

/* Retourne la séquence de la protéine requête.
 *
 * @return Référence constante vers la séquence.
 */
const string& query::get_seq() const{
    return this->sequence;
}

/* Retourne l’identifiant de la protéine requête.
 *
 * @return Référence constante vers l’identifiant.
 */
const string& query::get_id() const{
    return this->id;
}

/* Extrait l’identifiant et la séquence d’une protéine
 * à partir d’un fichier FASTA ne contenant qu’une seule protéine.
 * Cette fonction initialise l’état interne de l’objet query.
 *
 * @param filefasta Chemin vers le fichier FASTA.
 */
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
    fichier.close();
    this->id = id;
    this->sequence = sequence;
}

/* Recherche la protéine requête dans un vecteur de protéines.
 * Si une protéine possède la même séquence que la requête,
 * son identifiant est affiché.
 *
 * @param v Vecteur contenant les protéines de la base.
 */void query::findquery(const vector<Protein>& v){
    for(long unsigned int i = 0; i < v.size(); i++){
		//si la séquence correspond à celle d'une protéine
        if (this->sequence == v[i].getseq()){
            cout << v[i].getid() << endl; 
            break;
        }
    }
}
