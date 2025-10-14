#include <iostream>
#include <fstream>
using namespace std;
struct Prot {
    string id;
    string sequence;
};

// fonction qui permet d'afficher sur le terminal les id suivi des sequence de tout les éléments d'un fichier fasta
void printfastafile(string files){
    ifstream fichier(files);
    if (!fichier.good()) {
        std::cerr << "fichier non ouvert" << std::endl;
        return ;
    }
    string ligne, id, sequence;
    while(getline(fichier, ligne)){
        
        if (ligne[0] == '>'){
            if (!id.empty()){                  // si on arrive ds une ligne qui commence par > mais qu'il y a déjà un id enregistré c'est qu'on a fini de lire tout lé séquence précedente
                cout << id << " : " << sequence << endl ;
            }
            id = ligne.erase(0,1);  // on récupère l'id
            sequence.clear();       // on réinitalise le sequence
        }
        
        else{
            sequence += ligne;     // on actualise le sequence car elle peut prendre plusieurs lignes
        }
    }
    if (!id.empty()) {
        cout << id << " : " << sequence << endl;   // on s'occupe du dernier id qui n'a pas été traité
}

}

// fonction qui permet de récupérer l'id et le sequence d'un fichier Fasta ne contenant qu'UN SEUL élement comme P00533.fasta
Prot getIdandsequence(string files){
    ifstream fichier(files);
    if (!fichier.good()) {
        std::cerr << "fichier non ouvert" << std::endl;
        return {"",""};
    }
    string ligne, id, sequence;
    while(getline(fichier, ligne)){
        
        if (ligne[0] == '>'){
            id = ligne.erase(0,1);  // on récupère l'id
            sequence.clear();       // on réinitalise le sequence
        }
        
        else{
            sequence += ligne;     // on actualise le sequence car elle peut prendre plusieurs lignes
        }
    }
    Prot query;
    query.id = id;
    query.sequence = sequence;
    return query;
}

int main(int argc, char **argv){
    // printfastafile(argv[1]);
    Prot query = getIdandsequence(argv[1]);
    return 0;
}