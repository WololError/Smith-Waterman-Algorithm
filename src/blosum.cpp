#include "../headers/blosum.h"
//construit la matrice à partir du fichier en paramètre
Blosum::Blosum(const string& blosumfile) {
	//lit la taille de la matrice
    this->size = parseBlosumSize(blosumfile);
    //associe chaque acide aminé a son index
    this->indexMap = parseIndexMap(blosumfile);

    this->matrix.resize(this->size * this->size);
	
    ifstream file(blosumfile);
    if (!file) throw runtime_error("Blosum() : impossible d'ouvrir le fichier");
	
    string line;
    int i = 0; 
	//remplit la matrice en parcourant le fichier
    while (getline(file, line)) {
        if (line.empty() || line[0] == '#' || line[0] == ' ' ) {
            continue; 
        }
        else {
			//remplace la lettre au début pour ne garder que la valeur numérique
            line[0] = ' ';
            vector<int> row = linetovector(line);
            for (int j = 0; j < size; j++) {
                this->matrix[i * size + j] = row[j];
            }
            i++;
        }
    } 
    file.close();
}

//retourne le score BLOSUM entre deux acide aminés
int Blosum::Score(char acide1, char acide2) const {
    //si un acide aminé est inconnue on remplace par *
    if (this->indexMap.count(acide1) == 0) {
        acide1 = '*';
    }
    
    if (this->indexMap.count(acide2) == 0) {
        acide2 = '*';
    }

    int i = indexMap.at(acide1);
    int j = indexMap.at(acide2);

    return matrix[i * size + j];
}

//calcule la dimension de la matrice en comptant les lignes utiles du fichier
int Blosum::parseBlosumSize(const string& blosumfile) const{
    ifstream file(blosumfile);
    if (!file) throw runtime_error("parseBlosumSize() : Impossible d'ouvrir le fichier");
    
    string line;
    int dimension = 0;
    while (getline(file, line)) {
        if (line.empty() || line[0] == '#' || line[0] == ' ' ) {
            continue;
        } else {
            dimension++;
        }
    } 
    file.close();
    return dimension;
}

//construit une map qui associe chaque acide aminé à son index dans la matrice
unordered_map<char, int> Blosum::parseIndexMap(const string& blosumfile) const{
    ifstream file(blosumfile);
    if (!file) throw runtime_error("parseIndexMap() : impossible d'ouvrir le fichier");
    
    unordered_map<char, int> map;
    string line;
    int i = 0;
    while (getline(file, line)) {
        if (line.empty() || line[0] == '#' || line[0] == ' ' ) {
            continue;
        } else {
			//le premier charactère est l'acide aminé
            map[line[0]] = i;
            i++;
        }
    } 
    file.close();
    return map;
}

//convertit une ligne contenant plusieurs entiers en un vecteur d'entiers
vector<int> Blosum::linetovector(string& line) {

    istringstream iss(line);

    int num;
    vector<int> vectorofnum;

	//extrait chaque entier sur la ligne
    while (iss >> num)
        vectorofnum.push_back(num);

    return vectorofnum;
}

//affiche la matrice BLOSUM
void Blosum::printMatrix() const {
    
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            cout << " " << matrix[i * size + j] << " ";
        }
    cout << endl;
}
}
