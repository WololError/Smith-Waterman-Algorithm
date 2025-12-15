#include "../headers/blosum.h"
/* Construit un objet Blosum à partir d’un fichier BLOSUM.
 * Le constructeur initialise entièrement l’objet en :
 *  - déterminant la dimension de la matrice,
 *  - associant chaque acide aminé à un index,
 *  - remplissant la matrice de scores de substitution.
 *
 * @param blosumfile Chemin vers le fichier BLOSUM.
 */
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

/* Retourne le score de substitution BLOSUM entre deux acides aminés.
 * Si un acide aminé n’est pas reconnu, il est remplacé par '*'.
 *
 * @param acide1 Premier acide aminé.
 * @param acide2 Second acide aminé.
 *
 * @return Score BLOSUM correspondant à la paire d’acides aminés.
 */
int Blosum::Score(char acide1, char acide2) const {
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

/* Calcule la dimension de la matrice BLOSUM.
 * La taille est déterminée en comptant les lignes utiles du fichier,
 * en ignorant les lignes vides et les commentaires.
 *
 * @param blosumfile Chemin vers le fichier BLOSUM.
 *
 * @return Dimension de la matrice.
 */
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

/* Construit une table de correspondance entre chaque acide aminé
 * et son index dans la matrice BLOSUM.
 *
 * @param blosumfile Chemin vers le fichier BLOSUM.
 *
 * @return Map associant chaque acide aminé à son index.
 */
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

/* Convertit une ligne contenant des valeurs numériques en vecteur d’entiers.
 *
 * @param line Ligne du fichier BLOSUM à convertir.
 *
 * @return Vecteur contenant les valeurs de la ligne.
 */
vector<int> Blosum::linetovector(string& line) {

    istringstream iss(line);

    int num;
    vector<int> vectorofnum;

	//extrait chaque entier sur la ligne
    while (iss >> num)
        vectorofnum.push_back(num);

    return vectorofnum;
}

/* Affiche la matrice BLOSUM sur la sortie standard.
 */
void Blosum::printMatrix() const {
    
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            cout << " " << matrix[i * size + j] << " ";
        }
    cout << endl;
}
}
