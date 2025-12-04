#include "../headers/blosum.h"

Blosum::Blosum(const string& blosumfile) {
    this->size = parseBlosumSize(blosumfile);
    this->indexMap = parseIndexMap(blosumfile);

    this->matrix.resize(this->size * this->size);

    ifstream file(blosumfile);
    if (!file) throw runtime_error("Blosum() : impossible d'ouvrir le fichier");

    string line;
    int i = 0; 

    while (getline(file, line)) {
        
        if (line.empty() || line[0] == '#' || line[0] == ' ' ) {
            continue; 
        }
        else {
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
            map[line[0]] = i;
            i++;
        }
    } 
    file.close();
    return map;
}

vector<int> Blosum::linetovector(string& line) {

    istringstream iss(line);

    int num;
    vector<int> vectorofnum;

    while (iss >> num)
        vectorofnum.push_back(num);

    return vectorofnum;
}

void Blosum::printMatrix() const {
    
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            cout << " " << matrix[i * size + j] << " ";
        }
    cout << endl;
}
}