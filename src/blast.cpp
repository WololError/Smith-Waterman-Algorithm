#include "../headers/blast.h"


vector<char> Aminoacid = {'-', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'X', 'Y', 'Z', 'U', '*', 'O', 'J'};

/* Retourne le vecteur contenant les offsets des en-têtes dans le fichier .phr.
 *
 * @return Référence constante vers le vecteur des offsets des headers.
 */
const vector<uint32_t>& dataPin::get_ho() const{
    return this->header_offsets;
}

/* Retourne le vecteur contenant les offsets des séquences dans le fichier .psq.
 *
 * @return Référence constante vers le vecteur des offsets des séquences.
 */
const vector<uint32_t>& dataPin::get_so() const{
    return this->sequence_offsets;
}

/* Retourne le nombre total de protéines contenues dans la base de données.
 *
 * @return Nombre de protéines.
 */
int dataPin::get_nop() const{
    return this->numberOfprot;
}

/* Inverse l’ordre des octets d’un entier non signé sur 32 bits.
 * Cette fonction est utilisée pour convertir des valeurs encodées
 * en big endian vers le format little endian attendu par le processeur.
 *
 * @param val Valeur 32 bits à inverser.
 * @return Valeur après inversion des octets.
 */
uint32_t swap(uint32_t val) {
    return ((((val) & 0xff000000) >> 24)|
      (((val) & 0x00ff0000) >>  8) |
      (((val) & 0x0000ff00) <<  8) |
      (((val) & 0x000000ff) << 24));
}

/* Lit un fichier .pin et extrait les informations essentielles de la base BLAST.
 * La fonction récupère notamment le nombre de séquences ainsi que les offsets
 * des en-têtes et des séquences utilisés pour accéder aux fichiers .phr et .psq.
 *
 * @param filepin Chemin vers le fichier .pin à lire.
 */
void dataPin::read_pin(const string& filepin) {
    // ouvre le fichier en lecture binaire
    ifstream file(filepin, ios::binary);
    if (!file) throw runtime_error("Impossible d'ouvrir le fichier");

    // read() lit directement les octets du fichier binaire à partir de la position actuelle du curseur,
    // et seekg() sert à déplacer ce curseur là où on veut, on l'utulisera pr skip des donné pas utiles

    //la version, le type, le taille du titre, la taille de date, le nb de séquence et la longueur max d'une séquence sont 
    //encodé en 32 bits selon la documentation officiel.
    int32_t version, db_type, title_len, timestamp_len, n_sequences,max_seq_len;

    //on lit la version, on se retrouve donc aprés les 4 1er bytes, on fait reinterpret_cast<char*> car .read attend par défaut un char*
    file.read(reinterpret_cast<char*>(&version), 4);

    //on lit la type, on se retrouve donc aprés le 8 premiers bytes
    file.read(reinterpret_cast<char*>(&db_type), 4);

    // on lit la taille du titre, on se retrouve aprés les 12 premiers bytes
    file.read(reinterpret_cast<char*>(&title_len), 4);

    //on a besoin de connaitre la taille du titre pr savoir cmb de byte skip, on inverse donc title_len 
    //car il est encodé en big endian alors que le processeur s'attend à lire du little endian
    title_len = swap(title_len);

    //on utulise donc seekg pr placer le pointeur aprés le titre
    file.seekg(title_len, ios::cur);

    //on lit la taille de la date de création et l'inverse pr les même raison que title_len
    file.read(reinterpret_cast<char*>(&timestamp_len), 4);
    timestamp_len = swap(timestamp_len);

    //on passe le timestamp
    file.seekg(timestamp_len, ios::cur);

    //on lit le nb de séquence et on l'inverse pr le récuperer plus tard
    file.read(reinterpret_cast<char*>(&n_sequences), 4);
    n_sequences = swap(n_sequences);

    //on skip le nb de résidue qui est encodé sur 64bits, il ne nous servira pas
    int64_t residue_count;
    file.read(reinterpret_cast<char*>(&residue_count), 8);

    //on skip la longueur max d'une séquence
    file.read(reinterpret_cast<char*>(&max_seq_len), 4);

    //on crée un objet dataPin pour stocker les informations extraites du fichier .pin.
    //On initialise son nombre de protéines et on redimensionne les vecteurs
    //d’offsets des en-têtes et des séquences. 

    this->numberOfprot = n_sequences;
    this->header_offsets.resize(n_sequences + 1);
    this->sequence_offsets.resize(n_sequences + 1);

    // On lit les offsets des headers, (n_sequences + 1) valeurs de 4 octets chacune selon la documentation
    file.read(reinterpret_cast<char*>(&this->header_offsets[0]),
              (n_sequences + 1) * 4);

    //de même pr les séquences
    file.read(reinterpret_cast<char*>(&this->sequence_offsets[0]),
              (n_sequences + 1) * 4);

    // on inverse pr les même raisons que title_len
    int i = 0;
    while (i <= n_sequences){
        this->header_offsets[i] = swap(this->header_offsets[i]);
        this->sequence_offsets[i] = swap(this->sequence_offsets[i]);
        i++;
    }
}

/* Lit une séquence protéique depuis un fichier .psq entre deux offsets.
 * Chaque octet lu est converti en acide aminé à l’aide du tableau Aminoacid.
 *
 * @param file Flux du fichier .psq ouvert.
 * @param a Offset de début de la séquence.
 * @param b Offset de fin de la séquence.
 *
 * @return Chaîne de caractères représentant la séquence protéique.
 */
string read_sequence(ifstream& file, const int a,const int b){
    file.seekg(a, ios::beg);
    int size = llabs(b - a) ;
    string sequence;
    char byte;
    int value;
    //on lit byte par byte puis on convertie les valeurs en lettrez avec le tableau Aminoacid
    int i = 0;
    while (i < size - 1) {
        file.read(&byte, 1);
        value = static_cast<int>(static_cast<unsigned char>(byte));
        sequence.push_back(Aminoacid[value]);
        i++;
    }
    return sequence;
}

/* Lit un en-tête protéique depuis un fichier .phr entre deux offsets.
 * La fonction extrait les champs de type VisibleString et retourne
 * l’identifiant principal de la protéine.
 *
 * @param file Flux du fichier .phr ouvert.
 * @param a Offset de début de l’en-tête.
 * @param b Offset de fin de l’en-tête.
 *
 * @return Identifiant de la protéine sous forme de chaîne de caractères.
 */
string read_header(ifstream& file, const int a, const int b) {
    file.seekg(a);
    int size = b - a;
    vector<unsigned char> buffer(size);
    file.read(reinterpret_cast<char*>(&buffer[0]), size);

    string header;
    size_t i = 0;
    while (i < buffer.size()) {
        if (buffer[i] == 0x1A) { // 0x1A = VisibleString
            // longueur du texte (1 octet)
            int length = buffer[i + 1];
            string text;
            header.append(reinterpret_cast<char*>(&buffer[i + 2]), length);
            header.push_back(' ');
            // on ajoute ce texte à la sortie
            header += text + " ";
        }
        i++;
    }
    header = header.substr(0, header.find(' '));
    return header;
}