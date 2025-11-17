#include "../headers/blast.h"

vector<char> Aminoacid = {'-', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'X', 'Y', 'Z', 'U', '*', 'O', 'J'};

//fonction qui inverse un entier sur 32bits vu que la plupart des données sont en big endian et que le processeur s'attend à du little endian
uint32_t swap(uint32_t val) {
    return ((((val) & 0xff000000) >> 24)|
      (((val) & 0x00ff0000) >>  8) |
      (((val) & 0x0000ff00) <<  8) |
      (((val) & 0x000000ff) << 24));
}

//fonction qui return un objet dataPin contenant les info imporant d'un fichier .pin 
//( nb de séquence et offsets pr les fichier .psd et .phr )
dataPin read_pin(const string& filepin) {
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

    //selon la documentation officiel, des bytes NULL sont ajouté de sorte que le nombre de séquence 
    //commence à une offset multiple de 8 ce qui "augmente" le nb de bytes à skippe, 
    //on calcule donc ce nb de bytes NULL (padding) puis on skip jusqu'au nb de séquence avec seekg
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
    dataPin pindata;
    pindata.numberOfprot = n_sequences;
    pindata.header_offsets.resize(n_sequences + 1);
    pindata.sequence_offsets.resize(n_sequences + 1);

    // On lit les offsets des headers, (n_sequences + 1) valeurs de 4 octets chacune selon la documentation
    file.read(reinterpret_cast<char*>(&pindata.header_offsets[0]),
              (n_sequences + 1) * 4);
    //de même pr les séquences
    file.read(reinterpret_cast<char*>(&pindata.sequence_offsets[0]),
              (n_sequences + 1) * 4);

    // on inverse pr les même raisons que title_len
    for (int i = 0; i <= n_sequences; i++) {
        pindata.header_offsets[i] = swap(pindata.header_offsets[i]);
        pindata.sequence_offsets[i] = swap(pindata.sequence_offsets[i]);
    }
    return pindata;
}

string read_sequence(ifstream& file, const int a,const int b){
    file.seekg(a, ios::beg);
    int size = llabs(b - a) ;
    string sequence;
    char byte;
    int value;
    //on lit byte par byte puis on convertie les valeurs en lettrez avec le vector Aminoacid
    for(int i = 0; i < size - 1; i++){
        file.read(&byte, 1);
        value = static_cast<int>(static_cast<unsigned char>(byte));
        sequence += Aminoacid[value];
    }
    return sequence;
}

//fonction qui lit le fichier phr et qui renvoie l'identifint de la prot entre l'offset a et b
string read_header(ifstream& file, const int a, const int b) {
    file.seekg(a);
    int size = b - a;
    vector<unsigned char> buffer(size);
    file.read(reinterpret_cast<char*>(&buffer[0]), size);

    string header;
    for (size_t i = 0; i < buffer.size(); i++) {
        if (buffer[i] == 0x1A) { // 0x1A = VisibleString
            // longueur du texte (1 octet)
            int length = buffer[i + 1];
            string text;
            for (int j = 0; j < length; j++) {
                text += static_cast<char>(buffer[i + j + 2]);
            }
            // on ajoute ce texte à la sortie
            header += text + " ";
        }
    }
    header = header.substr(0, header.find(' '));
    return header;
}