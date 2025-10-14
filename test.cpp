#include <iostream>
#include <fstream>
#include <cstdint>
#include <stdexcept>
#include <vector>
#include <bit>
using namespace std;

inline uint32_t swap(uint32_t val) {
    return ((((val) & 0xff000000) >> 24)|
      (((val) & 0x00ff0000) >>  8) |
      (((val) & 0x0000ff00) <<  8) |
      (((val) & 0x000000ff) << 24));
}

struct Offsets {
    int numberOfprot;
    vector<uint32_t> header_offsets;
    vector<uint32_t> sequence_offsets;
};

Offsets read_offsets(const string& pin_path) {
    ifstream file(pin_path, ios::binary);
    if (!file) throw runtime_error("Impossible d'ouvrir le fichier");

    int32_t version, db_type, title_len, timestamp_len, n_sequences;

    file.read(reinterpret_cast<char*>(&version), 4);

    file.read(reinterpret_cast<char*>(&db_type), 4);

    file.read(reinterpret_cast<char*>(&title_len), 4);
    file.seekg(title_len, ios::cur);

    file.read(reinterpret_cast<char*>(&timestamp_len), 4);
    timestamp_len = swap(timestamp_len);
    int padding = (8 - (timestamp_len % 8)) % 8;
    file.seekg(timestamp_len + padding, ios::cur);

    file.read(reinterpret_cast<char*>(&n_sequences), 4);
    n_sequences = swap(n_sequences);

    int64_t residue_count;
    file.read(reinterpret_cast<char*>(&residue_count), 8);

    int32_t max_seq_len;
    file.read(reinterpret_cast<char*>(&max_seq_len), 4);

    Offsets offsets;
    offsets.numberOfprot = n_sequences;
    offsets.header_offsets.resize(n_sequences + 1);
    offsets.sequence_offsets.resize(n_sequences + 1);

    file.read(reinterpret_cast<char*>(offsets.header_offsets.data()),
              (n_sequences + 1) * 4);
    file.read(reinterpret_cast<char*>(offsets.sequence_offsets.data()),
              (n_sequences + 1) * 4);

    for (int i = 0; i <= n_sequences; i++) {
        offsets.header_offsets[i] = swap(offsets.header_offsets[i]);
        offsets.sequence_offsets[i] = swap(offsets.sequence_offsets[i]);
    }
    
    return offsets;
}

int main() {
    Offsets datapin = read_offsets("database/uniprot_sprot.fasta.pin");
    cout << datapin.numberOfprot << endl ;
}