#include <iostream>
#include <string>
#include "fingerprint_structure.hpp"
using namespace std;

#define MAXSIZE 50000

int main(int argc, char** argv) {
    if (argc < 2) {
        cerr << "Usage : ./read_translated_structure filename\n";
        return 0;
    }

    string filename = argv[1];
    int *ids;
    unsigned char *local_orientations = NULL, *local_coherences, *local_frequencies, *global_orientations, *global_frequencies;
    // cout << "Main address : " << (void*) local_orientations << endl;
    int count = read_translated_structure(filename, ids, local_orientations, local_coherences, local_frequencies, global_orientations, global_frequencies);
    for (int i=0 ; i<count ; i++) {
        cout << ids[i] << endl;
    }

    // cout << "Main address : " << (void*) local_orientations << endl;

    for (int i=0 ; i<count ; i++) {
        for (int j=0 ; j<36 ; j++) {
            cout << byte_to_orientation(local_orientations[i*36+j]) << " ";
        }
        cout << endl;
    }

    for (int i=0 ; i<count ; i++) {
        for (int j=0 ; j<36 ; j++) {
            cout << byte_to_coherence(local_coherences[i*36+j]) << " ";
        }
        cout << endl;
    }

    for (int i=0 ; i<count ; i++) {
        for (int j=0 ; j<36 ; j++) {
            cout << byte_to_frequency(local_frequencies[i*36+j]) << " ";
        }
        cout << endl;
    }

    for (int i=0 ; i<count ; i++) {
        cout << byte_to_orientation(global_orientations[i]) << endl;
    }

    for (int i=0 ; i<count ; i++) {
        cout << byte_to_frequency(global_frequencies[i]) << endl;
    }

    delete[] ids;
    delete[] local_orientations;
    delete[] local_coherences;
    delete[] local_frequencies;
    delete[] global_orientations;
    delete[] global_frequencies;

    // cerr << "Finished making " << get_last_id_from_file(new_filename) << " indexes\n";
    return 0;
}

// g++ -o read_translated_structure read_translated_structure.cpp fingerprint_structure.cpp -std=c++11