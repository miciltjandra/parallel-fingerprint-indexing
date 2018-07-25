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
    
    FILE *f;
    char fname[filename.length()+1];
    strcpy(fname, filename.c_str());
    f = fopen(fname, "rb+");

    if (f == NULL) {
        fprintf(stderr, "\nError opening file %s\n", fname);
        return 0;
    } else {
        printf("Successful opening %sa\n", fname);
    }

    int count = 0;
    
    fread(&count, sizeof(int), 1, f);
    cout << "Count : " << count << endl;

    // Read Local Orientations
    unsigned char local_orientations[count*36];
    fread(local_orientations, sizeof(unsigned char), count*36, f);
    for (int i=0 ; i<count ; i++) {
        for (int j=0 ; j<36 ; j++) {
            cout << byte_to_orientation(local_orientations[i*36+j]) << " ";
        }
        cout << endl;
    }

    // Read Local Coherences
    unsigned char local_coherences[count*36];
    fread(local_coherences, sizeof(unsigned char), count*36, f);
    for (int i=0 ; i<count ; i++) {
        for (int j=0 ; j<36 ; j++) {
            cout << byte_to_coherence(local_coherences[i*36+j]) << " ";
        }
        cout << endl;
    }

    // Read Local Frequencies
    unsigned char local_frequencies[count*36];
    fread(local_frequencies, sizeof(unsigned char), count*36, f);
    for (int i=0 ; i<count ; i++) {
        for (int j=0 ; j<36 ; j++) {
            cout << byte_to_frequency(local_frequencies[i*36+j]) << " ";
        }
        cout << endl;
    }

    // Read Global Orientations
    unsigned char global_orientations[count];
    fread(global_orientations, sizeof(unsigned char), count, f);
    for (int i=0 ; i<count ; i++) {
        cout << byte_to_orientation(global_orientations[i]) << endl;
    }

    // Read Global Frequencies
    unsigned char global_frequencies[count];
    fread(global_frequencies, sizeof(unsigned char), count, f);
    for (int i=0 ; i<count ; i++) {
        cout << byte_to_frequency(global_frequencies[i]) << endl;
    }
    
    fclose(f);

    // cerr << "Finished making " << get_last_id_from_file(new_filename) << " indexes\n";
    return 0;
}

// g++ -o read_translated_structure read_translated_structure.cpp fingerprint_structure.cpp -std=c++11