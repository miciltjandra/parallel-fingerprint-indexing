#include <iostream>
#include <string>
#include "fingerprint_structure.hpp"
using namespace std;

#define MAXSIZE 50000

int read_translated_structure(string filename, int* &ids, unsigned char* &dlocal_orientations, unsigned char* &local_coherences, unsigned char* &local_frequencies, unsigned char* &global_orientations, unsigned char* &global_frequencies) {
    FILE *f;
    char fname[filename.length()+1];
    strcpy(fname, filename.c_str());
    f = fopen(fname, "rb+");

    if (f == NULL) {
        fprintf(stderr, "\nError opening file %s\n", fname);
        exit(1);
    } else {
        printf("Successful opening %sa\n", fname);
    }

    int count = 0;
    
    // Read count
    fread(&count, sizeof(int), 1, f);
    cout << "Count : " << count << endl;

    // Read IDs
    ids = new int[count];
    fread(ids, sizeof(int), count, f);

    // Read Local Orientations
    dlocal_orientations = new unsigned char[count*36];
    fread(dlocal_orientations, sizeof(unsigned char), count*36, f);

    // Read Local Coherences
    local_coherences = new unsigned char[count*36];
    fread(local_coherences, sizeof(unsigned char), count*36, f);

    // Read Local Frequencies
    local_frequencies = new unsigned char[count*36];
    fread(local_frequencies, sizeof(unsigned char), count*36, f);

    // Read Global Orientations
    global_orientations = new unsigned char[count];
    fread(global_orientations, sizeof(unsigned char), count, f);

    // Read Global Frequencies
    global_frequencies = new unsigned char[count];
    fread(global_frequencies, sizeof(unsigned char), count, f);

    // for (int i=0 ; i<count ; i++) {
    //     cout << ids[i] << endl;
    // }

    /*for (int i=0 ; i<count ; i++) {
        for (int j=0 ; j<36 ; j++) {
            cout << byte_to_orientation(local_orientations[i*36+j]) << " ";
        }
        cout << endl;
    }*/

    // cout << "Function address : " << (void*) dlocal_orientations << endl;

    /*for (int i=0 ; i<count ; i++) {
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
    }*/
    
    fclose(f);
    return count;
}

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