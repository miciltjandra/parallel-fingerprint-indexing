#include <iostream>
#include <string>
#include "fingerprint_structure.hpp"
using namespace std;

#define MAXSIZE 50000

int main(int argc, char** argv) {
    if (argc < 3) {
        cerr << "Usage : ./translate_structure index-filename new-filename\n";
        return 0;
    }

    string filename = argv[1];
    string new_filename = argv[2];
    
    FILE *f;
    char fname[new_filename.length()+1];
    strcpy(fname, new_filename.c_str());
    f = fopen(fname, "wb+");

    // Read the fingerprint to be translated
    std::vector<struct fingerprint> fp;
    int count = read_from_file(fp, filename);
    cout << count << endl;

    fwrite(&count, sizeof(int), 1, f);

    // Write local orientation
    for (int i=0 ; i<count ; i++) {
        for (int j=0 ; j<36 ; j++) {
            cout << byte_to_orientation(fp[i].local_orientation[j]) << " ";
        }
        cout << endl;
        fwrite(fp[i].local_orientation, sizeof(unsigned char), 36, f);
    }

    // Write local coherence
    for (int i=0 ; i<count ; i++) {
        for (int j=0 ; j<36 ; j++) {
            cout << byte_to_coherence(fp[i].local_coherence[j]) << " ";
        }
        cout << endl;
        fwrite(fp[i].local_coherence, sizeof(unsigned char), 36, f);
    }

    // Write local frequency
    for (int i=0 ; i<count ; i++) {
        for (int j=0 ; j<36 ; j++) {
            cout << byte_to_frequency(fp[i].local_frequency[j]) << " ";
        }
        cout << endl;
        fwrite(fp[i].local_frequency, sizeof(unsigned char), 36, f);
    }

    // Write global orientation
    for (int i=0 ; i<count ; i++) {
        cout << byte_to_orientation(fp[i].avg_orientation) << endl;
        fwrite(&fp[i].avg_orientation, sizeof(unsigned char), 1, f);
    }

    // Write global frequency
    for (int i=0 ; i<count ; i++) {
        cout << byte_to_frequency(fp[i].avg_frequency) << endl;
        fwrite(&fp[i].avg_frequency, sizeof(unsigned char), 1, f);
    }

    fclose(f);

    cerr << "Finished making " << new_filename << endl;
    return 0;
}

// g++ -o translate_structure translate_structure.cpp fingerprint_structure.cpp -std=c++11