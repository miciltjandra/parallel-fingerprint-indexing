#include <iostream>
#include <string>
#include "fingerprint_structure.hpp"
using namespace std;

#define MAXSIZE 50000

int main(int argc, char** argv) {
    if (argc < 2) {
        cerr << "Usage : ./edit_index_data filename\n";
        return 0;
    }
    string filename = argv[1];

    vector<struct fingerprint> loads;
    int count = read_from_file(loads, filename);
    
    int last_id = loads[count-1].id;
    cout << "Last id " << last_id << endl;

    int id;
    cout << "Insert ID\n";
    cin >> id;
    struct fingerprint indexes[1];
    indexes[0].id = last_id+1;
    cout << "Insert 36 local orientation\n";
    float buf;
    for (int j=0 ; j<36 ; j++) {
        cin >> buf;
        indexes[0].local_orientation[j] = orientation_to_byte(buf);
    }
    cout << "Insert 36 local coherence\n";
    for (int j=0 ; j<36 ; j++) {
        cin >> buf;
        indexes[0].local_coherence[j] = coherence_to_byte(buf);
    }
    cout << "Insert 36 local frequency\n";
    for (int j=0 ; j<36 ; j++) {
        cin >> buf;
        indexes[0].local_frequency[j] = frequency_to_byte(buf);
    }
    cout << "Insert avg orientation\n";
    cin >> buf;
    indexes[0].avg_orientation = orientation_to_byte(buf);
    cout << "Insert avg frequency\n";
    cin >> buf;
    indexes[0].avg_frequency = frequency_to_byte(buf);
    
    FILE *f;
    char fname[filename.length()+1];
    strcpy(fname, filename.c_str());
    f = fopen(fname, "ab+");
    if (f == NULL) {
        fprintf(stderr, "\nError opening file %s\n", fname);
        return 0;
    }
    fseek(f, 0, SEEK_END);
    printf("\nFTELL %ld\n\n", ftell(f));
    if (ftell(f) < sizeof(struct fingerprint)) {
        printf("File empty/corrupt\n");
        return 0;
    }
    fseek(f, -sizeof(struct fingerprint), SEEK_END);
    printf("\nFTELL %ld\n\n", ftell(f));
    // struct fingerprint fp;
    // fread(&fp, sizeof(struct fingerprint), 1, f);
    // print_fingerprint_struct(fp);
    // return fp.id;
    fwrite(&indexes[0], sizeof(struct fingerprint), 1, f);
    fclose(f);

    return 0;
}

// g++ -o edit_index_data edit_index_data.cpp fingerprint_structure.cpp -std=c++11