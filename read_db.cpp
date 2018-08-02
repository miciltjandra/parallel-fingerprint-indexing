#include <iostream>
#include <string>
#include "fingerprint_structure.hpp"
using namespace std;

int main(int argc, char** argv) {
    if (argc < 2) {
        cerr << "Usage : ./read_db filename\n";
        return 0;
    }
    string filename = argv[1];
    cerr << "Reading fingerprints indexes from " << filename << endl;

    vector<struct fingerprint> loads;
    int count = read_from_file(loads, filename);
    
    for (int i=0 ; i<count ; i++) {
        print_fingerprint_struct(loads[i]);
        if (loads[i].id != i+1) {
            cerr << "Array index " << i << " , fingerprint ID " << loads[i].id << endl;
            return 0;
        }
    }

    printf("LOADS size %d\n", (int)loads.size());
    printf("Loaded fingerprint : %d\n", count);
    return 0;
}

// g++ -o read_db read_db.cpp fingerprint_structure.cpp -std=c++11
