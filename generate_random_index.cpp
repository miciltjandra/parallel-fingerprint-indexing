#include <iostream>
#include <string>
#include "fingerprint_structure.hpp"
using namespace std;

int main(int argc, char** argv) {
    if (argc < 3) {
        cerr << "Usage : ./generate_random_index number-of-index filename\n";
        return 0;
    }
    int num_index = stoi(argv[1]);
    string filename = argv[2];
    cerr << "Making " << num_index << " indexes and saving to " << filename << endl;
    struct fingerprint indexes[num_index];
    for (int i=0 ; i<num_index ; i++) {
        indexes[i].id = i+1;
        for (int j=0 ; j<36 ; j++) {
            indexes[i].local_orientation[j] = rand()%256;
            indexes[i].local_coherence[j] = rand()%256;
            indexes[i].local_frequency[j] = rand()%256;
        }
        indexes[i].avg_orientation = rand()%256;
        indexes[i].avg_frequency = rand()%256;
    }
    
    /*for (int i=0 ; i<num_index ; i++) {
        print_fingerprint_struct(indexes[i]);
    }*/

    save_to_file(num_index, indexes, filename);
}

// g++ -o generate_random_index generate_random_index.cpp fingerprint_structure.cpp -std=c++11