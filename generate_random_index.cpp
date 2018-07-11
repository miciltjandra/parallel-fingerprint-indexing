#include <iostream>
#include <string>
#include "fingerprint_structure.hpp"
using namespace std;

#define MAXSIZE 50000

void generate(int num_index, string filename) {
    int last_id = get_last_id_from_file(filename);
    struct fingerprint indexes[num_index];
    for (int i=0 ; i<num_index ; i++) {
        indexes[i].id = last_id+1+i;
        for (int j=0 ; j<36 ; j++) {
            indexes[i].local_orientation[j] = rand()%256;
            indexes[i].local_coherence[j] = rand()%256;
            indexes[i].local_frequency[j] = rand()%256;
        }
        indexes[i].avg_orientation = rand()%256;
        indexes[i].avg_frequency = rand()%256;
    }
    save_to_file(num_index, indexes, filename);
}

int main(int argc, char** argv) {
    if (argc < 3) {
        cerr << "Usage : ./generate_random_index number-of-index filename\n";
        return 0;
    }
    int num_index = stoi(argv[1]);
    string filename = argv[2];
    cerr << "Making " << num_index << " indexes and saving to " << filename << endl;
    
    while (num_index > 0) {
        int num = min(num_index, MAXSIZE);
        generate(num, filename);
        num_index -= num;
    }    
    /*for (int i=0 ; i<num_index ; i++) {
        print_fingerprint_struct(indexes[i]);
    }*/
    cerr << "Finished making " << get_last_id_from_file(filename) << " indexes\n";
    return 0;
}

// g++ -o generate_random_index generate_random_index.cpp fingerprint_structure.cpp -std=c++11