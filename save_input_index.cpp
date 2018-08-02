#include <iostream>
#include <string>
#include "fingerprint_structure.hpp"
using namespace std;

#define MAXSIZE 50000

int main(int argc, char** argv) {
    if (argc < 2) {
        cerr << "Usage : ./save_input_index filename\n";
        return 0;
    }
    string filename = argv[1];
    
    int last_id = get_last_id_from_file(filename);
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
    
    save_to_file(1, indexes, filename);

    cerr << "Finished saving ID " << get_last_id_from_file(filename) << " index\n";
    return 0;
}

// g++ -o save_input_index save_input_index.cpp fingerprint_structure.cpp -std=c++11