#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <chrono>
#include "fingerprint_structure.hpp"
using namespace std;

#define ARRAY_SIZE 36

// Constant weights
const float w1 = 0.16f;
const float w2 = 0.37f;
const float w3 = 0.16f;
const float w4 = 0.31f;

float calculate_s1(const vector<float> &local_orie_1, const vector<float> &local_coherence_1, const vector<float> &local_orie_2, const vector<float> &local_coherence_2) {
    float sum1 = 0.0f, sum2 = 0.0f, sum3 = 0.0f;
    for (int i=0 ; i<ARRAY_SIZE ; i++) {
        float s = local_coherence_1[i] * local_coherence_2[i];
        float d = local_orie_1[i] - local_orie_2[i];

        d = d * M_PI/180.0f * 2;

        sum1 += (s * cos(d));
        sum2 += (s * sin(d));
        sum3 += s;
    }
    float result = sqrt(pow(sum1,2)+pow(sum2,2))/sum3;

    return result;
}

float calculate_s2(const vector<float> &local_freq_1, const vector<float> &local_freq_2) {
    float sum1 = 0.0f, sum2 = 0.0f;
    for (int i=0 ; i<ARRAY_SIZE ; i++) {
        sum1 += abs(local_freq_1[i]-local_freq_2[i]);
        sum2 += local_freq_1[i]+local_freq_2[i];
    }
    float result = 1 - (sum1/sum2);
    return result;
}

float calculate_s3(const float &avg_freq_1, const float &avg_freq_2) {
    float result = 1 - (abs(avg_freq_1-avg_freq_2)/max(avg_freq_1, avg_freq_2));
    return result;
}

float calculate_s4(const float &local_orie_1, const float &local_orie_2) {
    float result = 1-(abs(local_orie_1-local_orie_2)/M_PI);
    return result;
}

float calculate_s(const float &s1, const float &s2, const float &s3, const float &s4) {
    float result = w1*s1 + w2*s2 + w3*s3 + w4*s4;
    return result;
}

void get_top_fingerprints(const struct fingerprint &fp, const vector<struct fingerprint> &db, vector< pair<float, int> > &results) {
    int best_core_idx = 0;
    float best_core_s1 = 0;
    vector<float> fp_local_orie, fp_local_cohe, fp_local_freq;
    get_fingerprint_local_values(fp, fp_local_orie, fp_local_cohe, fp_local_freq);
    
    float fp_avg_orie = get_fingerprint_average_orientation(fp);
    float fp_avg_freq = get_fingerprint_average_frequency(fp);
    int n = db.size();
    for (int i=0 ; i<n ; i++) {
        int current_id = db[i+1].id;
        vector<float> db_local_orie, db_local_cohe, stub;
        get_fingerprint_local_values(db[i], db_local_orie, db_local_cohe, stub);

        float s1 = calculate_s1(fp_local_orie, fp_local_cohe, db_local_orie, db_local_cohe);
        if (s1 > best_core_s1) {
            best_core_idx = i;
            best_core_s1 = s1;
        }

        // Last core for a fingerprint
        if (i==n-1 || (db[i+1].id%5 == 1)) {
            cout << "Best core " << best_core_idx << endl;
            vector<float> db_local_freq;
            get_fingerprint_local_values(db[best_core_idx], stub, stub, db_local_freq);

            float db_avg_o = get_fingerprint_average_orientation(db[best_core_idx]);
            float db_avg_f = get_fingerprint_average_frequency(db[best_core_idx]);

            float s2 = calculate_s2(fp_local_freq, db_local_freq);

            float s3 = calculate_s3(fp_avg_freq, db_avg_f);

            float s4 = calculate_s4(fp_avg_orie, db_avg_o);

            float s = calculate_s(best_core_s1,s2,s3,s4);

            results.push_back(make_pair(s, db[best_core_idx].id));
            best_core_idx = i+1;
            best_core_s1 = 0;
        }
    }
}

int main(int argc, char** argv) {
    cout << sizeof(struct fingerprint) << endl;
    if (argc < 3) {
        cerr << "Usage : ./indexing fingerprint-to-be-searched fingerprint-db\n";
        return 0;
    }
    string fp_filename = argv[1];
    string db_filename = argv[2];
    cerr << "FP " << fp_filename << " DB " << db_filename << endl;

    // Read the fingerprint to be searched
    vector<struct fingerprint> fp;
    int count_fp = read_from_file(fp, fp_filename);

    // Read the database
    vector<struct fingerprint> db;
    int count_db = read_from_file(db, db_filename);
    cerr << "Fingerprint database count : " << count_db << endl;

    // Start timer
    auto timer_start = chrono::steady_clock::now();

    // Read fingerprint FP value to local vectors
    vector<float> local_orie, local_cohe, local_freq;
    get_fingerprint_local_values(fp[0], local_orie, local_cohe, local_freq);
    float avg_o = get_fingerprint_average_orientation(fp[0]);
    float avg_f = get_fingerprint_average_frequency(fp[0]);

    // Vector for sorting s value
    vector< pair<float, int> > best_matches;

    // This is code for checking every core and assume every core is from different fingerprint
    /*for (int i=0 ; i<count_db ; i++) {
        // cout << "DB fingerprint ID " << db[i].id << endl;
        vector<float> db_local_orie, db_local_cohe, db_local_freq;
        get_fingerprint_local_values(db[i], db_local_orie, db_local_cohe, db_local_freq);

        float db_avg_o = get_fingerprint_average_orientation(db[i]);
        float db_avg_f = get_fingerprint_average_frequency(db[i]);

        float s1 = calculate_s1(local_orie, local_cohe, db_local_orie, db_local_cohe);

        float s2 = calculate_s2(local_freq, db_local_freq);

        float s3 = calculate_s3(avg_f, db_avg_f);

        float s4 = calculate_s4(avg_o, db_avg_o);

        float s = calculate_s(s1,s2,s3,s4);

        best_matches.push_back(make_pair(s, db[i].id));
    }
    sort(best_matches.rbegin(), best_matches.rend());
    cout << "\nBest matches\n";
    for (int i=0 ; i<best_matches.size() ; i++) {
        cout << "ID " << best_matches[i].second << "-"<< best_matches[i].second/5+(best_matches[i].second%5!=0) <<"\t: " << best_matches[i].first << endl;
    }

    auto timer_end = chrono::steady_clock::now();
    chrono::duration<double> diff = timer_end - timer_start;
    cout << "Time to get indexing result for " << count_db << " fingerprints in DB : " << diff.count()  << endl;
    */
   
    // This is code for checking that fingerprint may have more than 1 core
    // Core used is the one with best S1 value
    get_top_fingerprints(fp[0], db, best_matches);
    sort(best_matches.rbegin(), best_matches.rend());
    cout << "\nBest match\n";
    for (int i=0 ; i<best_matches.size() ; i++) {
        cout << "ID " << best_matches[i].second << "-"<< best_matches[i].second/5 <<"\t: " << best_matches[i].first;
        cout << endl;
    }
    auto timer_end = chrono::steady_clock::now();
    chrono::duration<double> diff = timer_end - timer_start;
    cout << "Time to get indexing result for " << count_db << " fingerprints in DB : " << diff.count()  << endl;
    
    return 0;
}

// g++ -o indexing indexing.cpp fingerprint_structure.cpp -std=c++11