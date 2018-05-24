#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
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

void get_top_fingerprints(const struct fingerprint &fp, vector<struct fingerprint> db, vector< pair<float, int> > results, int t) {
    int best_core_idx = 0;
    float best_core_s1 = 0;
    vector<float> fp_local_orie, fp_local_cohe, fp_local_freq;
    get_fingerprint_local_values(fp, fp_local_orie, fp_local_cohe, fp_local_freq);
    float fp_avg_orie = get_fingerprint_average_orientation(fp);
    float fp_avg_freq = get_fingerprint_average_frequency(fp);
    int n = db.size();
    for (int i=0 ; i<n ; i++) { 
        int current_id = db[i+1].id;
        vector<float> db_local_orie, db_local_cohe, db_local_freq;
        get_fingerprint_local_values(db[i], db_local_orie, db_local_cohe, db_local_freq);

        float s1 = calculate_s1(fp_local_orie, fp_local_cohe, db_local_orie, db_local_cohe);
        if (s1 > best_core_s1) {
            best_core_idx = i;
            best_core_s1 = s1;
        }

        // Last core for a fingerprint
        if (i<n-1 && (db[i+1].id%5 == 1)) {
            get_fingerprint_local_values(db[best_core_idx], db_local_orie, db_local_cohe, db_local_freq);
            float db_avg_o = get_fingerprint_average_orientation(db[best_core_idx]);
            float db_avg_f = get_fingerprint_average_frequency(db[best_core_idx]);

            float s2 = calculate_s2(fp_local_freq, db_local_freq);

            cout << "s2 " << s2 << endl;

            float s3 = calculate_s3(fp_avg_freq, db_avg_f);

            cout << "s3 " << s3 << endl;

            float s4 = calculate_s4(fp_avg_orie, db_avg_o);

            cout << "s4 " << s4 << endl;

            float s = calculate_s(s1,s2,s3,s4);

            cout << "s " << s << endl;
            results.push_back(make_pair(s, db[best_core_idx].id));
        }
    }
}

int main(int argc, char** argv) {
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
    cerr << count_fp << endl;

    vector<float> local_orie, local_cohe, local_freq;
    get_fingerprint_local_values(fp[0], local_orie, local_cohe, local_freq);
    float avg_o = get_fingerprint_average_orientation(fp[0]);
    float avg_f = get_fingerprint_average_frequency(fp[0]);

    // Read the database
    vector<struct fingerprint> db;
    int count_db = read_from_file(db, db_filename);
    cerr << count_db << endl;

    // Vector for sorting s value
    vector< pair<float, int> > best_matches;

    for (int i=0 ; i<count_db ; i++) {
        cout << "DB fingerprint ID " << db[i].id << endl;
        vector<float> db_local_orie, db_local_cohe, db_local_freq;
        get_fingerprint_local_values(db[i], db_local_orie, db_local_cohe, db_local_freq);
        float db_avg_o = get_fingerprint_average_orientation(db[0]);
        float db_avg_f = get_fingerprint_average_frequency(db[0]);

        float s1 = calculate_s1(local_orie, local_cohe, db_local_orie, db_local_cohe);

        cout << "s1 " << s1 << endl;

        float s2 = calculate_s2(local_freq, db_local_freq);

        cout << "s2 " << s2 << endl;

        float s3 = calculate_s3(avg_f, db_avg_f);

        cout << "s3 " << s3 << endl;

        float s4 = calculate_s4(avg_o, db_avg_o);

        cout << "s4 " << s4 << endl;

        float s = calculate_s(s1,s2,s3,s4);

        cout << "s " << s << endl;
        best_matches.push_back(make_pair(s, db[i].id));
    }
    sort(best_matches.rbegin(), best_matches.rend());
    cout << "\nBest matches\n";
    for (int i=0 ; i<best_matches.size() ; i++) {
        cout << "ID " << best_matches[i].second << "\t: " << best_matches[i].first << endl;
    }
    return 0;
}

/*int main() {
    // Example values
    // fp.tif CORE 104 152
    float l_orie_1[] = {0, 97.7305, 96.0037, 104.154, 104.27, 99.6121, 76.1126, 0, 0, 98.8438, 94.4882, 97.5149, 90.0553, 74.8645, 60.5941, 53.0763, 0, 0, 90.7041, 87.4423, 77.35, 64.5631, 57.3232, 0, 0, 86.3477, 77.309, 65.6758, 58.5865, 55.7453, 0, 0, 73.9436, 65.2782, 62.7022, 80.4487};
    float l_cohe_1[] = {0, 0.856395, 0.656309, 0.29507, 0.812482, 0.945255, 0.739501, 0, 0, 0.612811, 0.707421, 0.38111, 0.801206, 0.914547, 0.626329, 0.580194, 0, 0, 0.504631, 0.468045, 0.737631, 0.739284, 0.523021, 0, 0, 0.385378, 0.719162, 0.500714, 0.693436, 0.54261, 0, 0, 0.559602, 0.322596, 0.609416, 0.641271};
    float l_freq_1[] = {0, 0.181545, 0.259259, 0.16, 0.157895, 0.125, 0.222222, 0, 0, 0.183776, 0.272727, 0.15, 0.173913, 0.0555556, 0.285714, 0.142857, 0, 0, 0.158353, 0.0714286, 0.173913, 0.125, 0.1, 0, 0, 0.148289, 0.133333, 0.147289, 0.142857, 0.138313, 0, 0, 0.136441, 0.134879, 0.04, 0.110869};
    float avg_orie_1 = -0.125954;
    float avg_freq_1 = 0.142895;

    //fp.tif CORE 104 168
    float l_orie_2[] = {0, 98.8438, 94.4882, 97.5149, 90.0553, 74.8645, 60.5941, 0, 0, 0, 90.7041, 87.4423, 77.35, 64.5631, 57.3232, 53.6184, 0, 0, 86.3477, 77.309, 65.6758, 58.5865, 55.7453, 0, 0, 0, 73.9436, 65.2782, 62.7022, 62.2227, 0, 0, 80.4487, 70.5397, 68.9653, 0};
    float l_cohe_2[] = {0, 0.612811, 0.707421, 0.38111, 0.801206, 0.914547, 0.626329, 0, 0, 0, 0.504631, 0.468045, 0.737631, 0.739284, 0.523021, 0.554236, 0, 0, 0.385378, 0.719162, 0.500714, 0.693436, 0.54261, 0, 0, 0, 0.559602, 0.322596, 0.609416, 0.525523, 0, 0, 0.641271, 0.686779, 0.508802, 0};
    float l_freq_2[] = {0, 0.183776, 0.272727, 0.15, 0.173913, 0.0555556, 0.285714, 0, 0, 0, 0.158353, 0.0714286, 0.173913, 0.125, 0.1, 0.14295, 0, 0, 0.148289, 0.133333, 0.147289, 0.142857, 0.138313, 0, 0, 0, 0.136441, 0.134879, 0.04, 0.130762, 0, 0, 0.110869, 0.110092, 0.109249, 0};
    float avg_orie_2 = -0.125954;
    float avg_freq_2 = 0.142895;

    //fp2.tif CORE 104 136
    float l_orie_3[] = {92.8416, 82.703, 75.2967, 65.567, 59.879, 56.9597, 147.818, 0, 0, 78.7643, 70.4503, 63.0479, 59.4192, 58.2259, 60.4374, 69.9172, 0, 83.8517, 76.6113, 69.5137, 67.4031, 70.097, 74.9326, 0, 0, 89.5193, 81.9179, 80.0141, 86.7053, 93.444, 0, 0, 102.074, 100.163, 105.52, 0};
    float l_cohe_3[] = {0.556538, 0.511061, 0.570864, 0.42559, 0.713579, 0.648699, 0.62761, 0, 0, 0.591479, 0.520984, 0.505018, 0.758203, 0.895867, 0.716803, 0.40558, 0, 0.549759, 0.411015, 0.575119, 0.615031, 0.722703, 0.544702, 0, 0, 0.46137, 0.496499, 0.670199, 0.751085, 0.416534, 0, 0, 0.55776, 0.496624, 0.600075, 0};
    float l_freq_3[] = {0.141691, 0.139967, 0.0833333, 0.143637, 0.111111, 0.25, 0.153846, 0, 0, 0.142564, 0.0833333, 0.111111, 0.140395, 0.0909091, 0.0869565, 0.13509, 0, 0.125677, 0.130276, 0.130836, 0.131476, 0.111111, 0.117647, 0, 0, 0.111715, 0.115348, 0.117691, 0.119609, 0.125048, 0, 0, 0.101149, 0.101822, 0.102356, 0};
    float avg_orie_3 = 1.22509;
    float avg_freq_3 = 0.130515;

    vector<float>orie_1(l_orie_1, l_orie_1 + ARRAY_SIZE);
    vector<float>cohe_1(l_cohe_1, l_cohe_1 + ARRAY_SIZE);
    vector<float>freq_1(l_freq_1, l_freq_1 + ARRAY_SIZE);

    vector<float>orie_2(l_orie_2, l_orie_2 + ARRAY_SIZE);
    vector<float>cohe_2(l_cohe_2, l_cohe_2 + ARRAY_SIZE);
    vector<float>freq_2(l_freq_2, l_freq_2 + ARRAY_SIZE);

    vector<float>orie_3(l_orie_3, l_orie_3 + ARRAY_SIZE);
    vector<float>cohe_3(l_cohe_3, l_cohe_3 + ARRAY_SIZE);
    vector<float>freq_3(l_freq_3, l_freq_3 + ARRAY_SIZE);

    float s1 = calculate_s1(orie_1, cohe_1, orie_3, cohe_3);

    cout << "s1 " << s1 << endl;

    float s2 = calculate_s2(freq_1, freq_3);

    cout << "s2 " << s2 << endl;

    float s3 = calculate_s3(avg_freq_1, avg_freq_3);

    cout << "s3 " << s3 << endl;

    float s4 = calculate_s4(avg_orie_1, avg_orie_3);

    cout << "s4 " << s4 << endl;

    float s = calculate_s(s1,s2,s3,s4);

    cout << "s " << s << endl;

    return 0;
}*/

// g++ -o indexing indexing.cpp fingerprint_structure.cpp