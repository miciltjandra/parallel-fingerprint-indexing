#include "fingerprint_structure.hpp"
#define ARRAY_SIZE 36
using namespace std;

int main() {
    float a = 104.3;
    printf("Orientation %f\n", a);
    unsigned char b = orientation_to_byte(a);
    float c = byte_to_orientation(b);

    float coherence = 0.67;
    printf("\nCoherence %f\n", coherence);
    unsigned char bc = coherence_to_byte(coherence);
    float cc = byte_to_coherence(bc);

    float period = 3.0f;
    printf("\nPeriod %f\n", period);
    unsigned char bp = period_to_byte(period);
    float cp = byte_to_period(bp);

    float freq = 1.0f/12;
    printf("\nFrequency %f\n", freq);
    unsigned char bf = frequency_to_byte(freq);
    float cf = byte_to_frequency(bf);

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

    struct fingerprint fingerprints[3];
    fingerprints[0] = make_fingerprint_struct(1, orie_1, cohe_1, freq_1, avg_orie_1, avg_freq_1);
    fingerprints[1] = make_fingerprint_struct(2, orie_2, cohe_2, freq_2, avg_orie_2, avg_freq_2);
    fingerprints[2] = make_fingerprint_struct(3, orie_3, cohe_3, freq_3, avg_orie_3, avg_freq_3);

    for (int i=0 ; i<3 ; i++)
        print_fingerprint_struct(fingerprints[i]);

    vector<float> local_orie, local_cohe, local_freq;
    get_fingerprint_local_values(fingerprints[0], local_orie, local_cohe, local_freq);
    printf("Orie\n");
    for (int i=0 ; i<36 ; i++) {
        printf("%f ", local_orie[i]);
    }
    printf("\nCohe\n");
    for (int i=0 ; i<36 ; i++) {
        printf("%f ", local_cohe[i]);
    }
    printf("\nFreq\n");
    for (int i=0 ; i<36 ; i++) {
        printf("%f ", local_freq[i]);
    }
    printf("\n");

    float avg_o = get_fingerprint_average_orientation(fingerprints[0]);
    float avg_f = get_fingerprint_average_frequency(fingerprints[0]);

    // printf("Avg orie %f\n", avg_o);
    // printf("Avg freq %f\n", avg_f);

    // printf("\n\nSAVING TO FILE\n");
    // save_to_file(3, fingerprints, "fingerprint_db");

    vector<struct fingerprint> loads;
    printf("\n\nLOADING FROM FILE\n");
    int count = read_from_file(loads, "fingerprint_db");
    printf("LOADS size %d\n", (int)loads.size());

    printf("Loaded fingerprint : %d\n", count);
    for (int i=0 ; i<count ; i++)
        print_fingerprint_struct(loads[i]);

    int last_id = get_last_id_from_file("fingerprint_db");
    printf("Last id %d\n", last_id);
    for (int i=0 ; i<14 ; i++)
        printf("Last id %d Next fp id %d\n", i, get_new_fingerprint_id(i));
    return 0;
}