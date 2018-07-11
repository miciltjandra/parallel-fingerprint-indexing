#include "fingerprint_structure.hpp"

unsigned char orientation_to_byte(float orientation) {
    float fresult = orientation/orientation_unit;
    unsigned char result = (char)fresult;
    // printf("Result %f\n", fresult);
    // printf("Result %d\n", (int)result);
    // printf("Result %c\n", result);
    return result;
}

float byte_to_orientation(unsigned char c) {
    float result = orientation_unit*(int)c;
    // printf("Result %f\n", result);
    return result;
}

unsigned char coherence_to_byte(float coherence) {
    float fresult = coherence * coherence_unit;
    unsigned char result = (char)fresult;
    // printf("Result %f\n", fresult);
    // printf("Result %d\n", (int)result);
    // printf("Result %c\n", result);
    return result;
}

float byte_to_coherence(unsigned char c) {
    float result = (float)c/coherence_unit;
    // printf("Result %f\n", result);
    return result;
}

unsigned char period_to_byte(float period) {
    float fresult = period/period_unit;
    unsigned char result = (char)fresult;
    // printf("Result %f\n", fresult);
    // printf("Result %d\n", (int)result);
    // printf("Result %c\n", result);
    return result;
}

float byte_to_period(unsigned char c) {
    float result = period_unit*(int)c;
    // printf("Result %f\n", result);
    return result;
}

unsigned char frequency_to_byte(float frequency) {
    if (frequency == 0) {
        return period_to_byte(frequency);
    } else {
        return period_to_byte(1.0f/frequency);
    }
}

float byte_to_frequency(unsigned char c) {
    float result = byte_to_period(c);
    if (result == 0) return result;
    else return 1/result;
}

struct fingerprint make_fingerprint_struct(int id, std::vector<float> local_orientation, std::vector<float> local_coherence, std::vector<float> local_frequency, float avg_orie, float avg_freq) {
    struct fingerprint result;
    result.id = id;
    for (int i=0 ; i<36 ; i++) {
        result.local_orientation[i] = orientation_to_byte(local_orientation[i]);
        result.local_coherence[i] = coherence_to_byte(local_coherence[i]);
        result.local_frequency[i] = frequency_to_byte(local_frequency[i]);
    }
    result.avg_orientation = orientation_to_byte(avg_orie+90.0f);
    result.avg_frequency = frequency_to_byte(avg_freq);
    return result;
}

void print_fingerprint_struct(struct fingerprint fp) {
    printf("ID %d\n", fp.id);
    printf("Local orientation\n");
    for (int i=0 ; i<36 ; i++) {
        printf("%d ", (int)fp.local_orientation[i]);
    }
    printf("\n");

    printf("Local coherence\n");
    for (int i=0 ; i<36 ; i++) {
        printf("%d ", (int)fp.local_coherence[i]);
    }
    printf("\n");

    printf("Local frequency\n");
    for (int i=0 ; i<36 ; i++) {
        printf("%d ", (int)fp.local_frequency[i]);
    }
    printf("\n");

    printf("Avg orientation : %d\n", (int)fp.avg_orientation);
    printf("Avg frequency : %d\n", (int)fp.avg_frequency);
}

void get_fingerprint_local_values(struct fingerprint fp, std::vector<float> &local_orientation, std::vector<float> &local_coherence, std::vector<float> &local_frequency) {
    for (int i=0 ; i<36 ; i++) {
        local_orientation.push_back(byte_to_orientation(fp.local_orientation[i]));
        local_coherence.push_back(byte_to_coherence(fp.local_coherence[i]));
        local_frequency.push_back(byte_to_frequency(fp.local_frequency[i]));
    }
}

float get_fingerprint_average_orientation(struct fingerprint fp) {
    return byte_to_orientation(fp.avg_orientation);
}

float get_fingerprint_average_frequency(struct fingerprint fp) {
    return byte_to_frequency(fp.avg_frequency);
}

void save_to_file(int size, struct fingerprint fps[], std::string filename) {
    FILE *f;
    char fname[filename.length()+1];
    strcpy(fname, filename.c_str());
    f = fopen(fname, "ab+");

    fwrite(&fps[0], sizeof(struct fingerprint), size, f);

    // for (int i=0 ; i<size ; i++)
    //     print_fingerprint_struct(fps[i]);

    // for (int i=0 ; i<size ; i++) {
    //     fprintf(f, "%d", fps[i].id);
    //     for (int j=0 ; j<36 ; j++) {
    //         fprintf(f, "%c", fps[i].local_orientation[j]);
    //     }
    //     for (int j=0 ; j<36 ; j++) {
    //         fprintf(f, "%c", fps[i].local_coherence[j]);
    //     }
    //     for (int j=0 ; j<36 ; j++) {
    //         fprintf(f, "%c", fps[i].local_frequency[j]);
    //     }
    //     fprintf(f, "%c", fps[i].avg_orientation);
    //     fprintf(f, "%c", fps[i].avg_frequency);
    // }
    fclose(f);
}

int read_from_file(std::vector<struct fingerprint> &fps, std::string filename) {
    FILE *f;
    char fname[filename.length()+1];
    strcpy(fname, filename.c_str());
    f = fopen(fname, "rb");
    if (f == NULL) {
        fprintf(stderr, "\nError opening file %s\n", fname);
        return 0;
    } else {
        printf("Successful opening %sa\n", fname);
    }

    int count = 0;
    
    struct fingerprint fp;
    while(fread(&fp, sizeof(struct fingerprint), 1, f)) {
        count++;
        fps.push_back(fp);
    }
    fclose(f);
    return count;
}

int get_last_id_from_file(std::string filename) {
    FILE *f;
    char fname[filename.length()+1];
    strcpy(fname, filename.c_str());
    f = fopen(fname, "rb");
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
    struct fingerprint fp;
    fread(&fp, sizeof(struct fingerprint), 1, f);
    print_fingerprint_struct(fp);
    return fp.id;
}

int get_new_fingerprint_id(int last_id) {
    if (last_id%5 == 0) {
        return last_id+1;
    } else {
        return last_id+6-(last_id%5);
    }
}
