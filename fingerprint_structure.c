#include "fingerprint_structure.h"

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

struct fingerprint make_fingerprint_struct(int id, float local_orientation[36], float local_coherence[36], float local_frequency[36], float avg_orie, float avg_freq) {
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

void get_fingerprint_local_values(struct fingerprint fp, float local_orientation[36], float local_coherence[36], float local_frequency[36]) {
    for (int i=0 ; i<36 ; i++) {
        local_orientation[i] = byte_to_orientation(fp.local_orientation[i]);
        local_coherence[i] = byte_to_coherence(fp.local_coherence[i]);
        local_frequency[i] = byte_to_frequency(fp.local_frequency[i]);
    }
}

float get_fingerprint_average_orientation(struct fingerprint fp) {
    return byte_to_orientation(fp.avg_orientation);
}

float get_fingerprint_average_frequency(struct fingerprint fp) {
    return byte_to_frequency(fp.avg_frequency);
}

void save_to_file(int size, struct fingerprint fps[]) {
    FILE *f;
    f = fopen("fingerprint_db", "ab+");

    fwrite(&fps[0], sizeof(struct fingerprint), size, f);

    for (int i=0 ; i<size ; i++)
        print_fingerprint_struct(fps[i]);

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

int read_from_file(struct fingerprint fps[]) {
    FILE *f;
    f = fopen("fingerprint_db", "rb");
    if (f == NULL) {
        fprintf(stderr, "\nError opening file\n");
        return 0;
    }

    int count = 0;
    
    while(fread(&fps[count], sizeof(struct fingerprint), 1, f)) {
        count++;
    }
    fclose(f);
    return count;
}
