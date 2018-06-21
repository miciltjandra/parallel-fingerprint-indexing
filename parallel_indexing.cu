#include <stdio.h>
#include <iostream>
#include "fingerprint_structure.hpp"
using namespace std;

const int BLOCKSIZE = 36;
__global__ void test(float *orie) {
    int i = threadIdx.x;
    orie[i] += i;
}

__host__ __device__ unsigned char dperiod_to_byte(float period) {
    float fresult = period/period_unit;
    unsigned char result = (char)fresult;
    return result;
}

__host__ __device__ float dbyte_to_period(unsigned char c) {
    float result = period_unit*(int)c;
    return result;
}

__host__ __device__ unsigned char dfrequency_to_byte(float frequency) {
    if (frequency == 0) {
        return dperiod_to_byte(frequency);
    } else {
        return dperiod_to_byte(1.0f/frequency);
    }
}

__host__ __device__ float dbyte_to_frequency(unsigned char c) {
    float result = dbyte_to_period(c);
    if (result == 0) return result;
    else return 1/result;
}

__device__ float dbyte_to_coherence(unsigned char c) {
    float result = (float)c/coherence_unit;
    return result;
}

__device__ float dbyte_to_orientation(unsigned char c) {
    float result = orientation_unit*(int)c;
    return result;
}

__global__ void teststruct(fingerprint* s) {
    int i = threadIdx.x;
    s->local_frequency[i] = dfrequency_to_byte(dbyte_to_frequency(s->local_frequency[i])+0.1);
}

__global__ void calculate_s1(fingerprint* db, fingerprint* fp, float* result) {
    int j = blockIdx.x;
    int i = threadIdx.x;
    // (db+j)->local_frequency[i] = dfrequency_to_byte(dbyte_to_frequency((db+j)->local_frequency[i])+0.1);
    __shared__ float ss, scos, ssin;
    float s = dbyte_to_coherence(fp->local_coherence[i])*dbyte_to_coherence((db+j)->local_coherence[i]);
    float d = M_PI/180.0f * 2 * (dbyte_to_orientation(fp->local_orientation[i])-dbyte_to_orientation((db+j)->local_orientation[i]));
    float tcos = s*cos(d);
    float tsin = s*sin(d);
    __syncthreads();
    if (i == 0) {
        ss = 0;
        scos = 0;
        ssin = 0;
    }
    __syncthreads();
    atomicAdd(&ss, s);
    atomicAdd(&scos, tcos);
    atomicAdd(&ssin, tsin);
    if (i == 0) {
        result[j] = sqrt(pow(scos,2)+pow(ssin,2))/s;
    }
}

int main(int argc, char** argv) {
    if (argc < 3) {
        cerr << "Usage : ./parallel_indexing fingerprint-to-be-searched fingerprint-db\n";
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

    //Test struct
    /*fingerprint f = fp[0];
    fingerprint r;
    fingerprint *d_f;
    cudaMalloc((void **)&d_f, sizeof(fingerprint));
    cudaMemcpy(d_f, &f, sizeof(f), cudaMemcpyHostToDevice);
    teststruct<<<1,36>>>(d_f);
    cudaMemcpy(&r, d_f, sizeof(f), cudaMemcpyDeviceToHost);
    
    for (int i=0 ; i<36 ; i++) {
        cout << dbyte_to_frequency(f.local_frequency[i]) << " " << dbyte_to_frequency(f.local_frequency[i]) << endl;
    }*/

    // Test S1
    fingerprint *d_fp, *d_db;
    d_fp = &fp[0];
    fingerprint r[3];
    float result[3];
    float *d_result;
    cudaMalloc((void **)&d_db, 3*sizeof(fingerprint));
    cudaMalloc((void **)&d_result, 3*sizeof(float));
    cudaMemcpy(d_db, &db[0], 3*sizeof(fingerprint), cudaMemcpyHostToDevice);
    calculate_s1<<<3,36>>>(d_db, d_fp, d_result);
    cudaMemcpy(&result[0], d_result, 3*sizeof(float), cudaMemcpyDeviceToHost);

    for (int i=0 ; i<3 ; i++) {
        cout << i << " : ID " << db[i].id << " " << result[i] << endl;
    }
    /*cudaMemcpy(&r, d_db, 3*sizeof(fingerprint), cudaMemcpyDeviceToHost);
    for (int j=0 ; j<3 ; j++) {
        cout << "\n\nFingerprint " << j << endl;
        for (int i=0 ; i<36 ; i++) {
            cout << dbyte_to_frequency(dfrequency_to_byte(dbyte_to_frequency(db[j].local_frequency[i])+0.1)) << " " << dbyte_to_frequency(r[j].local_frequency[i]) << endl;
        }
    }*/

    return 0;
}

// nvcc -o parallel_indexing parallel_indexing.cu fingerprint_structure.cpp
