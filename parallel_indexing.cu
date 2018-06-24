#include <stdio.h>
#include <iostream>
#include "fingerprint_structure.hpp"
using namespace std;

const int BLOCKSIZE = 36;
// Constant weights
const float w1 = 0.16f;
const float w2 = 0.37f;
const float w3 = 0.16f;
const float w4 = 0.31f;

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

__global__ void calculate_s1(fingerprint* db, fingerprint* fp, float* result, float* hs, float* hcos, float* hsin) {
    __shared__ float ss, scos, ssin;
    int j = blockIdx.x;
    int i = threadIdx.x;
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (i == 0) {
        ss = 0;
        scos = 0;
        ssin = 0;
    }
    __syncthreads();
    // (db+j)->local_frequency[i] = dfrequency_t_byte(dbyte_to_frequency((db+j)->local_frequency[i])+0.1);
    hs[idx] = i;
    float s = dbyte_to_coherence(fp->local_coherence[i])*dbyte_to_coherence((db+j)->local_coherence[i]);
    float d = M_PI/180.0f * 2 * (dbyte_to_orientation(fp->local_orientation[i])-dbyte_to_orientation((db+j)->local_orientation[i]));
    float tcos = s*cos(d);
    float tsin = s*sin(d);

    hs[idx] = tsin;
    atomicAdd(&ss, s);
    atomicAdd(&scos, tcos);
    atomicAdd(&ssin, tsin);
    __syncthreads();
    if (i == 0) {
        hcos[j] = scos;
        hsin[j] = ssin;
        hs[j] = ss;
        result[j] = sqrt(pow(scos,2)+pow(ssin,2))/ss;
    }
}

__global__ void calculate_s2(fingerprint* db, fingerprint* fp, float* result, float* hs) {
    __shared__ float s_addition, s_absdiff;
    int j = blockIdx.x;
    int i = threadIdx.x;
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    float t_addition = dbyte_to_frequency(fp->local_frequency[i]) + dbyte_to_frequency((db+j)->local_frequency[i]);
    float t_absdiff = abs(dbyte_to_frequency(fp->local_frequency[i]) - dbyte_to_frequency((db+j)->local_frequency[i]));
    hs[idx] = t_absdiff;
    atomicAdd(&s_addition, t_addition);
    atomicAdd(&s_absdiff, t_absdiff);
    if (i == 0) {
        result[j] = 1 - (s_absdiff/s_addition);
    }
}

__global__ void calculate_s3(fingerprint* db, fingerprint* fp, float* result) {
    int j = blockIdx.x;
    result[j] = 1 - (abs(dbyte_to_frequency(fp->avg_frequency)-dbyte_to_frequency((db+j)->avg_frequency))/max(dbyte_to_frequency(fp->avg_frequency), dbyte_to_frequency((db+j)->avg_frequency)));
}

__global__ void calculate_s4(fingerprint* db, fingerprint* fp, float* result) {
    int j = blockIdx.x;
    result[j] = 1-(abs(dbyte_to_orientation(fp->avg_orientation)-dbyte_to_orientation((db+j)->avg_orientation))/M_PI);
}

__global__ void calculate_s(float* s1, float* s2, float*s3, float* s4, float* result) {
    int i = threadIdx.x;
    result[i] = w1*s1[i] + w2*s2[i] + w3*s3[i] + w4*s4[i];
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
    fingerprint r[3];
    float s1_result[3], s2_result[3], s3_result[3], s4_result[3];
    float result[3];
    float *d_result;
    float scos[3], ssin[3];
    float *s;
    s = new float[3*36];
    for (int i=0 ; i<108 ; i++) s[i] = -1;
    float *d_s, *d_scos, *d_ssin;
    cudaMalloc((void **)&d_fp, sizeof(fingerprint));
    cudaMalloc((void **)&d_db, 3*sizeof(fingerprint));
    cudaMalloc((void **)&d_result, 3*sizeof(float));
    cudaMalloc((void **)&d_s, 3*36*sizeof(float));
    cudaMalloc((void **)&d_scos, 3*sizeof(float));
    cudaMalloc((void **)&d_ssin, 3*sizeof(float));
    cudaMemcpy(d_db, &db[0], 3*sizeof(fingerprint), cudaMemcpyHostToDevice);
    cudaMemcpy(d_fp, &fp[0], sizeof(fingerprint), cudaMemcpyHostToDevice);
    calculate_s1<<<3,36>>>(d_db, d_fp, d_result, d_s, d_scos, d_ssin);
    cudaMemcpy(&s1_result[0], d_result, 3*sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(s, d_s, 3*36*sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(&scos[0], d_scos, 3*sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(&ssin[0], d_ssin, 3*sizeof(float), cudaMemcpyDeviceToHost);

    for (int i=0 ; i<108 ; i++) {
        cout << i << " " << s[i] << endl;
    }
    for (int i=0 ; i<3 ; i++) {
        cout << i << " : ID " << db[i].id << endl;
        cout << "s = " << s[i] << endl;
        cout << "scos = " << scos[i] << endl;
        cout << "ssin = " << ssin[i] << endl;
        cout << "result = " << s1_result[i] << endl;
    }

    // Test S2
    calculate_s2<<<3,36>>>(d_db, d_fp, d_result, d_s);
    cudaMemcpy(s, d_s, 3*36*sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(&s2_result[0], d_result, 3*sizeof(float), cudaMemcpyDeviceToHost);
    cout << "\n\nS2\n";
    for (int i=0 ; i<108 ; i++) {
        cout << i << " " << s[i] << endl;
    }
    for (int i=0 ; i<3 ; i++) {
        cout << i << " : ID " << db[i].id << endl;
        cout << "result = " << s2_result[i] << endl;
    }

    // Test S3
    calculate_s3<<<3,1>>>(d_db, d_fp, d_result);
    cudaMemcpy(&s3_result[0], d_result, 3*sizeof(float), cudaMemcpyDeviceToHost);
    cout << "\n\nS3\n";
    for (int i=0 ; i<3 ; i++) {
        cout << i << " : ID " << db[i].id << endl;
        cout << "result = " << s3_result[i] << endl;
    }

    // Test S4
    calculate_s4<<<3,1>>>(d_db, d_fp, d_result);
    cudaMemcpy(&s4_result[0], d_result, 3*sizeof(float), cudaMemcpyDeviceToHost);
    cout << "\n\nS4\n";
    for (int i=0 ; i<3 ; i++) {
        cout << i << " : ID " << db[i].id << endl;
        cout << "result = " << s4_result[i] << endl;
    }

    // Test S
    // Copy S1-S4 to device
    float *d_s1_result, *d_s2_result, *d_s3_result, *d_s4_result;
    cudaMalloc((void **)&d_s1_result, 3*sizeof(float));
    cudaMalloc((void **)&d_s2_result, 3*sizeof(float));
    cudaMalloc((void **)&d_s3_result, 3*sizeof(float));
    cudaMalloc((void **)&d_s4_result, 3*sizeof(float));
    cudaMemcpy(d_s1_result, &s1_result[0], 3*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_s2_result, &s2_result[0], 3*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_s3_result, &s3_result[0], 3*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_s4_result, &s4_result[0], 3*sizeof(float), cudaMemcpyHostToDevice);
    calculate_s<<<1,3>>>(d_s1_result, d_s2_result, d_s3_result, d_s4_result, d_result);
    cudaMemcpy(&result[0], d_result, 3*sizeof(float), cudaMemcpyDeviceToHost);
    cout << "\n\nS4\n";
    for (int i=0 ; i<3 ; i++) {
        cout << i << " : ID " << db[i].id << endl;
        cout << "result = " << result[i] << endl;
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
