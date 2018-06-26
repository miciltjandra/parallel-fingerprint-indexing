#include <stdio.h>
#include <iostream>
#include <chrono>
#include "fingerprint_structure.hpp"
using namespace std;

const int BLOCKSIZE = 36;

// Constant weights
const float w1 = 0.16f;
const float w2 = 0.37f;
const float w3 = 0.16f;
const float w4 = 0.31f;

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

__global__ void calculate_s1(fingerprint* db, fingerprint* fp, float* result) {
    __shared__ float ss, scos, ssin;
    int j = blockIdx.x;
    int i = threadIdx.x;
    // int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (i == 0) {
        ss = 0;
        scos = 0;
        ssin = 0;
    }
    __syncthreads();
    // (db+j)->local_frequency[i] = dfrequency_t_byte(dbyte_to_frequency((db+j)->local_frequency[i])+0.1);
    float s = dbyte_to_coherence(fp->local_coherence[i])*dbyte_to_coherence((db+j)->local_coherence[i]);
    float d = M_PI/180.0f * 2 * (dbyte_to_orientation(fp->local_orientation[i])-dbyte_to_orientation((db+j)->local_orientation[i]));
    float tcos = s*cos(d);
    float tsin = s*sin(d);

    atomicAdd(&ss, s);
    atomicAdd(&scos, tcos);
    atomicAdd(&ssin, tsin);
    __syncthreads();
    if (i == 0) {
        result[j] = sqrt(pow(scos,2)+pow(ssin,2))/ss;
    }
}

__global__ void calculate_s2(fingerprint* db, fingerprint* fp, float* result) {
    __shared__ float s_addition, s_absdiff;
    int j = blockIdx.x;
    int i = threadIdx.x;
    // int idx = blockIdx.x*blockDim.x + threadIdx.x;
    float t_addition = dbyte_to_frequency(fp->local_frequency[i]) + dbyte_to_frequency((db+j)->local_frequency[i]);
    float t_absdiff = abs(dbyte_to_frequency(fp->local_frequency[i]) - dbyte_to_frequency((db+j)->local_frequency[i]));
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

    vector<float> local_orie, local_cohe, local_freq;
    get_fingerprint_local_values(fp[0], local_orie, local_cohe, local_freq);
    float avg_o = get_fingerprint_average_orientation(fp[0]);
    float avg_f = get_fingerprint_average_frequency(fp[0]);

    // Read the database
    vector<struct fingerprint> db;
    int count_db = read_from_file(db, db_filename);
    cerr << "Fingerprint database count : " << count_db << endl;

    auto timer_start = chrono::steady_clock::now();
    // Test S1
    fingerprint *d_fp, *d_db;
    float s1_result[count_db], s2_result[count_db], s3_result[count_db], s4_result[count_db];
    float result[count_db];
    float *d_result;
    cudaMalloc((void **)&d_fp, sizeof(fingerprint));
    cudaMalloc((void **)&d_db, count_db*sizeof(fingerprint));
    cudaMalloc((void **)&d_result, count_db*sizeof(float));

    cudaMemcpy(d_db, &db[0], count_db*sizeof(fingerprint), cudaMemcpyHostToDevice);
    cudaMemcpy(d_fp, &fp[0], sizeof(fingerprint), cudaMemcpyHostToDevice);
    calculate_s1<<<count_db,BLOCKSIZE>>>(d_db, d_fp, d_result);
    cudaMemcpy(&s1_result[0], d_result, count_db*sizeof(float), cudaMemcpyDeviceToHost);

    // for (int i=0 ; i<count_db ; i++) {
    //     cout << i << " : ID " << db[i].id << endl;
    //     cout << "result = " << s1_result[i] << endl;
    // }

    // Test S2
    calculate_s2<<<count_db,BLOCKSIZE>>>(d_db, d_fp, d_result);
    cudaMemcpy(&s2_result[0], d_result, count_db*sizeof(float), cudaMemcpyDeviceToHost);
    // cout << "\n\nS2\n";
    // for (int i=0 ; i<count_db ; i++) {
    //     cout << i << " : ID " << db[i].id << endl;
    //     cout << "result = " << s2_result[i] << endl;
    // }

    // Test S3
    calculate_s3<<<count_db,1>>>(d_db, d_fp, d_result);
    cudaMemcpy(&s3_result[0], d_result, count_db*sizeof(float), cudaMemcpyDeviceToHost);
    // cout << "\n\nS3\n";
    // for (int i=0 ; i<count_db ; i++) {
    //     cout << i << " : ID " << db[i].id << endl;
    //     cout << "result = " << s3_result[i] << endl;
    // }

    // Test S4
    calculate_s4<<<count_db,1>>>(d_db, d_fp, d_result);
    cudaMemcpy(&s4_result[0], d_result, count_db*sizeof(float), cudaMemcpyDeviceToHost);
    // cout << "\n\nS4\n";
    // for (int i=0 ; i<count_db ; i++) {
    //     cout << i << " : ID " << db[i].id << endl;
    //     cout << "result = " << s4_result[i] << endl;
    // }

    // Test S
    // Copy S1-S4 to device
    float *d_s1_result, *d_s2_result, *d_s3_result, *d_s4_result;
    cudaMalloc((void **)&d_s1_result, count_db*sizeof(float));
    cudaMalloc((void **)&d_s2_result, count_db*sizeof(float));
    cudaMalloc((void **)&d_s3_result, count_db*sizeof(float));
    cudaMalloc((void **)&d_s4_result, count_db*sizeof(float));
    cudaMemcpy(d_s1_result, &s1_result[0], count_db*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_s2_result, &s2_result[0], count_db*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_s3_result, &s3_result[0], count_db*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_s4_result, &s4_result[0], count_db*sizeof(float), cudaMemcpyHostToDevice);
    calculate_s<<<1,count_db>>>(d_s1_result, d_s2_result, d_s3_result, d_s4_result, d_result);
    cudaMemcpy(&result[0], d_result, count_db*sizeof(float), cudaMemcpyDeviceToHost);
    cout << "\n\nS\n";
    for (int i=0 ; i<count_db ; i++) {
        cout << i << " : ID " << db[i].id << endl;
        cout << "result = " << result[i] << endl;
    }
    auto timer_end = chrono::steady_clock::now();
    chrono::duration<double> diff = timer_end - timer_start;
    cout << "Time to get indexing result for " << count_db << " fingerprints in DB : " << diff.count()  << endl;
    return 0;
}

// nvcc -o parallel_indexing parallel_indexing.cu fingerprint_structure.cpp -std=c++11
