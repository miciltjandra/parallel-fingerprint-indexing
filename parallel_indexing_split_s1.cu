#include <stdio.h>
#include <iostream>
#include <algorithm>
#include <chrono>
#include <thrust/sort.h>
#include "fingerprint_structure.hpp"

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

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

__global__ void calculate_s1_preparation(fingerprint* db, fingerprint* fp, float* s_sum, float* cos_sum, float* sin_sum) {
    __shared__ float ss, scos, ssin;
    int j = blockIdx.x;
    int i = threadIdx.x;
    if (i == 0) {
        ss = 0;
        scos = 0;
        ssin = 0;
    }
    __syncthreads();
    float s = dbyte_to_coherence(fp->local_coherence[i])*dbyte_to_coherence((db+j)->local_coherence[i]);
    float d = M_PI/180.0f * 2 * (dbyte_to_orientation(fp->local_orientation[i])-dbyte_to_orientation((db+j)->local_orientation[i]));
    float tcos = s*cos(d);
    float tsin = s*sin(d);

    atomicAdd(&ss, s);
    atomicAdd(&scos, tcos);
    atomicAdd(&ssin, tsin);
    __syncthreads();
    if (i == 0) {
        s_sum[j] = ss;
        cos_sum[j] = scos;
        sin_sum[j] = ssin;
    }
}

__global__ void calculate_s1(float* s_sum, float* cos_sum, float* sin_sum, float* result, int count) {
    int j = blockIdx.x*blockDim.x + threadIdx.x;
    if (j<count) result[j] = sqrt(pow(cos_sum[j],2)+pow(sin_sum[j],2))/s_sum[j];
}

__global__ void get_best_core_s1(fingerprint* db, float* result, int* mapping) {
    int i = blockIdx.x;
    if ((db+i)->id%5 == 1) {
        int max_idx = i;
        for (int j=1 ; j<5 ; j++) {
            if ((db+i+j)->id%5 == 1) break;
            else {
                if (result[i+j] > result[max_idx]) {
                    max_idx = i+j;
                }
            }
        }
        mapping[((db+i)->id-1)/5] = max_idx;
    }
}

__global__ void calculate_s2(fingerprint* db, fingerprint* fp, float* result, int* mapping) {
    __shared__ float s_addition, s_absdiff;
    int j = mapping[blockIdx.x];
    int i = threadIdx.x;
    if (i == 0) {
        s_addition = 0.0f;
        s_absdiff = 0.0f;
    }
    float t_addition = dbyte_to_frequency(fp->local_frequency[i]) + dbyte_to_frequency((db+j)->local_frequency[i]);
    float t_absdiff = abs(dbyte_to_frequency(fp->local_frequency[i]) - dbyte_to_frequency((db+j)->local_frequency[i]));
    atomicAdd(&s_addition, t_addition);
    atomicAdd(&s_absdiff, t_absdiff);
    __syncthreads();
    if (i == 0) {
        result[blockIdx.x] = 1 - (s_absdiff/s_addition);
    }
}

__global__ void calculate_s3(fingerprint* db, fingerprint* fp, float* result, int* mapping) {
    int j = mapping[blockIdx.x];
    result[blockIdx.x] = 1 - (abs(dbyte_to_frequency(fp->avg_frequency)-dbyte_to_frequency((db+j)->avg_frequency))/max(dbyte_to_frequency(fp->avg_frequency), dbyte_to_frequency((db+j)->avg_frequency)));
}

__global__ void calculate_s4(fingerprint* db, fingerprint* fp, float* result, int* mapping) {
    int j = mapping[blockIdx.x];
    result[blockIdx.x] = 1-(abs(dbyte_to_orientation(fp->avg_orientation)-dbyte_to_orientation((db+j)->avg_orientation))/180.0f);
}

__global__ void calculate_s(float* s1, float* s2, float*s3, float* s4, float* result, int* mapping) {
    int i = blockIdx.x;
    result[i] = w1*s1[mapping[i]] + w2*s2[i] + w3*s3[i] + w4*s4[i];
}

__global__ void get_top_fingerprints(float* s, float* result, int* mapping) {
    int i = threadIdx.x;
    result[i] = s[mapping[i]];
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage : ./parallel_indexing fingerprint-to-be-searched fingerprint-db\n";
        return 0;
    }

    std::string fp_filename = argv[1];
    std::string db_filename = argv[2];

    // Read the fingerprint to be searched
    std::vector<struct fingerprint> fp;
    int count_fp = read_from_file(fp, fp_filename);

    // Read the database
    std::vector<struct fingerprint> db;
    int count_db = read_from_file(db, db_filename);
    std::cerr << "Fingerprint core database count : " << count_db << std::endl;

    std::cerr << "Last fingerprint ID : " << db[count_db-1].id << std::endl;
    int count_db_fingerprint = (db[count_db-1].id-1)/5+1;
    std::cerr << "Fingerprint database count : " << count_db_fingerprint << std::endl;

    auto timer_start = std::chrono::steady_clock::now();

    // Preparing memory
    fingerprint *d_fp, *d_db;
    std::vector<float> result(count_db_fingerprint, 0);
    float *d_s1_result, *d_s2_result, *d_s3_result, *d_s4_result, *d_result;
    
    cudaMalloc((void **)&d_fp, sizeof(fingerprint));
    cudaMalloc((void **)&d_db, count_db*sizeof(fingerprint));

    cudaMalloc((void **)&d_s1_result, count_db*sizeof(float));
    cudaMalloc((void **)&d_s2_result, count_db_fingerprint*sizeof(float));
    cudaMalloc((void **)&d_s3_result, count_db_fingerprint*sizeof(float));
    cudaMalloc((void **)&d_s4_result, count_db_fingerprint*sizeof(float));
    cudaMalloc((void **)&d_result, count_db_fingerprint*sizeof(float));

    //Mapping for fingerprint to fingerprint core idx
    int *d_mapping;
    cudaMalloc((void **)&d_mapping, count_db_fingerprint*sizeof(int));

    cudaMemcpy(d_db, &db[0], count_db*sizeof(fingerprint), cudaMemcpyHostToDevice);
    cudaMemcpy(d_fp, &fp[0], sizeof(fingerprint), cudaMemcpyHostToDevice);

    //Additional Memory for S1
    float *d_s_sum, *d_cos_sum, *d_sin_sum;
    cudaMalloc((void **)&d_s_sum, count_db*sizeof(float));
    cudaMalloc((void **)&d_cos_sum, count_db*sizeof(float));
    cudaMalloc((void **)&d_sin_sum, count_db*sizeof(float));

    // S1
    calculate_s1_preparation<<<count_db,BLOCKSIZE>>>(d_db, d_fp, d_s_sum, d_cos_sum, d_sin_sum);
    calculate_s1<<<(count_db/256)+1, 256>>>(d_s_sum, d_cos_sum, d_sin_sum, d_s1_result, count_db);
    get_best_core_s1<<<count_db, 1>>>(d_db, d_s1_result, d_mapping);
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );

    std::vector<int> mapping(count_db_fingerprint, 0);
    cudaMemcpy(&mapping[0], d_mapping, count_db_fingerprint*sizeof(int), cudaMemcpyDeviceToHost);
    // std::vector<float> s1_result;
    // s1_result.resize(count_db, 0);
    // cudaMemcpy(&s1_result[0], d_s1_result, count_db*sizeof(float), cudaMemcpyDeviceToHost);

    // S2
    // Only calculate for 1 core per fingerprint using mapping
    calculate_s2<<<count_db_fingerprint,BLOCKSIZE>>>(d_db, d_fp, d_s2_result, d_mapping);
    // cudaMemcpy(&s2_result[0], d_s2_result, count_db_fingerprint*sizeof(float), cudaMemcpyDeviceToHost);

    // S3
    calculate_s3<<<count_db_fingerprint,1>>>(d_db, d_fp, d_s3_result,d_mapping);
    // cudaMemcpy(&s3_result[0], d_s3_result, count_db*sizeof(float), cudaMemcpyDeviceToHost);

    // S4
    calculate_s4<<<count_db_fingerprint,1>>>(d_db, d_fp, d_s4_result, d_mapping);
    // cudaMemcpy(&s4_result[0], d_s4_result, count_db*sizeof(float), cudaMemcpyDeviceToHost);

    // S
    calculate_s<<<count_db_fingerprint, 1>>>(d_s1_result, d_s2_result, d_s3_result, d_s4_result, d_result, d_mapping);
    // cudaMemcpy(&result[0], d_result, count_db_fingerprint*sizeof(float), cudaMemcpyDeviceToHost);

    // ID for identifying fingerprint during sort
    int* ids = new int[count_db_fingerprint];
    for (int i=0 ; i<count_db_fingerprint ; i++) {
        ids[i] = db[mapping[i]].id;
    }
    int* d_ids;
    cudaMalloc((void **)&d_ids, count_db_fingerprint*sizeof(int));
    cudaMemcpy(d_ids, &ids[0], count_db_fingerprint*sizeof(int), cudaMemcpyHostToDevice);
    
    auto sort_start = std::chrono::steady_clock::now();
    thrust::sort_by_key(thrust::device, d_result, d_result+count_db_fingerprint, d_ids);
    auto sort_end = std::chrono::steady_clock::now();

    cudaMemcpy(&result[0], d_result, count_db_fingerprint*sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(&ids[0], d_ids, count_db_fingerprint*sizeof(int), cudaMemcpyDeviceToHost);
    
    /*for (int i=count_db_fingerprint-1 ; i>=0 ; i--) {
        std::cout << "ID " << ids[i] << "-"<< ids[i]/5 <<"\t: " << result[i];
        std::cout << std::endl;
    }*/
    auto timer_end = std::chrono::steady_clock::now();
    std::chrono::duration<double> diff = timer_end - timer_start;
    std::chrono::duration<double> sort_time = sort_end - sort_start;
    std::cerr << "Time to get indexing result for " << count_db << " fingerprints in DB : " << diff.count()  << std::endl;
    std::cerr << "Time for sorting " << sort_time.count() << std::endl;

    cudaFree(d_fp);
    cudaFree(d_db);
    cudaFree(d_result);
    cudaFree(d_mapping);
    cudaFree(d_s1_result);
    cudaFree(d_s2_result);
    cudaFree(d_s3_result);
    cudaFree(d_s4_result);
    cudaFree(d_ids);
    cudaFree(d_s_sum);
    cudaFree(d_cos_sum);
    cudaFree(d_sin_sum);

    return 0;
}

// nvcc -o parallel_indexing_split_s1 parallel_indexing_split_s1.cu fingerprint_structure.cpp -std=c++11 -lineinfo
