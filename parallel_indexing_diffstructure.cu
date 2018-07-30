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

__global__ void calculate_s1_preparation(fingerprint* fp, unsigned char* orientations, unsigned char* coherences, float* s_sum, float* cos_sum, float* sin_sum, int count) {
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    int local_i = i%36;
    // printf("%d\n", i);
    if (i >= count*36) return;

    float s = dbyte_to_coherence(fp->local_coherence[local_i])*dbyte_to_coherence(coherences[i]);
    float d = M_PI/180.0f * 2 * (dbyte_to_orientation(fp->local_orientation[local_i])-dbyte_to_orientation(orientations[i]));
    float tcos = s*cos(d);
    float tsin = s*sin(d);

    // printf("%d %d %f %f\n", i, local_i, dbyte_to_coherence(fp->local_coherence[local_i]), dbyte_to_coherence(coherences[i]));

    atomicAdd(&s_sum[i/36], s);
    atomicAdd(&cos_sum[i/36], tcos);
    atomicAdd(&sin_sum[i/36], tsin);
    __syncthreads();
}

__global__ void calculate_s1(float* s_sum, float* cos_sum, float* sin_sum, float* result, int count) {
    int j = blockIdx.x*blockDim.x + threadIdx.x;
    if (j<count) {
        result[j] = sqrt(pow(cos_sum[j],2)+pow(sin_sum[j],2))/s_sum[j];
        // printf("S1 %d %f\n", j, result[j]);
    }
}


__global__ void get_best_core_s1(int* core_ids, float* result, int* mapping, int count) {
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= count) return;
    // printf("%d\n", i);
    if (core_ids[i]%5 == 1) {
        int max_idx = i;
        for (int j=1 ; j<5 ; j++) {
            if (i+j == count || core_ids[i+j]%5 == 1) break;
            else {
                if (result[i+j] > result[max_idx]) {
                    max_idx = i+j;
                }
            }
        }
        // printf("%d %d\n", (core_ids[i]-1)/5, max_idx);
        mapping[(core_ids[i]-1)/5] = max_idx;
    }
}

__global__ void calculate_s2_preparation(fingerprint* fp, unsigned char* frequencies, float* s_addition, float* s_absdiff, int count) {
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    int local_i = i%36;
    if (i >= count*36) return;

    float t_addition = dbyte_to_frequency(fp->local_frequency[local_i]) + dbyte_to_frequency(frequencies[i]);
    float t_absdiff = abs(dbyte_to_frequency(fp->local_frequency[local_i]) - dbyte_to_frequency(frequencies[i]));
    atomicAdd(&s_addition[i/36], t_addition);
    atomicAdd(&s_absdiff[i/36], t_absdiff);
    __syncthreads();
}

__global__ void calculate_s2(float* s_addition, float* s_absdiff, float* result, int count) {
    int j = blockIdx.x*blockDim.x + threadIdx.x;
    if (j<count) {
        result[j] = 1 - (s_absdiff[j]/s_addition[j]);
        // printf("S2 %d %f\n", j, result[j]);
    }
}


__global__ void calculate_s3(fingerprint* fp, unsigned char* global_frequencies, float* result, int count) {
    int j = blockIdx.x*blockDim.x + threadIdx.x;
    if (j >= count) return;
    result[j] = 1 - (abs(dbyte_to_frequency(fp->avg_frequency)-dbyte_to_frequency(global_frequencies[j]))/max(dbyte_to_frequency(fp->avg_frequency), dbyte_to_frequency(global_frequencies[j])));
    // printf("S3 %d %f\n", j, result[j]);
}

__global__ void calculate_s4(fingerprint* fp, unsigned char* global_orientations, float* result, int count) {
    int j = blockIdx.x*blockDim.x + threadIdx.x;
    if (j >= count) return;
    // printf("S4 %d %f %f\n", j, dbyte_to_orientation(fp->avg_orientation), dbyte_to_orientation(global_orientations[j]));
    result[j] = 1-(abs(dbyte_to_orientation(fp->avg_orientation)-dbyte_to_orientation(global_orientations[j]))/180.0f);
    // printf("S4 %d %f\n", j, result[j]);
}

__global__ void calculate_s(float* s1, float* s2, float*s3, float* s4, float* result, int* mapping) {
    int i = mapping[blockIdx.x];
    // printf("S %d\n", i);
    result[blockIdx.x] = w1*s1[i] + w2*s2[i] + w3*s3[i] + w4*s4[i];
    // printf("S %d %f\n", i, result[blockIdx.x]);
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
    int *core_ids;
    unsigned char *local_orientations, *local_coherences, *local_frequencies, *global_orientations, *global_frequencies;
    int count_db = read_translated_structure(db_filename, core_ids, local_orientations, local_coherences, local_frequencies, global_orientations, global_frequencies);
    std::cerr << "Fingerprint core database count : " << count_db << std::endl;

    std::cerr << "Last fingerprint ID : " << core_ids[count_db-1] << std::endl;
    int count_db_fingerprint = (core_ids[count_db-1]-1)/5+1;
    std::cerr << "Fingerprint database count : " << count_db_fingerprint << std::endl;
    /*for (int i=0 ; i<count_db ; i++) {
        for (int j=0 ; j<36 ; j++) {
            std::cout << byte_to_orientation(local_orientations[i*36+j]) << " ";
        }
        std::cout << std::endl;
    }

    for (int i=0 ; i<count_db ; i++) {
        for (int j=0 ; j<36 ; j++) {
            std::cout << byte_to_coherence(local_coherences[i*36+j]) << " ";
        }
        std::cout << std::endl;
    }

    for (int i=0 ; i<count_db ; i++) {
        for (int j=0 ; j<36 ; j++) {
            std::cout << byte_to_frequency(local_frequencies[i*36+j]) << " ";
        }
        std::cout << std::endl;
    }

    for (int i=0 ; i<count_db ; i++) {
        std::cout << byte_to_orientation(global_orientations[i]) << std::endl;
    }

    for (int i=0 ; i<count_db ; i++) {
        std::cout << byte_to_frequency(global_frequencies[i]) << std::endl;
    }*/

    auto timer_start = std::chrono::steady_clock::now();

    // Preparing memory
    fingerprint *d_fp;
    std::vector<float> result(count_db_fingerprint, 0);
    float *d_s1_result, *d_s2_result, *d_s3_result, *d_s4_result, *d_result;

    int *d_core_ids;
    unsigned char *d_local_orientations, *d_local_coherences, *d_local_frequencies, *d_global_orientations, *d_global_frequencies;
    
    cudaMalloc((void **)&d_fp, sizeof(fingerprint));

    cudaMalloc((void **)&d_core_ids, count_db*sizeof(int));
    cudaMalloc((void **)&d_local_orientations, count_db*36*sizeof(unsigned char));
    cudaMalloc((void **)&d_local_coherences, count_db*36*sizeof(unsigned char));
    cudaMalloc((void **)&d_local_frequencies, count_db*36*sizeof(unsigned char));
    cudaMalloc((void **)&d_global_orientations, count_db*sizeof(unsigned char));
    cudaMalloc((void **)&d_global_frequencies, count_db*sizeof(unsigned char));

    cudaMalloc((void **)&d_s1_result, count_db*sizeof(float));
    cudaMalloc((void **)&d_s2_result, count_db*sizeof(float));
    cudaMalloc((void **)&d_s3_result, count_db*sizeof(float));
    cudaMalloc((void **)&d_s4_result, count_db*sizeof(float));
    cudaMalloc((void **)&d_result, count_db_fingerprint*sizeof(float));

    printf("Malloc done\n");

    //Mapping for fingerprint to fingerprint core idx
    int *d_mapping;
    cudaMalloc((void **)&d_mapping, count_db_fingerprint*sizeof(int));

    cudaMemcpy(d_fp, &fp[0], sizeof(fingerprint), cudaMemcpyHostToDevice);

    printf("Memcpy FP done\n");

    cudaMemcpy(d_core_ids, core_ids, count_db*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_local_orientations, local_orientations, count_db*36*sizeof(unsigned char), cudaMemcpyHostToDevice);
    cudaMemcpy(d_local_coherences, local_coherences, count_db*36*sizeof(unsigned char), cudaMemcpyHostToDevice);
    cudaMemcpy(d_local_frequencies, local_frequencies, count_db*36*sizeof(unsigned char), cudaMemcpyHostToDevice);
    cudaMemcpy(d_global_orientations, global_orientations, count_db*sizeof(unsigned char), cudaMemcpyHostToDevice);
    cudaMemcpy(d_global_frequencies, global_frequencies, count_db*sizeof(unsigned char), cudaMemcpyHostToDevice);

    printf("Memcpy done\n");

    //Additional Memory for S1
    float *d_s_sum, *d_cos_sum, *d_sin_sum;
    cudaMalloc((void **)&d_s_sum, count_db*sizeof(float));
    cudaMalloc((void **)&d_cos_sum, count_db*sizeof(float));
    cudaMalloc((void **)&d_sin_sum, count_db*sizeof(float));

    //Additional Memory for S2
    float *d_s_addition, *d_s_absdiff;
    cudaMalloc((void **)&d_s_addition, count_db*sizeof(float));
    cudaMalloc((void **)&d_s_absdiff, count_db*sizeof(float));

    printf("Malloc again done\n");

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );

    // Use streams for S1-S4
    cudaStream_t streams[4];
    cudaStreamCreate(&streams[0]);
    cudaStreamCreate(&streams[1]);
    cudaStreamCreate(&streams[2]);
    cudaStreamCreate(&streams[3]);

    // S1
    std::cerr << "Num of block : " << ((count_db*36)/256)+1 << std::endl;
    calculate_s1_preparation<<<((count_db*36)/256)+1, 256, 32, streams[0]>>>(d_fp, d_local_orientations, d_local_coherences, d_s_sum, d_cos_sum, d_sin_sum, count_db);
    std::cerr << "s1 prep done\n";
    calculate_s1<<<(count_db/256)+1, 256, 4, streams[0]>>>(d_s_sum, d_cos_sum, d_sin_sum, d_s1_result, count_db);
    std::cerr << "s1 done\n";
    get_best_core_s1<<<(count_db/256)+1, 256, 32, streams[0]>>>(d_core_ids, d_s1_result, d_mapping, count_db);
    std::cerr << "best core done\n";

    // std::vector<float> s1_result;
    // s1_result.resize(count_db, 0);
    // cudaMemcpy(&s1_result[0], d_s1_result, count_db*sizeof(float), cudaMemcpyDeviceToHost);

    // S2
    // Only calculate for 1 core per fingerprint using mapping
    calculate_s2_preparation<<<((count_db*36)/256)+1, 256, 32, streams[1]>>>(d_fp, d_local_frequencies, d_s_addition, d_s_absdiff, count_db);
    calculate_s2<<<(count_db/256)+1, 256, 4, streams[1]>>>(d_s_addition, d_s_absdiff, d_s2_result, count_db);
    // cudaMemcpy(&s2_result[0], d_s2_result, count_db_fingerprint*sizeof(float), cudaMemcpyDeviceToHost);

    // S3
    calculate_s3<<<(count_db/256)+1, 256, 4, streams[2]>>>(d_fp, d_global_frequencies, d_s3_result, count_db);
    // cudaMemcpy(&s3_result[0], d_s3_result, count_db*sizeof(float), cudaMemcpyDeviceToHost);

    // S4
    calculate_s4<<<(count_db/256)+1, 256, 4, streams[3]>>>(d_fp, d_global_orientations, d_s4_result, count_db);
    // cudaMemcpy(&s4_result[0], d_s4_result, count_db*sizeof(float), cudaMemcpyDeviceToHost);

    std::vector<int> mapping(count_db_fingerprint, 0);
    cudaMemcpy(&mapping[0], d_mapping, count_db_fingerprint*sizeof(int), cudaMemcpyDeviceToHost);

    cudaDeviceSynchronize();
    // S
    calculate_s<<<count_db_fingerprint, 1>>>(d_s1_result, d_s2_result, d_s3_result, d_s4_result, d_result, d_mapping);
    // cudaMemcpy(&result[0], d_result, count_db_fingerprint*sizeof(float), cudaMemcpyDeviceToHost);

    // ID for identifying fingerprint during sort
    int* ids = new int[count_db_fingerprint];
    for (int i=0 ; i<count_db_fingerprint ; i++) {
        ids[i] = core_ids[mapping[i]];
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
    cudaFree(d_result);
    cudaFree(d_mapping);
    cudaFree(d_s1_result);
    cudaFree(d_s2_result);
    cudaFree(d_s3_result);
    cudaFree(d_s4_result);
    cudaFree(d_ids);

    delete[] ids;
    delete[] local_orientations;
    delete[] local_coherences;
    delete[] local_frequencies;
    delete[] global_orientations;
    delete[] global_frequencies;


    return 0;
}

// nvcc -o parallel_indexing_diffstructure parallel_indexing_diffstructure.cu fingerprint_structure.cpp -std=c++11 -lineinfo
