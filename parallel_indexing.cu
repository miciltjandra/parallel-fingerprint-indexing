#include <stdio.h>
#include <iostream>
#include <algorithm>
#include <chrono>
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
        if (ss != 0) result[j] = sqrt(pow(scos,2)+pow(ssin,2))/ss;
        else result[j] = 0;
    }
    /*__syncthreads();
    if (i == 0) printf("Idx %d %f\n", j, result[j]);
    // First core of a fingerprint check maximum from all core
    if (i == 0 && (db+j)->id%5 == 1) {
        // printf("Index %d\n", j);
        int max_idx = j;
        for (int k=1 ; k<5 ; k++) {
            if ((db+j+k)->id%5 == 1) break;
            else {
                // printf("%d %d %f %f\n", max_idx, j+k, result[max_idx], result[j+k]);
                if (result[j+k] > result[max_idx]) {
                    // printf("%d %d\n", max_idx, j+k);
                    max_idx = j+k;
                }
            }
        }
        mapping[((db+j)->id-1)/5] = max_idx;
    }*/
}

__global__ void get_best_core_s1(fingerprint* db, float* result, int* mapping) {
    int i = blockIdx.x;
    // printf("Index %d %d\n", i, (db+i)->id);
    if ((db+i)->id%5 == 1) {
        int max_idx = i;
        for (int j=1 ; j<5 ; j++) {
            if ((db+i+j)->id%5 == 1) break;
            else {
                // printf("%d %d %f %f\n", max_idx, i+j, result[max_idx], result[i+j]);
                if (result[i+j] > result[max_idx]) {
                    // printf("%d %d\n", max_idx, i+j);
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
    // printf("S2 idx %d %d\n", blockIdx.x, j);
    // int idx = blockIdx.x*blockDim.x + threadIdx.x;
    float t_addition = dbyte_to_frequency(fp->local_frequency[i]) + dbyte_to_frequency((db+j)->local_frequency[i]);
    float t_absdiff = abs(dbyte_to_frequency(fp->local_frequency[i]) - dbyte_to_frequency((db+j)->local_frequency[i]));
    if (j == 16) {
        // printf("%d %f\n", i, t_addition);
        // printf("%d %f %f\n", i, dbyte_to_frequency(fp->local_frequency[i]), dbyte_to_frequency((db+j)->local_frequency[i]));
    }
    atomicAdd(&s_addition, t_addition);
    atomicAdd(&s_absdiff, t_absdiff);
    __syncthreads();
    if (i == 0) {
        // printf("%d %f %f\n", j, s_addition, s_absdiff);
        result[blockIdx.x] = 1 - (s_absdiff/s_addition);
        // printf("%d %f %f %f\n", blockIdx.x, s_absdiff, s_addition, result[blockIdx.x]);
    }
}

__global__ void calculate_s3(fingerprint* db, fingerprint* fp, float* result, int* mapping) {
    int j = mapping[blockIdx.x];
    result[blockIdx.x] = 1 - (abs(dbyte_to_frequency(fp->avg_frequency)-dbyte_to_frequency((db+j)->avg_frequency))/max(dbyte_to_frequency(fp->avg_frequency), dbyte_to_frequency((db+j)->avg_frequency)));
}

__global__ void calculate_s4(fingerprint* db, fingerprint* fp, float* result, int* mapping) {
    int j = mapping[blockIdx.x];
    result[blockIdx.x] = 1-(abs(dbyte_to_orientation(fp->avg_orientation)-dbyte_to_orientation((db+j)->avg_orientation))/M_PI);
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
    std::cerr << "FP " << fp_filename << " DB " << db_filename << std::endl;

    // Read the fingerprint to be searched
    std::vector<struct fingerprint> fp;
    int count_fp = read_from_file(fp, fp_filename);

    std::vector<float> local_orie, local_cohe, local_freq;
    get_fingerprint_local_values(fp[0], local_orie, local_cohe, local_freq);
    float avg_o = get_fingerprint_average_orientation(fp[0]);
    float avg_f = get_fingerprint_average_frequency(fp[0]);

    // Read the database
    std::vector<struct fingerprint> db;
    int count_db = read_from_file(db, db_filename);
    std::cerr << "Fingerprint core database count : " << count_db << std::endl;

    std::cerr << "Last fingerprint ID : " << db[count_db-1].id << std::endl;
    int count_db_fingerprint = (db[count_db-1].id-1)/5+1;
    std::cerr << "Fingerprint database count : " << count_db_fingerprint << std::endl;

    auto timer_start = std::chrono::steady_clock::now();
    // S1
    fingerprint *d_fp, *d_db;
    float s1_result[count_db], s2_result[count_db], s3_result[count_db], s4_result[count_db];
    float result[count_db];
    float *d_result;

    cudaMalloc((void **)&d_fp, sizeof(fingerprint));
    cudaMalloc((void **)&d_db, count_db*sizeof(fingerprint));
    cudaMalloc((void **)&d_result, count_db*sizeof(float));

    //Mapping for block idx to fingerprint core idx
    int *d_mapping;
    cudaMalloc((void **)&d_mapping, count_db_fingerprint*sizeof(int));

    cudaMemcpy(d_db, &db[0], count_db*sizeof(fingerprint), cudaMemcpyHostToDevice);
    cudaMemcpy(d_fp, &fp[0], sizeof(fingerprint), cudaMemcpyHostToDevice);
    calculate_s1<<<count_db,BLOCKSIZE>>>(d_db, d_fp, d_result);
    get_best_core_s1<<<count_db, 1>>>(d_db, d_result, d_mapping);
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );
    cudaMemcpy(&s1_result[0], d_result, count_db*sizeof(float), cudaMemcpyDeviceToHost);
    

    int mapping[count_db_fingerprint];
    memset(mapping, 0, sizeof(mapping));
    cudaMemcpy(&mapping[0], d_mapping, count_db_fingerprint*sizeof(int), cudaMemcpyDeviceToHost);

    int counter = 0;
    for (int i=0 ; i<count_db ; i++) {
        std::cout << i << " : ID " << db[i].id << " ";
        std::cout << "S1 = " << s1_result[i] << std::endl;
        if (i==count_db-1 || (db[i+1].id%5 == 1)) {
            std::cout << counter << " " << mapping[counter] << std::endl;
            counter++;
        }
    }

    // Check mapping
    std::cout << "BEST CORE\n";
    for (int i=0 ; i<count_db_fingerprint ; i++) {
        std::cout << i << " " << mapping[i] << std::endl;
    }

    // Test S2
    // Only calculate for 1 core per fingerprint
    // calculate_s2<<<count_db,BLOCKSIZE>>>(d_db, d_fp, d_result);
    calculate_s2<<<count_db_fingerprint,BLOCKSIZE>>>(d_db, d_fp, d_result, d_mapping);
    cudaMemcpy(&s2_result[0], d_result, count_db_fingerprint*sizeof(float), cudaMemcpyDeviceToHost);
    // std::cout << "\n\nS2\n";
    // for (int i=0 ; i<count_db ; i++) {
    //     std::cout << i << " : ID " << db[i].id << std::endl;
    //     std::cout << "result = " << s2_result[i] << std::endl;
    // }

    // Test S3
    // calculate_s3<<<count_db,1>>>(d_db, d_fp, d_result);
    calculate_s3<<<count_db_fingerprint,1>>>(d_db, d_fp, d_result,d_mapping);
    cudaMemcpy(&s3_result[0], d_result, count_db*sizeof(float), cudaMemcpyDeviceToHost);
    // std::cout << "\n\nS3\n";
    // for (int i=0 ; i<count_db ; i++) {
    //     std::cout << i << " : ID " << db[i].id << std::endl;
    //     std::cout << "result = " << s3_result[i] << std::endl;
    // }

    // Test S4
    // calculate_s4<<<count_db,1>>>(d_db, d_fp, d_result);
    calculate_s4<<<count_db_fingerprint,1>>>(d_db, d_fp, d_result, d_mapping);
    cudaMemcpy(&s4_result[0], d_result, count_db*sizeof(float), cudaMemcpyDeviceToHost);
    // std::cout << "\n\nS4\n";
    // for (int i=0 ; i<count_db ; i++) {
    //     std::cout << i << " : ID " << db[i].id << std::endl;
    //     std::cout << "result = " << s4_result[i] << std::endl;
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
    // calculate_s<<<1,count_db>>>(d_s1_result, d_s2_result, d_s3_result, d_s4_result, d_result);
    calculate_s<<<count_db_fingerprint, 1>>>(d_s1_result, d_s2_result, d_s3_result, d_s4_result, d_result, d_mapping);
    cudaMemcpy(&result[0], d_result, count_db*sizeof(float), cudaMemcpyDeviceToHost);
    // S
    // for (int i=0 ; i<count_db ; i++) {
    /*for (int i=0 ; i<count_db_fingerprint ; i++) {
        std::cout << i << " : ID " << db[i].id << std::endl;
        std::cout << "result = " << result[i] << std::endl;
    }*/

    // float *d_final_result;
    /* This is for when not used with mapping */
    // get_top_fingerprints<<<1,count_db_fingerprint>>>(d_result, d_final_result, d_mapping);
    // cudaMemcpy(&result[0], d_final_result, count_db_fingerprint*sizeof(float), cudaMemcpyDeviceToHost);
    // std::cout << "\n\nFinal Result\n";
    std::vector< std::pair<float, int> > best_matches;
    for (int i=0 ; i<count_db_fingerprint ; i++) {
    //     std::cout << "result = " << result[i] << std::endl;
        best_matches.push_back(std::make_pair(result[i], db[mapping[i]].id));
    }
    sort(best_matches.rbegin(), best_matches.rend());
    std::cout << "\nBest match\n";
    for (int i=0 ; i<best_matches.size() ; i++) {
        std::cout << "ID " << best_matches[i].second << "-"<< best_matches[i].second/5 <<"\t: " << best_matches[i].first;
        std::cout << std::endl;
    }
    auto timer_end = std::chrono::steady_clock::now();
    std::chrono::duration<double> diff = timer_end - timer_start;
    std::cout << "Time to get indexing result for " << count_db << " fingerprints in DB : " << diff.count()  << std::endl;

    // DEBUG
    // std::cout << "\nS1\n";
    // for (int i=0 ; i<count_db ; i++) {
    //     std::cout << s1_result[i] << std::endl;
    // }

    std::cout << "\nS1\n";
    for (int i=0 ; i<count_db_fingerprint ; i++) {
        std::cout << s1_result[mapping[i]] << std::endl;
    }

    std::cout << "\nS2\n";
    for (int i=0 ; i<count_db_fingerprint ; i++) {
        std::cout << s2_result[i] << std::endl;
    }
    
    std::cout << "\nS3\n";
    for (int i=0 ; i<count_db_fingerprint ; i++) {
        std::cout << s3_result[i] << std::endl;
    }

    std::cout << "\nS4\n";
    for (int i=0 ; i<count_db_fingerprint ; i++) {
        std::cout << s4_result[i] << std::endl;
    }

    std::cout << "\nS\n";
    for (int i=0 ; i<count_db_fingerprint ; i++) {
        std::cout << result[i] << std::endl;
    }

    cudaFree(d_fp);
    cudaFree(d_db);
    cudaFree(d_result);
    cudaFree(d_mapping);
    cudaFree(d_s1_result);
    cudaFree(d_s2_result);
    cudaFree(d_s3_result);
    cudaFree(d_s4_result);
    // cudaFree(d_final_result);

    return 0;
}

// nvcc -o parallel_indexing parallel_indexing.cu fingerprint_structure.cpp -std=c++11
