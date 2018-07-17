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

__global__ void calculate_s1(fingerprint* db, fingerprint* fp, float* s1, float* result) {
    __shared__ float ss, scos, ssin, s_addition, s_absdiff;
    int j = blockIdx.x;
    int i = threadIdx.x;
    // int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (i == 0) {
        ss = 0;
        scos = 0;
        ssin = 0;
        s_addition = 0.0f;
        s_absdiff = 0.0f;
    }
    __syncthreads();
    // (db+j)->local_frequency[i] = dfrequency_t_byte(dbyte_to_frequency((db+j)->local_frequency[i])+0.1);
    float s = dbyte_to_coherence(fp->local_coherence[i])*dbyte_to_coherence((db+j)->local_coherence[i]);
    float d = M_PI/180.0f * 2 * (dbyte_to_orientation(fp->local_orientation[i])-dbyte_to_orientation((db+j)->local_orientation[i]));
    float tcos = s*cos(d);
    float tsin = s*sin(d);
    float t_addition = dbyte_to_frequency(fp->local_frequency[i]) + dbyte_to_frequency((db+j)->local_frequency[i]);
    float t_absdiff = abs(dbyte_to_frequency(fp->local_frequency[i]) - dbyte_to_frequency((db+j)->local_frequency[i]));

    atomicAdd(&ss, s);
    atomicAdd(&scos, tcos);
    atomicAdd(&ssin, tsin);
    atomicAdd(&s_addition, t_addition);
    atomicAdd(&s_absdiff, t_absdiff);
    __syncthreads();
    if (i == 0) {
        float s2, s3, s4;
        if (ss != 0) s1[j] = sqrt(pow(scos,2)+pow(ssin,2))/ss;
        else s1[j] = 0;
        s2 = 1 - (s_absdiff/s_addition);
        s3 = 1 - (abs(dbyte_to_frequency(fp->avg_frequency)-dbyte_to_frequency((db+j)->avg_frequency))/max(dbyte_to_frequency(fp->avg_frequency), dbyte_to_frequency((db+j)->avg_frequency)));
        s4 = 1-(abs(dbyte_to_orientation(fp->avg_orientation)-dbyte_to_orientation((db+j)->avg_orientation))/180.0f);
        result[j] = w1*s1[j] + w2*s2 + w3*s3 + w4*s4;
    }
}

__global__ void get_best_core_s1(fingerprint* db, float* s1, float* result, int* mapping) {
    int i = blockIdx.x;
    if ((db+i)->id%5 == 1) {
        int max_idx = i;
        for (int j=1 ; j<5 ; j++) {
            if ((db+i+j)->id%5 == 1) break;
            else {
                if (s1[i+j] > s1[max_idx]) {
                    max_idx = i+j;
                }
            }
        }
        mapping[((db+i)->id-1)/5] = max_idx;
        float best_result = result[max_idx];
        __syncthreads();
        result[((db+i)->id-1)/5] = best_result;
    }
}

int main(int argc, char** argv) {
    if (argc < 3) {
        // std::cerr << "Usage : ./parallel_indexing fingerprint-to-be-searched fingerprint-db\n";
        return 0;
    }

    std::string fp_filename = argv[1];
    std::string db_filename = argv[2];
    // std::cerr << "FP " << fp_filename << " DB " << db_filename << std::endl;

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
    // float s1_result[count_db], s2_result[count_db], s3_result[count_db], s4_result[count_db];
    // float result[count_db];
    std::vector<float> result(count_db_fingerprint, 0);
    float *d_result;
    float *d_s1_result;
    // std::cout << "Starting cudaMalloc\n";
    cudaMalloc((void **)&d_s1_result, count_db*sizeof(float));
    
    cudaMalloc((void **)&d_fp, sizeof(fingerprint));
    cudaMalloc((void **)&d_db, count_db*sizeof(fingerprint));
    cudaMalloc((void **)&d_result, count_db*sizeof(float));

    //Mapping for block idx to fingerprint core idx
    int *d_mapping;
    cudaMalloc((void **)&d_mapping, count_db_fingerprint*sizeof(int));

    // std::cout << "Starting cudaMemcpy\n";

    cudaMemcpy(d_db, &db[0], count_db*sizeof(fingerprint), cudaMemcpyHostToDevice);
    cudaMemcpy(d_fp, &fp[0], sizeof(fingerprint), cudaMemcpyHostToDevice);
    // auto timer_start = std::chrono::steady_clock::now();
    // std::cout << "Starting S1\n";
    calculate_s1<<<count_db,BLOCKSIZE>>>(d_db, d_fp, d_s1_result, d_result);
    get_best_core_s1<<<count_db, 1>>>(d_db, d_s1_result, d_result, d_mapping);
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );
    // cudaMemcpy(&s1_result[0], d_result, count_db*sizeof(float), cudaMemcpyDeviceToHost);
    

    // int mapping[count_db_fingerprint];
    // memset(mapping, 0, sizeof(mapping));
    std::vector<int> mapping(count_db_fingerprint, 0);
    cudaMemcpy(&mapping[0], d_mapping, count_db_fingerprint*sizeof(int), cudaMemcpyDeviceToHost);
    std::vector<float> s1_result;
    s1_result.resize(count_db, 0);
    cudaMemcpy(&s1_result[0], d_s1_result, count_db*sizeof(float), cudaMemcpyDeviceToHost);

    // float* t_result = new float[count_db_fingerprint];
    // std::vector<float> t_result(count_db_fingerprint, 0);
    int* ids = new int[count_db_fingerprint];
    for (int i=0 ; i<count_db_fingerprint ; i++) {
        ids[i] = db[mapping[i]].id;
    }
    int* d_ids;
    cudaMalloc((void **)&d_ids, count_db_fingerprint*sizeof(int));
    cudaMemcpy(d_ids, &ids[0], count_db_fingerprint*sizeof(int), cudaMemcpyHostToDevice);
    // std::vector<int> ids(count_db_fingerprint, 0);
    // std::vector< std::pair<float, int> > best_matches;
    // for (int i=0 ; i<count_db_fingerprint ; i++) {
    //     // std::cout << "result = " << result[i] << std::endl;
        // best_matches.push_back(std::make_pair(result[i], db[mapping[i]].id));
        // t_result[i] = result[i];
        // ids[i] = db[mapping[i]].id;
    // }
    // sort(best_matches.rbegin(), best_matches.rend());
    // thrust::sort(best_matches.rbegin(), best_matches.rend());
    // thrust::sort(result, result+count_db_fingerprint);
    // // std::cout << "Before sort\n";
    // for (int i=0 ; i<count_db_fingerprint ; i++) {
    //     // std::cout << i << " " << ids[i] << " " << t_result[i] << std::endl;
    // }
    auto sort_start = std::chrono::steady_clock::now();
    // thrust::sort_by_key(t_result, t_result+count_db_fingerprint, ids);
    thrust::sort_by_key(thrust::device, d_result, d_result+count_db_fingerprint, d_ids);
    auto sort_end = std::chrono::steady_clock::now();
    cudaMemcpy(&result[0], d_result, count_db_fingerprint*sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(&ids[0], d_ids, count_db_fingerprint*sizeof(int), cudaMemcpyDeviceToHost);
    // std::cout << "\nBest match\n";
    /*for (int i=0 ; i<best_matches.size() ; i++) {
        // // std::cout << "ID " << best_matches[i].second << "-"<< best_matches[i].second/5 <<"\t: " << best_matches[i].first;
        // std::cout << std::endl;
    }*/
    for (int i=count_db_fingerprint-1 ; i>=0 ; i--) {
        std::cout << "ID " << ids[i] << "-"<< ids[i]/5 <<"\t: " << result[i];
        std::cout << std::endl;
    }
    auto timer_end = std::chrono::steady_clock::now();
    std::chrono::duration<double> diff = timer_end - timer_start;
    std::chrono::duration<double> sort_time = sort_end - sort_start;
    std::cerr << "Time to get indexing result for " << count_db << " fingerprints in DB : " << diff.count()  << std::endl;
    std::cerr << "Time for sorting " << sort_time.count() << std::endl;

/*
    // DEBUG
    // // std::cout << "\nS1\n";
    // for (int i=0 ; i<count_db ; i++) {
    //     // std::cout << s1_result[i] << std::endl;
    // }
    std::vector<float> s2_result, s3_result, s4_result;
    s1_result.resize(count_db, 0);
    s2_result.resize(count_db_fingerprint, 0);
    s3_result.resize(count_db_fingerprint, 0);
    s4_result.resize(count_db_fingerprint, 0);
    cudaMemcpy(&s1_result[0], d_s1_result, count_db*sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(&s2_result[0], d_s2_result, count_db_fingerprint*sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(&s3_result[0], d_s3_result, count_db_fingerprint*sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(&s4_result[0], d_s4_result, count_db_fingerprint*sizeof(float), cudaMemcpyDeviceToHost);

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
*/
    cudaFree(d_fp);
    cudaFree(d_db);
    cudaFree(d_result);
    cudaFree(d_mapping);
    cudaFree(d_s1_result);
    // cudaFree(d_final_result);

    return 0;
}

// nvcc -o parallel_indexing_op1 parallel_indexing_op1.cu fingerprint_structure.cpp -std=c++11
