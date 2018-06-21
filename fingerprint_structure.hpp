#ifndef FINGERPRINT_STRUCTURE_H
#define FINGERPRINT_STRUCTURE_H

#include <stdio.h>
#include <vector>
#include <string>
#include <cstring>

// static const float orientation_unit = 0.703125;
static const float orientation_unit = 180.0f/256.0f;
static const float coherence_unit = 255.0f;
static const float period_unit = 0.1f;

struct fingerprint {
    int id;
    unsigned char local_orientation[36];
    unsigned char local_coherence[36];
    unsigned char local_frequency[36];
    unsigned char avg_orientation;
    unsigned char avg_frequency; 
};

__host__ __device__ unsigned char orientation_to_byte(float orientation);

__host__ __device__ float byte_to_orientation(unsigned char c);

__host__ __device__ unsigned char coherence_to_byte(float coherence);

__host__ __device__ float byte_to_coherence(unsigned char c);

__host__ __device__ unsigned char frequency_to_byte(float frequency);

__host__ __device__ float byte_to_frequency(unsigned char c);

__host__ __device__ unsigned char period_to_byte(float period);

__host__ __device__ float byte_to_period(unsigned char c);

__host__ __device__ struct fingerprint make_fingerprint_struct(int id, std::vector<float> local_orientation, std::vector<float> local_coherence, std::vector<float> local_frequency, float avg_orie, float avg_freq);

__host__ __device__ void print_fingerprint_struct(struct fingerprint fp);

__host__ __device__ void get_fingerprint_local_values(struct fingerprint fp, std::vector<float> &local_orientation, std::vector<float> &local_coherence, std::vector<float> &local_frequency);

__host__ __device__ float get_fingerprint_average_orientation(struct fingerprint fp);

__host__ __device__ float get_fingerprint_average_frequency(struct fingerprint fp);

__host__ __device__ void save_to_file(int size, struct fingerprint fps[], std::string filename);

__host__ __device__ int read_from_file(std::vector<struct fingerprint> &fps, std::string filename);

__host__ __device__ int get_last_id_from_file(std::string filename);

__host__ __device__ int get_new_fingerprint_id(int last_id);

#endif