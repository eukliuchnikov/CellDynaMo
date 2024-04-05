#pragma once
#include "langevin.cuh"

__global__ void pushing(float3* r, float3* f, int N, int* type, int kn, long long int step, int* n_force, float* f_force, Param* d_parameters);
