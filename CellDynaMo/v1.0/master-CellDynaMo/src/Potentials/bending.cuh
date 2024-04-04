#pragma once
#include "langevin.cuh"

__global__ void computeAngles(float3* r, float3* f, int N, int* type, int kn, float* kt_cos, float* dyn_cos, Param* d_parameters, float* length);
