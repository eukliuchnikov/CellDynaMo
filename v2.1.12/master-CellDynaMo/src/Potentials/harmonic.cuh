#pragma once
#include "langevin.cuh"

__global__ void computeHarmonic(float3* r, float3* f, int N, int* type, int kn, float rk, int* harmonicKinCount, int* harmonicKin, float* harmonicKinRadii, float* kt_radius, Param* d_parameters, float* length, int* connect1, int* connect2, float* len1, float* len2, float* f_out);
