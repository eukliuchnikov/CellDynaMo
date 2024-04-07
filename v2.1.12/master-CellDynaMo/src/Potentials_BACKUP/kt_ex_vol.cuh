#pragma once
#include "langevin.cuh"

__global__ void excl_vol_kt(float3* r, float3* f, int N, int* type, int kn, Param* d_parameters);
