#pragma once
#include "langevin.cuh"

__global__ void excl_vol_kt(float3* r, float3* f, int N, int* type, int kn, Param* d_parameters);
__global__ void chrom_ter(float3* r, float3* f, int N, int* type, int kn, Param* d_parameters, float radi);
