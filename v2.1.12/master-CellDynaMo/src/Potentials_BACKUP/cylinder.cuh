#pragma once
#include "langevin.cuh"

__global__ void cyl_cyl(float3* r, float3* f, int N, int* type, Param* d_parameters);
__global__ void cyl_sphere(float3* r, float3* f, int N, int* type, int kn, Param* d_parameters);
//__global__ void end_force(float3* r, float3* f, int N);
