#pragma once
#include "../Initializiation/visualization.h"

#include "harmonic.cuh"
#include "bending.cuh"		
#include "membrane.h"		
#include "pushing.cuh"
#include "pulling.cuh"
#include "kt_ex_vol.cuh"
#include "cylinder.cuh"		

__global__ void integrateGPU(float3* r, float3* f, int N, int* type, int kn, float3 l_f, float3 r_f, Param* d_parameters);
float3 kin_force(int index, float3* f, Param* d_parameters);
void KT_bonds_init(float* d_kin_cos, float* kin_cos, int* d_harmonicKinCount, int* a_harmonicKinCount, float* d_harmonicKinRadii, float* a_harmonicKinRadii, int* d_harmonicKin, int* a_harmonicKin, Param* d_parameters);
void cudaRANDinit();
void dynamics(DCD dcd, float time);
__global__ void CoM(float3* r, float3* f, int N, int* type, int kn);
