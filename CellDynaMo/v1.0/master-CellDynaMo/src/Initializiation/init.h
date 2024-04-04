#pragma once

#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>

#include "parameters.h"
#include "md.h"
#include "../Math/ran2.h"
#include "../Math/mat.h"
#include "visualization.h"
#include "mt_init.h"
#include "kt_init.h"
#include <cuda_runtime_api.h>
#include <cuda.h>

extern int seed;
//COMPONENTS
extern MT* pol[POLE_COUNT];
extern KT kt[KIN_COUNT];
extern float rk;
extern int kn;
extern float* kin_cos;
extern float* dyn_cos;
extern float* kin_rad;
extern int* harmonicKinCount[KIN_COUNT];
extern int* harmonicKin[KIN_COUNT];
extern float* harmonicKinRadii[KIN_COUNT];
extern int* a_harmonicKinCount;
extern int* a_harmonicKin;
extern float* a_harmonicKinRadii;

extern float* d_kin_cos;
extern float* d_dyn_cos;
extern float* d_kin_rad;
extern int* d_harmonicKinCount;
extern float* d_harmonicKinRadii;
extern int* d_harmonicKin;

extern int reactnum[35];

//MDS 
extern MDSystem mds;
//DCD
extern DCD dcd;
//FUNCTIONS
void MechanicalINIT();
