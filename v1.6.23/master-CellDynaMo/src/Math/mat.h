#pragma once
#include "math.h"
#include "../Initializiation/init.h"
#include "../Initializiation/initSV.h"
#include "../Potentials/langevin.cuh"

int computeGamma(float r);
float myDistance(float3 ri, float3 rj);
float length(float3 r);
float scalar(float3 ri, float3 rj);

float3 get_coord(int3 cell);
int3 get_cell(float3 r);
int COORintoARRAY(int x, int y, int z);
int3 ARRAYintoCOOR(int N);
float3 COORtransfer(float3 r);
float3 COORtransfer1(float3 r);
