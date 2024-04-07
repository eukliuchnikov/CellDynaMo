#pragma once
#include "../Math/mat.h"
#include "react.h"

void kinesin_attach(float3* r, int* type, int l, float3* center, int* connect_1);
void kinesin_detach(float3* r, int* type, int l, float3* center, int* connect_1, int* connect_2, float* len1, float* len2, float* att_l);
void kinesin_step(float3* r, int* type, float* length, int l, float3* center, int* connect_1, int* connect_2, float* len1, float* len2, float* f_out, float* att_l);
