#pragma once
#include "../Math/mat.h"
#include "react.h"

void diffusion(float3* r0, int* cell_type, int* Anum, int* Bnum, int D, int* neighbornum, int* neighbors, int* Aneighbornum, int* Aneighbors, float3* SVcenter);
