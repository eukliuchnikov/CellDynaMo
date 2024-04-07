#pragma once
#include "init.h"

using namespace std;

std::vector<kin> generateKinetochoreHollow();
void KT_seeds(KT* kt, float3* r, vector <kin> kinetochore, int* type);
void KT_bonds(KT* kt, int** harmonicKinCount, int** harmonicKin, float** harmonicKinRadii);
void KT_rotation(KT* kt, Config& config);
void KT_rotationRAND(KT* kt, Config& config);
