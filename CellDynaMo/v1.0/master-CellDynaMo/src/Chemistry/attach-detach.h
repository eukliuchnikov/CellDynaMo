#pragma once
#include "../Math/mat.h"
#include "react.h"

void AttachDettach(int attType, FILE* outp_att, int l, float3* r, int* type, int reactnum[35], int mu, int reactType, float time, int* connector, float* att_l, float* len, float* link_l);
