#pragma once
#include "init.h"
#include "../Math/mat.h"
#include "visualization.h"
#include "implicit.h"
#include "diff_update.h"
#include "configreader.h"

extern int* aurora_B_conc;
extern int* aurora_A_conc;
extern int* phosphatase_conc;

extern int* NDC80_0_conc;
extern int* NDC80_1_conc;
extern int* NDC80_2_conc;
extern int* NDC80_3_conc;
extern int* NDC80_4_conc;
extern int* NDC80_5_conc;
extern int* NDC80_6_conc;
extern int* NDC80_7_conc;

extern int* NDC80_0_ATT_conc;
extern int* NDC80_1_ATT_conc;
extern int* NDC80_2_ATT_conc;
extern int* NDC80_3_ATT_conc;
extern int* NDC80_4_ATT_conc;
extern int* NDC80_5_ATT_conc;
extern int* NDC80_6_ATT_conc;
extern int* NDC80_7_ATT_conc;

extern int* mt_pol_conc;
extern int* mt_depol_conc;
extern int* mt_det_conc;
extern int* mt_att_conc;

extern int* kinesin_att_conc;
extern int* kinesin_det_conc;

extern int* SV_centers_before;
extern int* SV_centers_after;

extern int* cell_type;

extern float3* SVcenter;

extern int* neighbornum_X;
extern int* neighbornum_Y;
extern int* neighbornum_Z;
extern int* neighbornum;
extern int* neighbors;

extern int* Aneighbornum_X;
extern int* Aneighbornum_Y;
extern int* Aneighbornum_Z;
extern int* Aneighbornum;
extern int* Aneighbors;

extern int* number_AB;

void implicitINIT();
