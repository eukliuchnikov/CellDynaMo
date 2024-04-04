#pragma once
#include "initSV.h"

void membrane_impl(int* cell_type);
void mt_dynamic(int* cell_type, float3* r, int* mt_pol_conc, int* mt_depol_conc, int* mt_det_conc, int* mt_att_conc, int* type);
void pole_initial(int* cell_type, float3* r);
void NDC80_dynamic(int* cell_type, float3* r, int* NDC80_0_conc, int* NDC80_1_conc, int* NDC80_2_conc, int* NDC80_3_conc, int* NDC80_4_conc, int* NDC80_5_conc, int* NDC80_6_conc, int* NDC80_7_conc, int* NDC80_0_ATT_conc, int* NDC80_1_ATT_conc, int* NDC80_2_ATT_conc, int* NDC80_3_ATT_conc, int* NDC80_4_ATT_conc, int* NDC80_5_ATT_conc, int* NDC80_6_ATT_conc, int* NDC80_7_ATT_conc);
void count_initial(float3* r0, int* cell_type, int* aurora_B_conc, int* aurora_A_conc, int N, int* number_AB);
void count_dyn(float3* r0, int* cell_type, int* aurora_B_conc, int N, int* number_AB);
