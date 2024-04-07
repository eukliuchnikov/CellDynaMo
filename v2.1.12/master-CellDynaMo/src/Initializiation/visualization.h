#include "init.h"
#include "../TCO_Writers/pdbio.h"
#include "../TCO_Writers/xyzio.h"
#include "../TCO_Writers/psfio.h"
#include "../TCO_Writers/dcdio.h"

#include "initSV.h"

void membrane_visualization();
int initPSF(KT* kt);
void DCD_step(DCD dcd);

void map(int* cell_type, int* aurora_B_conc, int* aurora_A_conc, int* NDC80_0_conc, char* gillespie, int number);
void Au_out(float time, int* cell_type, FILE* outp, int* aurora_B_conc, int* aurora_A_conc);
void ndc_out(float time, FILE* outp, int* NDC80_0_conc, int* NDC80_1_conc, int* NDC80_2_conc, int* NDC80_3_conc, int* NDC80_4_conc, int* NDC80_5_conc, int* NDC80_6_conc, int* NDC80_7_conc);
void att_out(float time, FILE* outp, int triger, int mt_id, int kt_id);
