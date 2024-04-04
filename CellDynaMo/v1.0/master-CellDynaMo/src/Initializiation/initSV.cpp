#include "initSV.h"

int* aurora_B_conc;
int* aurora_A_conc;
int* phosphatase_conc;

int* NDC80_0_conc;
int* NDC80_1_conc;
int* NDC80_2_conc;
int* NDC80_3_conc;
int* NDC80_4_conc;
int* NDC80_5_conc;
int* NDC80_6_conc;
int* NDC80_7_conc;

int* NDC80_0_ATT_conc;
int* NDC80_1_ATT_conc;
int* NDC80_2_ATT_conc;
int* NDC80_3_ATT_conc;
int* NDC80_4_ATT_conc;
int* NDC80_5_ATT_conc;
int* NDC80_6_ATT_conc;
int* NDC80_7_ATT_conc;

int* mt_pol_conc;
int* mt_depol_conc;

int* mt_det_conc;
int* mt_att_conc;

int* SV_centers_before;
int* SV_centers_after;

int* cell_type;

float3* SVcenter;

int* number_AB;

int* neighbornum_X;
int* neighbornum_Y;
int* neighbornum_Z;
int* neighbornum;
int* neighbors;

int* Aneighbornum_X;
int* Aneighbornum_Y;
int* Aneighbornum_Z;
int* Aneighbornum;
int* Aneighbors;

void implicitINIT(){
	aurora_B_conc = (int*)calloc(SVn, sizeof(int));
	aurora_A_conc = (int*)calloc(SVn, sizeof(int));
	phosphatase_conc = (int*)calloc(SVn, sizeof(int));

	NDC80_0_conc = (int*)calloc(SVn, sizeof(int));
	NDC80_1_conc = (int*)calloc(SVn, sizeof(int));
	NDC80_2_conc = (int*)calloc(SVn, sizeof(int));
	NDC80_3_conc = (int*)calloc(SVn, sizeof(int));
	NDC80_4_conc = (int*)calloc(SVn, sizeof(int));
	NDC80_5_conc = (int*)calloc(SVn, sizeof(int));
	NDC80_6_conc = (int*)calloc(SVn, sizeof(int));
	NDC80_7_conc = (int*)calloc(SVn, sizeof(int));

	NDC80_0_ATT_conc = (int*)calloc(SVn, sizeof(int));
	NDC80_1_ATT_conc = (int*)calloc(SVn, sizeof(int));
	NDC80_2_ATT_conc = (int*)calloc(SVn, sizeof(int));
	NDC80_3_ATT_conc = (int*)calloc(SVn, sizeof(int));
	NDC80_4_ATT_conc = (int*)calloc(SVn, sizeof(int));
	NDC80_5_ATT_conc = (int*)calloc(SVn, sizeof(int));
	NDC80_6_ATT_conc = (int*)calloc(SVn, sizeof(int));
	NDC80_7_ATT_conc = (int*)calloc(SVn, sizeof(int));

	mt_pol_conc = (int*)calloc(SVn, sizeof(int));
	mt_depol_conc = (int*)calloc(SVn, sizeof(int));

	mt_det_conc = (int*)calloc(SVn, sizeof(int));
	mt_att_conc = (int*)calloc(SVn, sizeof(int));

	SV_centers_before = (int*)calloc(KIN_COUNT/2, sizeof(int));
	SV_centers_after = (int*)calloc(KIN_COUNT/2, sizeof(int));

	cell_type = (int*)calloc(SVn, sizeof(int));

	SVcenter = (float3*)calloc(SVn, sizeof(float3));

	neighbornum_X = (int*)calloc(SVn, sizeof(int));
	neighbornum_Y = (int*)calloc(SVn, sizeof(int));
	neighbornum_Z = (int*)calloc(SVn, sizeof(int));
	neighbornum = (int*)calloc(SVn, sizeof(int));
	neighbors = (int*)calloc(SVn*6, sizeof(int));

	Aneighbornum_X = (int*)calloc(SVn, sizeof(int));
	Aneighbornum_Y = (int*)calloc(SVn, sizeof(int));
	Aneighbornum_Z = (int*)calloc(SVn, sizeof(int));
	Aneighbornum = (int*)calloc(SVn, sizeof(int));
	Aneighbors = (int*)calloc(SVn*6, sizeof(int));

	number_AB = (int*)calloc(1, sizeof(int));
	
	membrane_impl(cell_type);
	mt_dynamic(cell_type, mds.h_r, mt_pol_conc, mt_depol_conc, mt_det_conc, mt_att_conc, mds.h_type);
	pole_initial(cell_type, mds.h_r);
	NDC80_dynamic(cell_type, mds.h_r, NDC80_0_conc, NDC80_1_conc, NDC80_2_conc, NDC80_3_conc, NDC80_4_conc, NDC80_5_conc, NDC80_6_conc, NDC80_7_conc, NDC80_0_ATT_conc, NDC80_1_ATT_conc, NDC80_2_ATT_conc, NDC80_3_ATT_conc, NDC80_4_ATT_conc, NDC80_5_ATT_conc, NDC80_6_ATT_conc, NDC80_7_ATT_conc);
	count_initial(mds.h_r, cell_type, aurora_B_conc, aurora_A_conc, mds.N, number_AB);
	diffusionBar(cell_type, neighbors, neighbornum,  neighbornum_X,  neighbornum_Y,  neighbornum_Z,  SVcenter);
	AdiffusionBar(cell_type, Aneighbors, Aneighbornum,  Aneighbornum_X,  Aneighbornum_Y,  Aneighbornum_Z,  SVcenter);

	/*char gillespie[1024];
	int number = 0;
	if (MAP == 1){
		map(cell_type, aurora_B_conc, aurora_A_conc, NDC80_0_conc, gillespie, number);
	}*/
}