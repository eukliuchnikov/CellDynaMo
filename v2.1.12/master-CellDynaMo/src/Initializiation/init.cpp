#include "init.h"

using namespace std;

int seed = -rseed;
//COMPONENTS
MT* pol[POLE_COUNT];
KT kt[KIN_COUNT];
float rk;
int kn;
float* kin_cos;
float* dyn_cos;
float* kin_rad;
int* harmonicKinCount[KIN_COUNT];
int* harmonicKin[KIN_COUNT];
float* harmonicKinRadii[KIN_COUNT];

int* a_harmonicKinCount;
int* a_harmonicKin;
float* a_harmonicKinRadii;

float* d_kin_cos;
float* d_kin_rad;
float* d_dyn_cos;
int* d_harmonicKinCount;
float* d_harmonicKinRadii;
int* d_harmonicKin;

//MDS 
MDSystem mds;
//DCD
DCD dcd;
//INITIATION
void MechanicalINIT(){
	membrane_visualization();
	vector <kin> kinetochore = generateKinetochoreHollow();

	kn = kinetochore.size()/KIN_COUNT;
	//MDS initialization
	mds.N = POLE_COUNT*MT_NUM*MAX_MT_LENGTH + POLE_COUNT + kinetochore.size() + 2*POLE_COUNT*MT_NUM;	
	mds.h_r = (float3*)calloc(mds.N, sizeof(float3));
	mds.h_f = (float3*)calloc(mds.N, sizeof(float3));
	mds.h_type = (int*)calloc(mds.N, sizeof(int));
	mds.h_length = (float*)calloc(mds.N, sizeof(float));
	mds.h_connector = (int*)calloc(mds.N, sizeof(int));
	mds.h_att_l = (float*)calloc(mds.N, sizeof(float));
	mds.h_att_n = (int*)calloc(mds.N, sizeof(int));
	mds.h_att_f = (float*)calloc(mds.N, sizeof(float));
	mds.h_link_l = (float*)calloc(mds.N, sizeof(float));
	mds.h_connect_1 = (int*)calloc(mds.N, sizeof(int));
	mds.h_connect_2 = (int*)calloc(mds.N, sizeof(int));
	mds.h_len_1 = (float*)calloc(mds.N, sizeof(float));
	mds.h_len_2 = (float*)calloc(mds.N, sizeof(float));
	mds.h_dom_center = (float3*)calloc(mds.N, sizeof(float3));
	mds.h_f_out = (float*)calloc(mds.N, sizeof(float));

	cudaMalloc((void**)&mds.d_r, mds.N*sizeof(float3));
	cudaMalloc((void**)&mds.d_f, mds.N*sizeof(float3));
	cudaMalloc((void**)&mds.d_type, mds.N*sizeof(int));
	cudaMalloc((void**)&mds.d_length, mds.N*sizeof(float));
	cudaMalloc((void**)&mds.d_connector, mds.N*sizeof(int));
	cudaMalloc((void**)&mds.d_att_l, mds.N*sizeof(float));
	cudaMalloc((void**)&mds.d_att_n, mds.N*sizeof(int));
	cudaMalloc((void**)&mds.d_att_f, mds.N*sizeof(float));
	cudaMalloc((void**)&mds.d_link_l, mds.N*sizeof(float));
	cudaMalloc((void**)&mds.d_connect_1, mds.N*sizeof(int));
	cudaMalloc((void**)&mds.d_connect_2, mds.N*sizeof(int));
	cudaMalloc((void**)&mds.d_len_1, mds.N*sizeof(float));
	cudaMalloc((void**)&mds.d_len_2, mds.N*sizeof(float));
	cudaMalloc((void**)&mds.d_dom_center, mds.N*sizeof(float3));
	cudaMalloc((void**)&mds.d_f_out, mds.N*sizeof(float));

	int i;
	for (i = 0; i < mds.N; i++){
		mds.h_f[i].x = 0.0;
		mds.h_f[i].y = 0.0;
		mds.h_f[i].z = 0.0;
		mds.h_type[i] = 0;
		mds.h_length[i] = 0.0;
		mds.h_connector[i] = 0;
		mds.h_connect_1[i] = 0;
		mds.h_connect_2[i] = 0;
		mds.h_len_1[i] = 0.0;
		mds.h_len_2[i] = 0.0;
		mds.h_dom_center[i].x = 0.0;
		mds.h_dom_center[i].y = 0.0;
		mds.h_dom_center[i].z = 0.0;
		mds.h_f_out[i] = 0.0;
	}
	poles(mds.h_r, mds.N, mds.h_type);
	MT_seeds(pol, mds.h_type, mds.h_length);
	KT_seeds(kt, mds.h_r, kinetochore, mds.h_type);
	kinesin_init(mds.h_r, mds.h_type);

	kin_cos =  (float*)calloc(KIN_COUNT*kn, sizeof(float));
	dyn_cos =  (float*)calloc(KIN_COUNT*kn, sizeof(float));
	kin_rad =  (float*)calloc(KIN_COUNT*kn, sizeof(float));
	int h;
	for (h = 0; h < KIN_COUNT; h++){
		for (i = 0; i < kn; i++){
			kin_cos[h*kn + i] = kt[h].cos[i];
			dyn_cos[h*kn + i] = kt[h].cos[i];
			kin_rad[h*kn + i] = kt[h].radius[i];
		}
		if (h % 2 == 0){
			int n_chrom = 2*floor((CONTOUR_LENGTH/2 - KIN_RADIUS)/(2*KIN_RADIUS) + 0.5);
			mds.h_type[POLE_COUNT*MT_NUM*MAX_MT_LENGTH + (h + 1)*kn - n_chrom - 1] = LEFT_KT;
			if (CHROM_ARMS == 1){
				for (int j = 0; j < n_chrom; j++){
					mds.h_type[POLE_COUNT*MT_NUM*MAX_MT_LENGTH + (h + 1)*kn - n_chrom + j] = CHROM;
				}
			}
		}
		else	{
			int n_chrom = 2*floor((CONTOUR_LENGTH/2 - KIN_RADIUS)/(2*KIN_RADIUS) + 0.5);
			mds.h_type[POLE_COUNT*MT_NUM*MAX_MT_LENGTH + (h + 1)*kn - n_chrom - 1] = RIGHT_KT;
			if (CHROM_ARMS == 1){
				for (int j = 0; j < n_chrom; j++){
					mds.h_type[POLE_COUNT*MT_NUM*MAX_MT_LENGTH + (h + 1)*kn - n_chrom + j] = CHROM;
				}
			}
		}
		
		harmonicKinCount[h] = (int*)calloc(kt[0].n, sizeof(int));
		harmonicKin[h] = (int*)malloc(((kt[0].n)*maxHarmonicKinPerMonomer)*sizeof(int));
		harmonicKinRadii[h] = (float*)malloc(((kt[0].n)*maxHarmonicKinPerMonomer)*sizeof(float));
	}
	KT_bonds(kt, harmonicKinCount, harmonicKin, harmonicKinRadii);
	a_harmonicKinCount = (int*)calloc(KIN_COUNT*kn, sizeof(int));
	a_harmonicKinRadii = (float*)calloc(KIN_COUNT*kn*maxHarmonicKinPerMonomer, sizeof(float));
	a_harmonicKin = (int*)calloc(KIN_COUNT*kn*maxHarmonicKinPerMonomer, sizeof(int));

	for (h = 0; h < KIN_COUNT; h++){
		for (i = 0; i < kn; i++){
			a_harmonicKinCount[h*kn + i] = harmonicKinCount[h][i];
			int m;
			for (m = 0; m < maxHarmonicKinPerMonomer; m++){
				a_harmonicKinRadii[h*kn*maxHarmonicKinPerMonomer + i*maxHarmonicKinPerMonomer + m] = harmonicKinRadii[h][i*maxHarmonicKinPerMonomer + m];
				a_harmonicKin[h*kn*maxHarmonicKinPerMonomer + i*maxHarmonicKinPerMonomer + m] = harmonicKin[h][i*maxHarmonicKinPerMonomer + m];
				//printf("%d\t%d\t%d\n", h*kn*maxHarmonicKinPerMonomer + i*maxHarmonicKinPerMonomer + m, a_harmonicKin[h*kn*maxHarmonicKinPerMonomer + i*maxHarmonicKinPerMonomer + m], a_harmonicKinCount[h*kn + i]);
			} 
		}
	}
	KT_rotationRAND(kt, config);

	initPSF(kt);

	mds.h_att_n[0] = 0; //total left
    mds.h_att_n[1] = 0;	//total right
    mds.h_att_f[0] = 0.0;	//force left
    mds.h_att_f[1] = 0.0;	//force right

    mds.h_att_n[2] = 0; //push left
    mds.h_att_n[3] = 0;	//push right
    mds.h_att_n[4] = 0; //pull left
    mds.h_att_n[5] = 0;	//pull right
}
