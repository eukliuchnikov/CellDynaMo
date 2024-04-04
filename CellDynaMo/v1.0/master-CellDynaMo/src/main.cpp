/*
 * mitosis model
 *
 *  Created on: Mar 15, 2018
 *	Edited on:	Mar 30, 2019
 *      Author: kliuchnikov
 */

#include "Initializiation/init.h"
#include "Initializiation/initSV.h"
#include "Potentials/langevin.cuh"
#include "Chemistry/react.h"
#include <time.h>

int main(){
//PART 1: SYSTEM CONFIGURATION
	loadConfig(config);
	cudaSetDevice(DEVICE);
	parameters(config);
	if (MEM_SHAPE == 0){
		printf("INCORRECT SHAPE PARAMETER\n");
		exit(0);
	}	
    if (RING == 1 && RAND_PLACE == 1){
        printf("DETERMINE INITIAL POSITION. BOTH, RANDOM AND RING POSITIONING SELECTED\n");
		exit(0);
    }
	clock_t timeMC;
	clock_t timeLD;
	clock_t timer;
	timeMC = clock();
	timeLD = clock();
	printf("Monte-Carlo time step is %f\n", tau);
	char gillespie[1024];	
	FILE *outp_au;
	FILE *outp_ndc;
	FILE *outp_att;
//PART 2: Physical initialisation
	MechanicalINIT();
//PART 3: Subvolume initialisation
	implicitINIT();
//PART 4: DCD-file creation
	createDCD(&dcd, mds.N, 11, 0, 1.0, 1, 0, 0, 0, 0);
  	const char* dcd_name = config.outputdcd.c_str();
  	char s2[config.outputdcd.size()+1];
  	strcpy(s2, dcd_name);
	dcdOpenWrite(&dcd, dcd_name);
	dcdWriteHeader(dcd);
	//Memory allocation for DCD-frames
	dcd.frame.X = (float*)calloc(mds.N, sizeof(float));
	dcd.frame.Y = (float*)calloc(mds.N, sizeof(float));
	dcd.frame.Z = (float*)calloc(mds.N, sizeof(float));
//PART 5: HARMONIC INIT
	cudaMalloc((void**)&d_kin_cos, KIN_COUNT*kn*sizeof(float));
    cudaMemcpy(d_kin_cos, kin_cos, KIN_COUNT*kn*sizeof(float), cudaMemcpyHostToDevice);
    cudaMalloc((void**)&d_dyn_cos, KIN_COUNT*kn*sizeof(float));
    cudaMemcpy(d_dyn_cos, dyn_cos, KIN_COUNT*kn*sizeof(float), cudaMemcpyHostToDevice);
    cudaMalloc((void**)&d_kin_rad, KIN_COUNT*kn*sizeof(float));
    cudaMemcpy(d_kin_rad, kin_rad, KIN_COUNT*kn*sizeof(float), cudaMemcpyHostToDevice);

    cudaMalloc((void**)&d_harmonicKinCount, KIN_COUNT*kn*sizeof(int));
    cudaMemcpy(d_harmonicKinCount, a_harmonicKinCount, KIN_COUNT*kn*sizeof(int), cudaMemcpyHostToDevice);
    cudaMalloc((void**)&d_harmonicKinRadii, KIN_COUNT*kn*maxHarmonicKinPerMonomer*sizeof(float));
    cudaMemcpy(d_harmonicKinRadii, a_harmonicKinRadii, KIN_COUNT*kn*maxHarmonicKinPerMonomer*sizeof(float), cudaMemcpyHostToDevice);
    cudaMalloc((void**)&d_harmonicKin, KIN_COUNT*kn*maxHarmonicKinPerMonomer*sizeof(int));
    cudaMemcpy(d_harmonicKin, a_harmonicKin, KIN_COUNT*kn*maxHarmonicKinPerMonomer*sizeof(int), cudaMemcpyHostToDevice);
//PART 6: REACTIONS&DIFFUSION (MAIN CYCLE)
	float time = 0.0;
	float ct = DCD_FREQ;
	float ot = 10.0*DCD_FREQ;
	int steps_num = 0;	
	DCD_step(dcd);
	DCD_step(dcd);

    if (RAND_PLACE == 1){
        minimiz(dcd);
	    printf("Minimization done\n");
    }
    else if (RAND_PLACE == 0){
        //dynamics(dcd, time);
    }

	int n_chrom = 2*floor((CONTOUR_LENGTH/2 - KIN_RADIUS)/(2*KIN_RADIUS) + 0.5);
	while (time < config.total_t){
		//DCD_step(dcd);
		timer = clock();
		int minimiz = 0;
		//diffusion(mds.h_r, cell_type, aurora_B_conc, aurora_A_conc, config.aurora_diff, neighbornum, neighbors, Aneighbornum, Aneighbors, SVcenter);
		reaction(time, outp_att, reactnum, aurora_B_conc, aurora_A_conc, NDC80_0_conc, NDC80_1_conc, NDC80_2_conc, NDC80_3_conc, NDC80_4_conc, NDC80_5_conc, NDC80_6_conc, NDC80_7_conc, NDC80_0_ATT_conc, NDC80_1_ATT_conc, NDC80_2_ATT_conc, NDC80_3_ATT_conc, NDC80_4_ATT_conc, NDC80_5_ATT_conc, NDC80_6_ATT_conc, NDC80_7_ATT_conc, mt_pol_conc, mt_depol_conc, mt_det_conc, minimiz);
		NDC80_dynamic(cell_type, mds.h_r, NDC80_0_conc, NDC80_1_conc, NDC80_2_conc, NDC80_3_conc, NDC80_4_conc, NDC80_5_conc, NDC80_6_conc, NDC80_7_conc, NDC80_0_ATT_conc, NDC80_1_ATT_conc, NDC80_2_ATT_conc, NDC80_3_ATT_conc, NDC80_4_ATT_conc, NDC80_5_ATT_conc, NDC80_6_ATT_conc, NDC80_7_ATT_conc);
		mt_dynamic(cell_type, mds.h_r, mt_pol_conc, mt_depol_conc, mt_det_conc, mt_att_conc, mds.h_type);
		timer = clock() - timer;
		timeMC += timer;
		if (minimiz == 1){
			timer = clock();
			
			for (int m = 0; m < KIN_COUNT; m += 2){
				float3 rcom1 = mds.h_r[POLE_COUNT*MT_NUM*MAX_MT_LENGTH + (m + 1)*kt[0].n - n_chrom - 1];
				float3 rcom2 = mds.h_r[POLE_COUNT*MT_NUM*MAX_MT_LENGTH + (m + 2)*kt[1].n - n_chrom - 1];

				float3 kincenter;
				kincenter.x = (rcom1.x + rcom2.x)/2;
				kincenter.y = (rcom1.y + rcom2.y)/2;
				kincenter.z = (rcom1.z + rcom2.z)/2;

				int number = (m + 2)/2 - 1;

				int3 null = get_cell(COORtransfer1(kincenter));
				int n = COORintoARRAY(null.x, null.y, null.z);
				SV_centers_before[number] = n;
			}		
			dynamics(dcd, time);
			int diff_trigger = 0;
			for (int m = 0; m < KIN_COUNT; m += 2){
				float3 rcom1 = mds.h_r[POLE_COUNT*MT_NUM*MAX_MT_LENGTH + (m + 1)*kt[0].n - n_chrom - 1];
				float3 rcom2 = mds.h_r[POLE_COUNT*MT_NUM*MAX_MT_LENGTH + (m + 2)*kt[1].n - n_chrom - 1];

				float3 kincenter;
				kincenter.x = (rcom1.x + rcom2.x)/2;
				kincenter.y = (rcom1.y + rcom2.y)/2;
				kincenter.z = (rcom1.z + rcom2.z)/2;

				int number = (m + 2)/2 - 1;

				int3 null = get_cell(COORtransfer1(kincenter));
				int n = COORintoARRAY(null.x, null.y, null.z);
				SV_centers_after[number] = n;
				if (SV_centers_after[number] != SV_centers_before[number]){
					diff_trigger = 1;
				}
			}
			if (diff_trigger == 1){
				int NABK = number_AB[0];
				count_dyn(mds.h_r, cell_type, aurora_B_conc, mds.N, NABK);
			}
			//diffusionBar(cell_type, neighbors, neighbornum,  neighbornum_X,  neighbornum_Y,  neighbornum_Z,  SVcenter);
			//AdiffusionBar(cell_type, Aneighbors, Aneighbornum,  Aneighbornum_X,  Aneighbornum_Y,  Aneighbornum_Z,  SVcenter);
			timer = clock() - timer;
			timeLD += timer;
			float inp = (double)timeMC/CLOCKS_PER_SEC + (double)timeLD/CLOCKS_PER_SEC;
			int days = inp/(3600*24);
			int rest = inp - days*3600*24;
			int hours = rest/3600;
			rest = rest - hours*3600;
			int minutes = rest/60;
			int seconds = rest - minutes*60;

			float inpMC = (double)timeMC/CLOCKS_PER_SEC;
			int daysMC = inpMC/(3600*24);
			int restMC = inpMC - daysMC*3600*24;
			int hoursMC = restMC/3600;
			restMC = restMC - hoursMC*3600;
			int minutesMC = restMC/60;
			int secondsMC = restMC - minutesMC*60;

			float MC_perc = 100.0*inpMC/inp;

			float inpLD = (double)timeLD/CLOCKS_PER_SEC;
			int daysLD = inpLD/(3600*24);
			int restLD = inpLD - daysLD*3600*24;
			int hoursLD = restLD/3600;
			restLD = restLD - hoursLD*3600;
			int minutesLD = restLD/60;
			int secondsLD = restLD - minutesLD*60;

			float LD_perc = 100.0*inpLD/inp;

			int b_hours = time/3600;
			int b_rest = time - b_hours*3600;
			int b_minutes = b_rest/60;
			int b_seconds = b_rest - b_minutes*60;
			printf("***************************************************************\n");
			printf("Biological time: %f s (%d h %d min %d s)\n",time, b_hours, b_minutes, b_seconds);
			printf("Computational time: %d days %d h %d min %d s (Monte Carlo: %f%%, Langevin: %f%%)\n", days, hours, minutes, seconds, MC_perc, LD_perc);
			printf("(Monte Carlo time: %d days %d h %d min %d s (%f%%))\n", daysMC, hoursMC, minutesMC, secondsMC, MC_perc);
			printf("(Langevin time: %d days %d h %d min %d s (%f%%))\n", daysLD, hoursLD, minutesLD, secondsLD, LD_perc);
			printf("Number of Monte Carlo steps: %d\t (%f steps/sec)\n", steps_num, steps_num/inp);
			inp = (config.total_t - time)*((double)timeMC/CLOCKS_PER_SEC + (double)timeLD/CLOCKS_PER_SEC)/time;
	 	 	days = inp/(3600*24);
			rest = inp - days*3600*24;
			hours = rest/3600;
			rest = rest - hours*3600;
			minutes = rest/60;
			seconds = rest - minutes*60;
			printf("Estimated time left: %d days %d h %d min %d s\n", days, hours, minutes, seconds);
			printf("***************************************************************\n");
			NDC80_dynamic(cell_type, mds.h_r, NDC80_0_conc, NDC80_1_conc, NDC80_2_conc, NDC80_3_conc, NDC80_4_conc, NDC80_5_conc, NDC80_6_conc, NDC80_7_conc, NDC80_0_ATT_conc, NDC80_1_ATT_conc, NDC80_2_ATT_conc, NDC80_3_ATT_conc, NDC80_4_ATT_conc, NDC80_5_ATT_conc, NDC80_6_ATT_conc, NDC80_7_ATT_conc);
			mt_dynamic(cell_type, mds.h_r, mt_pol_conc, mt_depol_conc, mt_det_conc, mt_att_conc, mds.h_type);
		}
		minimiz = 0;
		if (time > ct){
			DCD_step(dcd);
			//int number = int(ct);
		    Au_out(time, cell_type, outp_au, aurora_B_conc, aurora_A_conc);
			/*if (MAP == 1){
				map(cell_type, aurora_B_conc, aurora_A_conc, NDC80_0_conc, gillespie, number);
			}*/
			ct += DCD_FREQ;
		}
		if (time > ot){
			ndc_out(time, outp_au, NDC80_0_conc, NDC80_1_conc, NDC80_2_conc, NDC80_3_conc, NDC80_4_conc, NDC80_5_conc, NDC80_6_conc, NDC80_7_conc);
			ot += 10.0*DCD_FREQ;
		}
		time += tau;
		steps_num++;
	}
	return 0;
}
