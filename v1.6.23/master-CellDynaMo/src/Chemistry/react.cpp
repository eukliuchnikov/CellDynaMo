/*
 *NSV Gillespie module
 *
 *  Created on: Mar 15, 2018
 *	Edited on:	Mar 30, 2019
 *      Author: kliuchnikov
 */

#include "react.h"
#include "../Initializiation/init.h"
#include "../Initializiation/initSV.h"

int rseed_react = rseed*3;
int rseed_mu = rseed + 12345;

void reaction(float time, FILE* outp_att, int reactnum[35], int* Nab, int* Naa, int* Nn0, int* Nn1, int* Nn2, int* Nn3, int* Nn4, int* Nn5, int* Nn6, int* Nn7, int* Nn0_A, int* Nn1_A, int* Nn2_A, int* Nn3_A, int* Nn4_A, int* Nn5_A, int* Nn6_A, int* Nn7_A, int* Nmtpol, int* Nmtdepol, int* Nmtdet, int & minimiz){
	float c[35];
	float hi[35];

	c[0] = 0.0;
	c[1] = K_PHOS;		//NOT REAL
	c[2] = K_PHOS;
	c[3] = K_PHOS;
	c[4] = K_PHOS;
	c[5] = K_PHOS;
	c[6] = K_PHOS;
	c[7] = K_PHOS;

	c[8] = K_DEPHOS;		//NOT REAL
	c[9] = K_DEPHOS;
	c[10] = K_DEPHOS;
	c[11] = K_DEPHOS;
	c[12] = K_DEPHOS;
	c[13] = K_DEPHOS;
	c[14] = K_DEPHOS;

	c[15] = K_GROWTH;
	c[16] = K_SHORTEN;
	c[17] = K_RESCUE;
	c[18] = K_CAT;

	c[19] = K_ATTACH*trans/(sv_size*sv_size*sv_size);
	c[20] = K_ATTACH*trans/(sv_size*sv_size*sv_size);
	c[21] = K_ATTACH*trans/(sv_size*sv_size*sv_size);
	c[22] = K_ATTACH*trans/(sv_size*sv_size*sv_size);
	c[23] = K_ATTACH*trans/(sv_size*sv_size*sv_size);
	c[24] = K_ATTACH*trans/(sv_size*sv_size*sv_size);
	c[25] = K_ATTACH*trans/(sv_size*sv_size*sv_size);
	c[26] = K_ATTACH*trans/(sv_size*sv_size*sv_size);

	c[27] = K_DETACH_0;
	c[28] = K_DETACH_1;
	c[29] = K_DETACH_2;
	c[30] = K_DETACH_3;
	c[31] = K_DETACH_4;
	c[32] = K_DETACH_5;
	c[33] = K_DETACH_6;
	c[34] = K_DETACH_7;

	hi[0] = 0.0;
	int l, n;
	int i, j, k;
	for (i = 1; i < num_x_cell + 1; i++){
		for (j = 1; j < num_y_cell + 1; j++){
			for (k = 1; k < num_z_cell + 1; k++){
				l = COORintoARRAY(i, j, k);
				int mu = 0;
		
				hi[1] = Nab[l]*Nn0[l];					//phos1
				hi[2] = Nab[l]*Nn1[l];					//phos2
				hi[3] = Nab[l]*Nn2[l];					//phos3
				hi[4] = Nab[l]*Nn3[l];					//phos4
				hi[5] = Nab[l]*Nn4[l];					//phos5
				hi[6] = Nab[l]*Nn5[l];					//phos6
				hi[7] = Nab[l]*Nn6[l];					//phos7
					
				hi[8] = Nn7[l];					//dephos7
				hi[9] = Nn6[l];					//dephos6
				hi[10] = Nn5[l];					//dephos5
				hi[11] = Nn4[l];					//dephos4
				hi[12] = Nn3[l];					//dephos3
				hi[13] = Nn2[l];					//dephos2
				hi[14] = Nn1[l];					//dephos1

				hi[15] = Nmtpol[l];					//MT add
				hi[16] = Nmtdepol[l];					//MT remove
				hi[17] = Nmtdepol[l];					//MT rescue change
				hi[18] = Nmtpol[l];					//MT catastrophe change

				hi[19] = Nmtdet[l]*Nn0[l];	//MT and NDC-80 association
				hi[20] = Nmtdet[l]*Nn1[l];
				hi[21] = Nmtdet[l]*Nn2[l];
				hi[22] = Nmtdet[l]*Nn3[l];
				hi[23] = Nmtdet[l]*Nn4[l];
				hi[24] = Nmtdet[l]*Nn5[l];
				hi[25] = Nmtdet[l]*Nn6[l];
				hi[26] = Nmtdet[l]*Nn7[l];
				
				hi[27] = Nn0_A[l];	//MT and NDC-80 (0 phosphates) dissociation
				hi[28] = Nn1_A[l];
				hi[29] = Nn2_A[l];
				hi[30] = Nn3_A[l];
				hi[31] = Nn4_A[l];
				hi[32] = Nn5_A[l];
				hi[33] = Nn6_A[l];
				hi[34] = Nn7_A[l];
				
				int r_id;
				float alpha0 = 0;
				for (r_id = 1; r_id < 35; r_id++){
					alpha0 += c[r_id]*hi[r_id];
				}
				float r1 = ran2(&rseed_react);
				if (r1 < 1 - exp(-alpha0*tau)){
					float r2 = ran2(&rseed_mu);
					float sum = 0;
					for (r_id = 1; r_id < 35; r_id++){
						sum += c[r_id]*hi[r_id];
						if (sum - c[r_id]*hi[r_id] < r2*alpha0 && r2*alpha0 < sum){
							mu = r_id;
						}
					}
					reactnum[mu]++;
					if (mu == 1){
						int reactType = NDC_0;
						int prodType = NDC_1;
						PhosDephos(Nn0, Nn1, reactType, prodType, l);
					}
					if (mu == 2){
						int reactType = NDC_1;
						int prodType = NDC_2;
						PhosDephos(Nn1, Nn2, reactType, prodType, l);
					}
					if (mu == 3){
						int reactType = NDC_2;
						int prodType = NDC_3;
						PhosDephos(Nn2, Nn3, reactType, prodType, l);
					}
					if (mu == 4){
						int reactType = NDC_3;
						int prodType = NDC_4;
						PhosDephos(Nn3, Nn4, reactType, prodType, l);
					}
					if (mu == 5){
						int reactType = NDC_4;
						int prodType = NDC_5;
						PhosDephos(Nn4, Nn5, reactType, prodType, l);
					}
					if (mu == 6){
						int reactType = NDC_5;
						int prodType = NDC_6;
						PhosDephos(Nn5, Nn6, reactType, prodType, l);
					}
					if (mu == 7){
						int reactType = NDC_6;
						int prodType = NDC_7;
						PhosDephos(Nn6, Nn7, reactType, prodType, l);
					}
///////////////////////////////////////////
					if (mu == 8){
						int reactType = NDC_7;
						int prodType = NDC_6;
						PhosDephos(Nn7, Nn6, reactType, prodType, l);
					}
					if (mu == 9){
						int reactType = NDC_6;
						int prodType = NDC_5;
						PhosDephos(Nn6, Nn5, reactType, prodType, l);
					}
					if (mu == 10){
						int reactType = NDC_5;
						int prodType = NDC_4;
						PhosDephos(Nn5, Nn4, reactType, prodType, l);
					}
					if (mu == 11){
						int reactType = NDC_4;
						int prodType = NDC_3;
						PhosDephos(Nn4, Nn3, reactType, prodType, l);
					}
					if (mu == 12){
						int reactType = NDC_3;
						int prodType = NDC_2;
						PhosDephos(Nn3, Nn2, reactType, prodType, l);
					}
					if (mu == 13){
						int reactType = NDC_2;
						int prodType = NDC_1;
						PhosDephos(Nn2, Nn1, reactType, prodType, l);
					}
					if (mu == 14){
						int reactType = NDC_1;
						int prodType = NDC_0;
						PhosDephos(Nn1, Nn0, reactType, prodType, l);
					}
///////////////////////////////////////////
					if (mu == 15){
						GROWTHfunc(l, mds.h_r, mds.h_type, minimiz, mds.h_length);
					}
					if (mu == 16){
						SHORTfunc(l, mds.h_r, mds.h_type, minimiz, mds.h_length, mds.h_att_l);
					}
					if (mu == 17){
						RESCfunc(l, mds.h_r, mds.h_type);
					}
					if (mu == 18){
						CATfunc(l, mds.h_r);
					}
///////////////////////////////////////////
					if (mu == 19){
						reactnum[mu]--;
						int attType = DETACH;
						int reactType = NDC_0; 
						AttachDettach(attType, outp_att, l, mds.h_r, mds.h_type, reactnum, mu, reactType, time, mds.h_connector, mds.h_att_l, mds.h_length, mds.h_link_l);
					}
					if (mu == 20){
						reactnum[mu]--;
						int attType = DETACH;
						int reactType = NDC_1; 
						AttachDettach(attType, outp_att, l, mds.h_r, mds.h_type, reactnum, mu, reactType, time, mds.h_connector, mds.h_att_l, mds.h_length, mds.h_link_l);
					}
					if (mu == 21){
						reactnum[mu]--;
						int attType = DETACH;
						int reactType = NDC_2; 
						AttachDettach(attType, outp_att, l, mds.h_r, mds.h_type, reactnum, mu, reactType, time, mds.h_connector, mds.h_att_l, mds.h_length, mds.h_link_l);
					}
					if (mu == 22){
						reactnum[mu]--;
						int attType = DETACH;
						int reactType = NDC_3; 
						AttachDettach(attType, outp_att, l, mds.h_r, mds.h_type, reactnum, mu, reactType, time, mds.h_connector, mds.h_att_l, mds.h_length, mds.h_link_l);
					}
					if (mu == 23){
						reactnum[mu]--;
						int attType = DETACH;
						int reactType = NDC_4; 
						AttachDettach(attType, outp_att, l, mds.h_r, mds.h_type, reactnum, mu, reactType, time, mds.h_connector, mds.h_att_l, mds.h_length, mds.h_link_l);
					}
					if (mu == 24){
						reactnum[mu]--;
						int attType = DETACH;
						int reactType = NDC_5; 
						AttachDettach(attType, outp_att, l, mds.h_r, mds.h_type, reactnum, mu, reactType, time, mds.h_connector, mds.h_att_l, mds.h_length, mds.h_link_l);
					}
					if (mu == 25){
						reactnum[mu]--;
						int attType = DETACH;
						int reactType = NDC_6; 
						AttachDettach(attType, outp_att, l, mds.h_r, mds.h_type, reactnum, mu, reactType, time, mds.h_connector, mds.h_att_l, mds.h_length, mds.h_link_l);
					}
					if (mu == 26){
						reactnum[mu]--;
						int attType = DETACH;
						int reactType = NDC_7; 
						AttachDettach(attType, outp_att, l, mds.h_r, mds.h_type, reactnum, mu, reactType, time, mds.h_connector, mds.h_att_l, mds.h_length, mds.h_link_l);
					}
					if (mu == 27){
						reactnum[mu]--;
						int attType = ATTACH;
						int reactType = NDC_0;
						AttachDettach(attType, outp_att, l, mds.h_r, mds.h_type, reactnum, mu, reactType, time, mds.h_connector, mds.h_att_l, mds.h_length, mds.h_link_l);
					}
					if (mu == 28){
						reactnum[mu]--;
						int attType = ATTACH; 
						int reactType = NDC_1; 
						AttachDettach(attType, outp_att, l, mds.h_r, mds.h_type, reactnum, mu, reactType, time, mds.h_connector, mds.h_att_l, mds.h_length, mds.h_link_l);
					}
					if (mu == 29){
						reactnum[mu]--;
						int attType = ATTACH;
						int reactType = NDC_2; 
						AttachDettach(attType, outp_att, l, mds.h_r, mds.h_type, reactnum, mu, reactType, time, mds.h_connector, mds.h_att_l, mds.h_length, mds.h_link_l);
					}
					if (mu == 30){
						reactnum[mu]--;
						int attType = ATTACH; 
						int reactType = NDC_3;
						AttachDettach(attType, outp_att, l, mds.h_r, mds.h_type, reactnum, mu, reactType, time, mds.h_connector, mds.h_att_l, mds.h_length, mds.h_link_l);
					}
					if (mu == 31){
						reactnum[mu]--;
						int attType = ATTACH;
						int reactType = NDC_4; 
						AttachDettach(attType, outp_att, l, mds.h_r, mds.h_type, reactnum, mu, reactType, time, mds.h_connector, mds.h_att_l, mds.h_length, mds.h_link_l);
					}
					if (mu == 32){
						reactnum[mu]--;
						int attType = ATTACH; 
						int reactType = NDC_5;
						AttachDettach(attType, outp_att, l, mds.h_r, mds.h_type, reactnum, mu, reactType, time, mds.h_connector, mds.h_att_l, mds.h_length, mds.h_link_l);
					}
					if (mu == 33){
						reactnum[mu]--;
						int attType = ATTACH; 
						int reactType = NDC_6;
						AttachDettach(attType, outp_att, l, mds.h_r, mds.h_type, reactnum, mu, reactType, time, mds.h_connector, mds.h_att_l, mds.h_length, mds.h_link_l);
					}
					if (mu == 34){
						reactnum[mu]--;
						int attType = ATTACH; 
						int reactType = NDC_7;
						AttachDettach(attType, outp_att, l, mds.h_r, mds.h_type, reactnum, mu, reactType, time, mds.h_connector, mds.h_att_l, mds.h_length, mds.h_link_l);
					}
				}
			} 
		}
	}
}