/*
 *attachment-detachment module
 *
 *  Created on: Mar 15, 2018
 *	Edited on:	Mar 30, 2019
 *      Author: kliuchnikov
 */

#include "attach-detach.h"

void AttachDettach(int attType, FILE* outp_att, int l, float3* r, int* type, int reactnum[35], int mu, int reactType, float time, int* connector, float* att_l, float* len, float* link_l){
	int h, m;
	int countMT = 0;
	int rseed_kin = rseed + 123;
	int n_chrom = 2*floor((CONTOUR_LENGTH/2 - KIN_RADIUS)/(2*KIN_RADIUS) + 0.5);
	for (h = 0; h < POLE_COUNT; h++){	
		for (m = 0; m < MT_NUM; m++){
			int i = h*MT_NUM*MAX_MT_LENGTH + m*MAX_MT_LENGTH + MAX_MT_LENGTH - 1;
			float3 ri = r[i];
			int3 mtindex = get_cell(COORtransfer1(ri));
			if (mtindex.x > num_x_cell || mtindex.y > num_y_cell || mtindex.z > num_z_cell || mtindex.x <= 0 || mtindex.y <= 0 && mtindex.z <= 0){
				printf("att: %d\t%d\t%d\n", i, h, m);
			}	
			int s = COORintoARRAY(mtindex.x, mtindex.y, mtindex.z);
			if (s == l && attType == DETACH && type[i] != PLUS_DET_INVALID && pol[h][m].NDCn < NDC_PER_MT){
				countMT++;
				pol[h][m].rand = countMT;
			}
			else	if (s == l && pol[h][m].NDCn > 0 && attType == ATTACH){
				countMT++;
				pol[h][m].rand = countMT;
			}
		}
	}
	int x_num = CORONA_X;
	int y_num = floor(x_num/RATIO);
	float c = countMT*ran2(&rseed_kin);
	for (h = 0; h < POLE_COUNT; h++){	
		for (m = 0; m < MT_NUM; m++){
			int i = h*MT_NUM*MAX_MT_LENGTH + m*MAX_MT_LENGTH + MAX_MT_LENGTH - 1;
			if (type[i] != PLUS_ATT && pol[h][m].NDCn > 0){
				pol[h][m].NDCn = 0;
				connector[i] = 0;
				float minR = 1000.0;
				int k_id, n_id;
				for (int j = 0; j < KIN_COUNT; j++){ 
					for (int ks = 0; ks < x_num*y_num; ks++){
						if (kt[j].cell[ks] == l && kt[j].attach[ks] == attType && kt[j].type[ks] == reactType){
							float3 rndc = r[POLE_COUNT*MT_NUM*MAX_MT_LENGTH + j*kt[j].n + ks];
	                        if (myDistance(r[i], rndc) < minR){
		                        minR = myDistance(r[i], rndc);
		                        k_id = j;
		                        n_id = ks;
	                        } 
						}
					}
				}
				kt[k_id].attach[n_id] = (-1)*attType + 1;
			}
			float3 ri = r[i];
			int3 mtindex = get_cell(COORtransfer1(ri));
			/*if (mtindex.x > num_x_cell || mtindex.y > num_y_cell || mtindex.z > num_z_cell || mtindex.x <= 0 || mtindex.y <= 0 && mtindex.z <= 0){
				printf("att: %d\t%d\t%d\n", i, h, m);
			}*/
			int s = COORintoARRAY(mtindex.x, mtindex.y, mtindex.z);
			if (s == l && attType == DETACH && type[i] != PLUS_DET_INVALID && pol[h][m].NDCn < NDC_PER_MT){
				if (pol[h][m].rand - 1 < c && c <= pol[h][m].rand){
					int j1;
					int j2;
					int countNDC = 0;
					//printf("ATTACHMENT: %d\n", i);
					if (pol[h][m].NDCn == 0){
							att_l[i] = len[i];
					}
					type[i] = PLUS_ATT;
					pol[h][m].state = MT_STATE_DEPOL;
					reactnum[mu]++;
					float minR = 1000.0;
					for (int j = 0; j < KIN_COUNT; j++){
						int k = POLE_COUNT*MT_NUM*MAX_MT_LENGTH + (j + 1)*kn - n_chrom - 1;
						float3 rj = r[k];
						float dxK = rj.x - ri.x;
						float dyK = rj.y - ri.y;
						float dzK = rj.z - ri.z;

						float drK = sqrtf(dxK*dxK + dyK*dyK + dzK*dzK);
						if (j % 2 == 0){
							j1 = POLE_COUNT*MT_NUM*MAX_MT_LENGTH + (j + 1)*kn - n_chrom - 1;
							j2 = POLE_COUNT*MT_NUM*MAX_MT_LENGTH + (j + 2)*kn - n_chrom - 1;
						}
						else	{
							j1 = POLE_COUNT*MT_NUM*MAX_MT_LENGTH + j*kn - n_chrom - 1;
							j2 = POLE_COUNT*MT_NUM*MAX_MT_LENGTH + (j + 1)*kn - n_chrom - 1;
						}
					
						float3 rj1 = r[j1];
				        float dx1 = rj1.x - ri.x;
				        float dy1 = rj1.y - ri.y;
				        float dz1 = rj1.z - ri.z;

				        float dr1 = sqrtf(dx1*dx1 + dy1*dy1 + dz1*dz1);

				        float3 rj2 = r[j2];
				        float dx2 = rj2.x - ri.x;
				        float dy2 = rj2.y - ri.y;
				        float dz2 = rj2.z - ri.z;

				        float dr2 = sqrtf(dx2*dx2 + dy2*dy2 + dz2*dz2);
				        if (dr1 < dr2 && dr1 < minR){
				        	minR = dr1;
				        	connector[i] = j1;
				        }
				        else	if (dr2 <= dr1 && dr2 < minR){
				        	minR = dr2;
				        	connector[i] = j2;
				        }
						for (int ks = 0; ks < x_num*y_num; ks++){
							if (kt[j].cell[ks] == l && kt[j].attach[ks] == attType && kt[j].type[ks] == reactType){
								countNDC++;
								kt[j].rand[ks] = countNDC;
							}
						}
				    }
				    int triger = ATTACH;
			        att_out(time, outp_att, triger, i, connector[i]);
				    if (pol[h][m].NDCn == 0){
				    	link_l[i] = minR;
				    }
				    //printf("att: %d\t%f\t%d\t%f\t%f\n", i, link_l[i], connector[i]);   
					pol[h][m].NDCn++;
					float c1 = countNDC*ran2(&rseed_kin);
					for (int j = 0; j < KIN_COUNT; j++){ 
						for (int ks = 0; ks < x_num*y_num; ks++){
							if (kt[j].cell[ks] == l && kt[j].attach[ks] == attType && kt[j].type[ks] == reactType){							
								if (kt[j].rand[ks] - 1 < c1 && c1 <= kt[j].rand[ks]){
									kt[j].attach[ks] = (-1)*attType + 1;
								}
							}
						}
					}
				}
			}	
			else	if (s == l && pol[h][m].NDCn > 0 && attType == ATTACH){
				if (pol[h][m].rand - 1 < c && c <= pol[h][m].rand){
					if (pol[h][m].NDCn == 1){
						int triger = DETACH;
				        att_out(time, outp_att, triger, i, connector[i]);
						type[i] = PLUS_DET;
						connector[i] = 0;
						link_l[i] = 0.0;
					}
					pol[h][m].NDCn--;
					reactnum[mu]++;
					int countNDC = 0;
					for (int j = 0; j < KIN_COUNT; j++){ 
						for (int ks = 0; ks < x_num*y_num; ks++){
							if (kt[j].cell[ks] == l && kt[j].attach[ks] == attType && kt[j].type[ks] == reactType){
								countNDC++;
								kt[j].rand[ks] = countNDC;
							}
						}
					}
					float c1 = countNDC*ran2(&rseed_kin);
				    //}
					for (int j = 0; j < KIN_COUNT; j++){
						for (int ks = 0; ks < x_num*y_num; ks++){
							if (kt[j].cell[ks] == l && kt[j].attach[ks] == attType && kt[j].type[ks] == reactType){
								if (kt[j].rand[ks] - 1 < c1 && c1 <= kt[j].rand[ks]){
									kt[j].attach[ks] = (-1)*attType + 1;
								}
							}
						}
					}
				}
			}
		}
	}
}