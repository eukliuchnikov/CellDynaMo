/*
 * implicit diffusion module
 *
 *  Created on: Mar 15, 2018
 *	Edited on:	Mar 30, 2019
 *      Author: kliuchnikov
 */

#include "diffusion.h"

void diffusion(float3* r0, int* cell_type, int* Anum, int* Bnum, int D, int* neighbornum, int* neighbors, int* Aneighbornum, int* Aneighbors, float3* SVcenter){
	int l, n;
	int i, j, k;
	int rseed_diffA = rseed*2;
	for (i = 1; i < num_x_cell + 1; i++){
		for (j = 1; j < num_y_cell + 1; j++){
			for (k = 1; k < num_z_cell + 1; k++){
				l = COORintoARRAY(i, j, k);
				if (cell_type[l] == DIFFUSION_SPACE || cell_type[l] == DIFFUSION_MIX_SPACE || cell_type[l] == AUA_SPACE){
					for (n = 0; n < Anum[l]; n++){
						float r = ran2(&rseed_diffA);
						float p_diff = D*tau/(sv_size*sv_size*sv_size);
						if (r < p_diff){
							float part = p_diff/neighbornum[l];
							int section = 0;
							while (section < neighbornum[l]){
								if (section*part <= r && r < (section + 1)*part){
									Anum[l]--;
									Anum[neighbors[6*l + section]]++;
								}
								section++;
							}
						}
					}
					for (n = 0; n < Bnum[l]; n++){
						float r = ran2(&rseed_diffA);
						float p_diff = D*tau/(sv_size*sv_size*sv_size);
						if (r < p_diff){
							float part = p_diff/Aneighbornum[l];
							int section = 0;
							while (section < Aneighbornum[l]){
								if (section*part <= r && r < (section + 1)*part){
									Bnum[l]--;
									Bnum[Aneighbors[6*l + section]]++;
								}
								section++;
							}
						}
					}
				}
				else	if ((cell_type[l] == INNER_SPACE || cell_type[l] == AUA_SPACE) && Anum[l] != 0){
					float min_dist = 10000.0;
					int min_l;
					for (int kt; kt < KIN_COUNT; kt += 2){
						int n_chrom = 2*floor((CONTOUR_LENGTH/2 - KIN_RADIUS)/(2*KIN_RADIUS) + 0.5);
						float3 kt1_center = r0[POLE_COUNT*MT_NUM*MAX_MT_LENGTH + (kt + 1)*kn - n_chrom - 1];
						float3 kt2_center = r0[POLE_COUNT*MT_NUM*MAX_MT_LENGTH + (kt + 2)*kn - n_chrom - 1];
						float3 kin_center;
						kin_center.x = (kt1_center.x + kt2_center.x)/2;
						kin_center.y = (kt1_center.y + kt2_center.y)/2;
						kin_center.z = (kt1_center.z + kt2_center.z)/2;

						float dist = myDistance(SVcenter[l], kin_center);
						if (dist < min_dist){
							min_dist = dist;
							int3 index = get_cell(COORtransfer1(kin_center));
							min_l = COORintoARRAY(index.x, index.y, index.z);
						}
					}
					Anum[min_l] += Anum[l];
					Anum[l] = 0;

					/*float r = ran2(&rseed_diffA);
					float p_diff = D*tau/(sv_size*sv_size*sv_size);
					float part = p_diff/neighbornum[l];
					int section = 0;
					while (section < neighbornum[l]){
						if (section*part <= r && r < (section + 1)*part){
							Anum[neighbors[6*l + section]] += Anum[l];
							Anum[l] = 0;
						}
						section++;
					}*/
				}
			}
		}
	}
}