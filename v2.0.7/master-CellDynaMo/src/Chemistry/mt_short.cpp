/*
 * MT's shortening module
 *
 *  Created on: Mar 15, 2018
 *	Edited on:	Feb 20, 2022
 *      Author: kliuchnikov
 */

#include "mt_short.h"

void SHORTfunc(int l, float3* r, int* type, int & minimiz, float* length, float* att_l){
	int m, j, i, q, mi;
	float p;
	int countMT = 0;
	int rseed_kin = rseed + 123;
	for (j = 0; j < POLE_COUNT; j++){
		for (m = 0; m < MT_NUM; m++){
			i = j*MT_NUM*MAX_MT_LENGTH + m*MAX_MT_LENGTH + MAX_MT_LENGTH - 1;
			float3 ri = r[i];
			int3 mtindex = get_cell(COORtransfer1(ri));
			if (mtindex.x > num_x_cell || mtindex.y > num_y_cell || mtindex.z > num_z_cell || mtindex.x <= 0 || mtindex.y <= 0 || mtindex.z <= 0){
				printf("short: %d\n", i);
			}
			q = COORintoARRAY(mtindex.x, mtindex.y, mtindex.z);
			if (l == q && pol[j][m].state == MT_STATE_DEPOL && pol[j][m].length > MIN_MT_LENGTH){
				countMT++;
				pol[j][m].rand = countMT;
			}
		}
	}
	float c = countMT*ran2(&rseed_kin);
	for (j = 0; j < POLE_COUNT; j++){
		for (m = 0; m < MT_NUM; m++){
			i = j*MT_NUM*MAX_MT_LENGTH + m*MAX_MT_LENGTH + MAX_MT_LENGTH - 1;
			float3 ri = r[i];
			int3 mtindex = get_cell(COORtransfer1(ri));
			if (mtindex.x > num_x_cell || mtindex.y > num_y_cell || mtindex.z > num_z_cell || mtindex.x <= 0 || mtindex.y <= 0 || mtindex.z <= 0){
				printf("short: %d\n", i);
			}
			q = COORintoARRAY(mtindex.x, mtindex.y, mtindex.z);
			if (l == q && pol[j][m].state == MT_STATE_DEPOL && pol[j][m].length > MIN_MT_LENGTH){
				if (pol[j][m].rand - 1 < c && c <= pol[j][m].rand){
					float sum_l = 0.0;
					float tip_diff = 0.0;
					for (mi = 1; mi < MAX_MT_LENGTH; mi ++){
						float dx = pol[j][m].r[MAX_MT_LENGTH - mi].x - pol[j][m].r[MAX_MT_LENGTH - mi - 1].x;
						float dy = pol[j][m].r[MAX_MT_LENGTH - mi].y - pol[j][m].r[MAX_MT_LENGTH - mi - 1].y;
						float dz = pol[j][m].r[MAX_MT_LENGTH - mi].z - pol[j][m].r[MAX_MT_LENGTH - mi - 1].z;
						float dr = sqrt(dx*dx + dy*dy + dz*dz);
						if (mi == 1){
							tip_diff += dr;
						}
						if (mi == 2){
							tip_diff -= dr;
						}
						sum_l += dr;
					}
					float delta = pol[j][m].length - sum_l;

					float dx0 = pol[j][m].r[0].x - r[mds.N - 2 + j].x;
					float dy0 = pol[j][m].r[0].y - r[mds.N - 2 + j].y;
					float dz0 = pol[j][m].r[0].z - r[mds.N - 2 + j].z;
					float dr0 = sqrt(dx0*dx0 + dy0*dy0 + dz0*dz0);
					dx0 = dx0/dr0;
					dy0 = dy0/dr0;
					dz0 = dz0/dr0;
					float dx1 = pol[j][m].r[1].x - pol[j][m].r[0].x;
					float dy1 = pol[j][m].r[1].y - pol[j][m].r[0].y;
					float dz1 = pol[j][m].r[1].z - pol[j][m].r[0].z;
					float dr1 = sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1);
					dx1 = dx1/dr1;
					dy1 = dy1/dr1;
					dz1 = dz1/dr1;

					if ((sum_l > att_l[i] + MT_RADIUS || pol[j][m].length > att_l[i] + MT_RADIUS) && type[i] == PLUS_ATT){
						pol[j][m].length = att_l[i];

						for (mi = 1; mi < MAX_MT_LENGTH; mi ++){
							pol[j][m].r[MAX_MT_LENGTH - mi].x = pol[j][m].r[0].x + ((MAX_MT_LENGTH - mi)*pol[j][m].length/(MAX_MT_LENGTH - 1))*dx1;
							pol[j][m].r[MAX_MT_LENGTH - mi].y = pol[j][m].r[0].y + ((MAX_MT_LENGTH - mi)*pol[j][m].length/(MAX_MT_LENGTH - 1))*dy1;
							pol[j][m].r[MAX_MT_LENGTH - mi].z = pol[j][m].r[0].z + ((MAX_MT_LENGTH - mi)*pol[j][m].length/(MAX_MT_LENGTH - 1))*dz1;
						}
						length[i] = pol[j][m].length;
						for (mi = 1; mi < MAX_MT_LENGTH; mi ++){
							length[i - mi] = pol[j][m].length;
						}

					}
					if (delta < -MT_RADIUS || tip_diff < -MT_RADIUS){
						pol[j][m].length = sum_l;
						length[i] = pol[j][m].length;
						for (mi = 1; mi < MAX_MT_LENGTH; mi ++){
							length[i - mi] = pol[j][m].length;
						}
						for (mi = 2; mi < MAX_MT_LENGTH; mi ++){
							pol[j][m].r[MAX_MT_LENGTH - mi].x = pol[j][m].r[0].x + ((MAX_MT_LENGTH - mi)*pol[j][m].length/(MAX_MT_LENGTH - 1))*dx1;
							pol[j][m].r[MAX_MT_LENGTH - mi].y = pol[j][m].r[0].y + ((MAX_MT_LENGTH - mi)*pol[j][m].length/(MAX_MT_LENGTH - 1))*dy1;
							pol[j][m].r[MAX_MT_LENGTH - mi].z = pol[j][m].r[0].z + ((MAX_MT_LENGTH - mi)*pol[j][m].length/(MAX_MT_LENGTH - 1))*dz1;
						}
					}

					for (mi = 1; mi < MAX_MT_LENGTH; mi ++){
						float dx = pol[j][m].r[MAX_MT_LENGTH - mi].x - pol[j][m].r[MAX_MT_LENGTH - mi - 1].x;
						float dy = pol[j][m].r[MAX_MT_LENGTH - mi].y - pol[j][m].r[MAX_MT_LENGTH - mi - 1].y;
						float dz = pol[j][m].r[MAX_MT_LENGTH - mi].z - pol[j][m].r[MAX_MT_LENGTH - mi - 1].z;
						float dr = sqrt(dx*dx + dy*dy + dz*dz);
						dx = 2.0*MT_RADIUS*dx/dr;
						dy = 2.0*MT_RADIUS*dy/dr;
						dz = 2.0*MT_RADIUS*dz/dr;

						float dx0 = (MAX_MT_LENGTH - mi)*dx/(MAX_MT_LENGTH - 1);
						float dy0 = (MAX_MT_LENGTH - mi)*dy/(MAX_MT_LENGTH - 1);
						float dz0 = (MAX_MT_LENGTH - mi)*dz/(MAX_MT_LENGTH - 1);
						float dr0 = sqrt(dx0*dx0 + dy0*dy0 + dz0*dz0);

						pol[j][m].r[MAX_MT_LENGTH - mi].x -= dx0;
						pol[j][m].r[MAX_MT_LENGTH - mi].y -= dy0;
						pol[j][m].r[MAX_MT_LENGTH - mi].z -= dz0;
					}
					pol[j][m].length -= 2.0*MT_RADIUS;
					length[i] = pol[j][m].length;
					for (mi = 1; mi < MAX_MT_LENGTH; mi ++){
						length[i - mi] = pol[j][m].length;
					}
					if(type[i] == PLUS_ATT){
						minimiz = 1;
					}
				}
			}
			else 	if (l == q && pol[j][m].state == MT_STATE_DEPOL && pol[j][m].length <= MIN_MT_LENGTH){
				pol[j][m].length = MIN_MT_LENGTH;
				float dx0 = pol[j][m].r[0].x - r[mds.N - 2 + j].x;
				float dy0 = pol[j][m].r[0].y - r[mds.N - 2 + j].y;
				float dz0 = pol[j][m].r[0].z - r[mds.N - 2 + j].z;
				float dr0 = sqrt(dx0*dx0 + dy0*dy0 + dz0*dz0);
				dx0 = dx0/dr0;
				dy0 = dy0/dr0;
				dz0 = dz0/dr0;

				for (mi = 1; mi < MAX_MT_LENGTH; mi ++){
					pol[j][m].r[MAX_MT_LENGTH - mi].x = r[mds.N - 2 + j].x + (POLE_RADIUS + (MAX_MT_LENGTH - mi)*pol[j][m].length/(MAX_MT_LENGTH - 1))*dx0;
					pol[j][m].r[MAX_MT_LENGTH - mi].y = r[mds.N - 2 + j].y + (POLE_RADIUS + (MAX_MT_LENGTH - mi)*pol[j][m].length/(MAX_MT_LENGTH - 1))*dy0;
					pol[j][m].r[MAX_MT_LENGTH - mi].z = r[mds.N - 2 + j].z + (POLE_RADIUS + (MAX_MT_LENGTH - mi)*pol[j][m].length/(MAX_MT_LENGTH - 1))*dz0;
				}
				length[i] = pol[j][m].length;
				for (mi = 1; mi < MAX_MT_LENGTH; mi ++){
					length[i - mi] = pol[j][m].length;
				}
				if(type[i] == PLUS_ATT){
					minimiz = 1;
				}
			}
		}
	}
    minimiz = 1;
}
