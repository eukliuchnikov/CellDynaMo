/*
 *MT' shortening module
 *
 *  Created on: Mar 15, 2018
 *	Edited on:	Mar 30, 2019
 *      Author: kliuchnikov
 */

#include "mt_short.h"

void SHORTfunc(int l, float3* r, int* type, int & minimiz, float* length, float* att_l){
	int m, j, i, q;
	float p;
	int countMT = 0;
	int rseed_kin = rseed + 123;
	for (j = 0; j < POLE_COUNT; j++){
		for (m = 0; m < MT_NUM; m++){
			i = j*MT_NUM*MAX_MT_LENGTH + m*MAX_MT_LENGTH + MAX_MT_LENGTH - 1;
			float3 ri = r[i];
			int3 mtindex = get_cell(COORtransfer1(ri));
			if (mtindex.x > num_x_cell || mtindex.y > num_y_cell || mtindex.z > num_z_cell || mtindex.x <= 0 || mtindex.y <= 0 && mtindex.z <= 0){
				printf("short: %d\n", i);
			}
			q = COORintoARRAY(mtindex.x, mtindex.y, mtindex.z);
			if (l == q && pol[j][m].state == MT_STATE_DEPOL && pol[j][m].length > MIN_MT_LENGTH*4*MT_RADIUS){
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
			if (mtindex.x > num_x_cell || mtindex.y > num_y_cell || mtindex.z > num_z_cell || mtindex.x <= 0 || mtindex.y <= 0 && mtindex.z <= 0){
				printf("short: %d\n", i);
			}
			q = COORintoARRAY(mtindex.x, mtindex.y, mtindex.z);
			if (l == q && pol[j][m].state == MT_STATE_DEPOL && pol[j][m].length > MIN_MT_LENGTH*4*MT_RADIUS){
				if (pol[j][m].rand - 1 < c && c <= pol[j][m].rand){
					float dlx1 = pol[j][m].r[2].x - pol[j][m].r[1].x;
					float dly1 = pol[j][m].r[2].y - pol[j][m].r[1].y;
					float dlz1 = pol[j][m].r[2].z - pol[j][m].r[1].z;
					float dl1 = sqrt(dlx1*dlx1 + dly1*dly1 + dlz1*dlz1);
					float dlx2 = pol[j][m].r[1].x - pol[j][m].r[0].x;
					float dly2 = pol[j][m].r[1].y - pol[j][m].r[0].y;
					float dlz2 = pol[j][m].r[1].z - pol[j][m].r[0].z;
					float dl2 = sqrt(dlx2*dlx2 + dly2*dly2 + dlz2*dlz2);
					float delta = pol[j][m].length - (dl1 + dl2);

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
					float dx2 = pol[j][m].r[2].x - pol[j][m].r[1].x;
					float dy2 = pol[j][m].r[2].y - pol[j][m].r[1].y;
					float dz2 = pol[j][m].r[2].z - pol[j][m].r[1].z;
					float dr2 = sqrt(dx2*dx2 + dy2*dy2 + dz2*dz2);
					dx2 = dx2/dr2;
					dy2 = dy2/dr2;
					dz2 = dz2/dr2;

					if ((dl1 + dl2 > att_l[i] + MT_RADIUS || pol[j][m].length > att_l[i] + MT_RADIUS) && type[i] == PLUS_ATT){
						pol[j][m].length = att_l[i];
						length[i] = pol[j][m].length;
						length[i - 1] = pol[j][m].length;
						length[i - 2] = pol[j][m].length;
						pol[j][m].r[1].x = pol[j][m].r[0].x + (pol[j][m].length/2)*dx1;
						pol[j][m].r[1].y = pol[j][m].r[0].y + (pol[j][m].length/2)*dy1;
						pol[j][m].r[1].z = pol[j][m].r[0].z + (pol[j][m].length/2)*dz1;
						pol[j][m].r[2].x = pol[j][m].r[1].x + (pol[j][m].length/2)*dx2;
						pol[j][m].r[2].y = pol[j][m].r[1].y + (pol[j][m].length/2)*dy2;
						pol[j][m].r[2].z = pol[j][m].r[1].z + (pol[j][m].length/2)*dz2;
					}
					if (delta < -MT_RADIUS || dl1 - dl2 < -MT_RADIUS){
						pol[j][m].length = dl1 + dl2;
						length[i] = pol[j][m].length;
						length[i - 1] = pol[j][m].length;
						length[i - 2] = pol[j][m].length;
						pol[j][m].r[1].x = pol[j][m].r[0].x + (pol[j][m].length/2)*dx1;
						pol[j][m].r[1].y = pol[j][m].r[0].y + (pol[j][m].length/2)*dy1;
						pol[j][m].r[1].z = pol[j][m].r[0].z + (pol[j][m].length/2)*dz1;
					}
					pol[j][m].r[2].x -= 2.0*MT_RADIUS*dx2;
					pol[j][m].r[2].y -= 2.0*MT_RADIUS*dy2;
					pol[j][m].r[2].z -= 2.0*MT_RADIUS*dz2;
					pol[j][m].length -= 2.0*MT_RADIUS;
					pol[j][m].r[1].x -= MT_RADIUS*dx1;
					pol[j][m].r[1].y -= MT_RADIUS*dy1;
					pol[j][m].r[1].z -= MT_RADIUS*dz1;
					length[i] = pol[j][m].length;
					length[i - 1] = pol[j][m].length;
					length[i - 2] = pol[j][m].length;
					if(type[i] == PLUS_ATT){
						minimiz = 1;
					}
				}
			}
			else 	if (l == q && pol[j][m].state == MT_STATE_DEPOL && pol[j][m].length <= MIN_MT_LENGTH*4*MT_RADIUS){
				pol[j][m].length = MIN_MT_LENGTH*4*MT_RADIUS;
				float dx0 = pol[j][m].r[0].x - r[mds.N - 2 + j].x;
				float dy0 = pol[j][m].r[0].y - r[mds.N - 2 + j].y;
				float dz0 = pol[j][m].r[0].z - r[mds.N - 2 + j].z;
				float dr0 = sqrt(dx0*dx0 + dy0*dy0 + dz0*dz0);
				dx0 = dx0/dr0;
				dy0 = dy0/dr0;
				dz0 = dz0/dr0;
				pol[j][m].r[1].x = r[mds.N - 2 + j].x + (POLE_RADIUS + pol[j][m].length/2)*dx0;
				pol[j][m].r[1].y = r[mds.N - 2 + j].y + (POLE_RADIUS + pol[j][m].length/2)*dy0;
				pol[j][m].r[1].z = r[mds.N - 2 + j].z + (POLE_RADIUS + pol[j][m].length/2)*dz0;
				pol[j][m].r[2].x = r[mds.N - 2 + j].x + (POLE_RADIUS + pol[j][m].length)*dx0;
				pol[j][m].r[2].y = r[mds.N - 2 + j].y + (POLE_RADIUS + pol[j][m].length)*dy0;
				pol[j][m].r[2].z = r[mds.N - 2 + j].z + (POLE_RADIUS + pol[j][m].length)*dz0;
				length[i] = pol[j][m].length;
				length[i - 1] = pol[j][m].length;
				length[i - 2] = pol[j][m].length;
				if(type[i] == PLUS_ATT){
					minimiz = 1;
				}
			}
		}
	}
}