/*
 *MT' catastrophe module
 *
 *  Created on: Mar 15, 2018
 *	Edited on:	Mar 30, 2019
 *      Author: kliuchnikov
 */

#include "mt_cat.h"

void CATfunc(int l, float3* r){
	int m, j, i, q;
	int countMT = 0;
	int rseed_kin = rseed + 123;
	for (j = 0; j < POLE_COUNT; j++){
		for (m = 0; m < MT_NUM; m++){
			i = j*MT_NUM*MAX_MT_LENGTH + m*MAX_MT_LENGTH + MAX_MT_LENGTH - 1;
			float3 ri = r[i];
			int3 mtindex = get_cell(COORtransfer1(ri));
			if (mtindex.x > num_x_cell || mtindex.y > num_y_cell || mtindex.z > num_z_cell || mtindex.x <= 0 || mtindex.y <= 0 && mtindex.z <= 0){
				printf("cat: %d\n", i);
			}	
			q = COORintoARRAY(mtindex.x, mtindex.y, mtindex.z);
			if (l == q && pol[j][m].state == MT_STATE_POL){
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
				printf("cat: %d\n", i);
			}	
			q = COORintoARRAY(mtindex.x, mtindex.y, mtindex.z);
			if (l == q && pol[j][m].state == MT_STATE_POL){
				if (pol[j][m].rand - 1 < c && c <= pol[j][m].rand){
					pol[j][m].state = MT_STATE_DEPOL;
				}
			}
		}
	}
}