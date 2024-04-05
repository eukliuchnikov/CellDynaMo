/*
 *MTs initiation module
 *
 *  Created on: Mar 15, 2018
 *	Edited on:	Mar 30, 2019
 *      Author: kliuchnikov
 */

#include "mt_init.h"

float3* pol_center;

void poles(float3* r, int N, int* type){
	pol_center = (float3*)malloc(2*sizeof(float3));
	pol_center = &r[N - 2];
	r[N - 2].x = -POLE_DISTANCE; 
	r[N - 2].y = 0; 
	r[N - 2].z = 0;
	type[N - 2] = LEFT_POLE;

	r[N - 1].x = POLE_DISTANCE;
	r[N - 1].y = 0;
	r[N - 1].z = 0;
	type[N - 1] = RIGHT_POLE;
}

void MT_seeds(MT* pol[], int* type, float* length){
	pol[0] = (MT*)calloc(MT_NUM, sizeof(MT));
	pol[1] = (MT*)calloc(MT_NUM, sizeof(MT));
	int j, m, mi;
	for (j = 0; j < POLE_COUNT; j++){
		for (m = 0; m < MT_NUM; m++){
			pol[j][m].r = &mds.h_r[j*MT_NUM*MAX_MT_LENGTH + MAX_MT_LENGTH*m];
			pol[j][m].length = (START_LENGTH*ran2(&rseed) + MIN_MT_LENGTH - 1)*2*MT_RADIUS;
			pol[j][m].state = MT_STATE_POL;
			pol[j][m].attach = DETACH;
			pol[j][m].NDCn = 0;
		}	
	} 
	for (j = 0; j < POLE_COUNT; j++){
		for (m = 0; m < MT_NUM; m++){
			for (mi = 0; mi < MAX_MT_LENGTH; mi++){
				int i = j*MT_NUM*MAX_MT_LENGTH + m*MAX_MT_LENGTH + mi;
				length[i] = pol[j][m].length;
				type[i] = MT_REG;
				if (i == j*MT_NUM*MAX_MT_LENGTH + m*MAX_MT_LENGTH + MAX_MT_LENGTH - 1){
					type[i] = PLUS_DET;
				}
				if (i == m*MAX_MT_LENGTH){
					type[i] = LEFT_MINUS;
				}
				if (i == MT_NUM*MAX_MT_LENGTH + m*MAX_MT_LENGTH){
					type[i] = RIGHT_MINUS;
				}
			}
		}
	}
	for (j = 0; j < POLE_COUNT; j++){
		for (m = 0; m < MT_NUM; m++){
			float theta = (1 - 2*j)*(M_PI/4 + 0.5*M_PI*ran2(&rseed));
			float phi = (1 - 2*j)*(0.5*M_PI*ran2(&rseed) - M_PI/4);
			pol[j][m].r[0].x = pol_center[j].x + POLE_RADIUS*sin(theta)*cos(phi);
			pol[j][m].r[0].y = pol_center[j].y + POLE_RADIUS*sin(theta)*sin(phi);
			pol[j][m].r[0].z = pol_center[j].z + POLE_RADIUS*cos(theta);
			pol[j][m].r[2].x = pol_center[j].x + (POLE_RADIUS + pol[j][m].length)*sin(theta)*cos(phi);
			pol[j][m].r[2].y = pol_center[j].y + (POLE_RADIUS + pol[j][m].length)*sin(theta)*sin(phi);
			pol[j][m].r[2].z = pol_center[j].z + (POLE_RADIUS + pol[j][m].length)*cos(theta);
			float mindist = 100.0;
			pol[j][m].r[1].x = pol_center[j].x + (POLE_RADIUS + pol[j][m].length/2)*sin(theta)*cos(phi);
			pol[j][m].r[1].y = pol_center[j].y + (POLE_RADIUS + pol[j][m].length/2)*sin(theta)*sin(phi);
			pol[j][m].r[1].z = pol_center[j].z + (POLE_RADIUS + pol[j][m].length/2)*cos(theta);
			int m1;
			/*for (m1 = 0; m1 < m; m1++){
				float dx = pol[j][m1].r[0].x - pol[j][m].r[0].x;
				float dy = pol[j][m1].r[0].y - pol[j][m].r[0].y;
				float dz = pol[j][m1].r[0].z - pol[j][m].r[0].z;
				float dr = sqrt(dx*dx + dy*dy + dz*dz);
				if (dr < mindist){
					mindist = dr;
				}
			}
			if (mindist < MT_RADIUS*2.0){
				m--;
			}*/		
		}
	}
}