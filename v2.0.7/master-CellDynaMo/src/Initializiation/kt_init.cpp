/*
 * chromosome initialization module
 *
 *  Created on: Mar 15, 2018
 *	Edited on:	Mar 30, 2019
 *      Author: kliuchnikov
 */

#include "kt_init.h"

int kinIndex1;
int kinIndex2;

void KT_seeds(KT* kt, float3* r, vector <kin> kinetochore, int* type){
	int i, h;
	for (h = 0; h < KIN_COUNT; h++){
		kt[h].n = kinetochore.size()/KIN_COUNT;
		kt[h].r = (float3*)malloc(kt[h].n*sizeof(float3));
		kt[h].state = (int*)malloc(kt[h].n*sizeof(int));
		kt[h].type = (int*)malloc(kt[h].n*sizeof(int));
		kt[h].cell = (int*)malloc(kt[h].n*sizeof(int));
		kt[h].react = (int*)malloc(kt[h].n*sizeof(int));
		kt[h].rand = (int*)malloc(kt[h].n*sizeof(int));
		kt[h].cos = (float*)malloc(kt[h].n*sizeof(float));
		kt[h].radius = (float*)malloc(kt[h].n*sizeof(float));	
		kt[h].attach = (int*)malloc(kt[h].n*sizeof(int));
		for (i = 0; i < kt[h].n; i++){
			r[POLE_COUNT*MT_NUM*MAX_MT_LENGTH + h*kt[h].n + i].x = kinetochore[i + h*kt[h].n].x;
			r[POLE_COUNT*MT_NUM*MAX_MT_LENGTH + h*kt[h].n + i].y = kinetochore[i + h*kt[h].n].y;
			r[POLE_COUNT*MT_NUM*MAX_MT_LENGTH + h*kt[h].n + i].z = kinetochore[i + h*kt[h].n].z;
			int x_num = CORONA_X;
			int y_num = floor(x_num/RATIO);
			if (h % 2 == 0){
				if (i < x_num*y_num){
					type[POLE_COUNT*MT_NUM*MAX_MT_LENGTH + h*kt[h].n + i] = LEFT_NDC;
					kt[h].state[i] = 0;
				}
				else	{
					type[POLE_COUNT*MT_NUM*MAX_MT_LENGTH + h*kt[h].n + i] = SHELL_LEFT;
					kt[h].state[i] = 0;
				}
			}
			else	{
				if (i < x_num*y_num){
					type[POLE_COUNT*MT_NUM*MAX_MT_LENGTH + h*kt[h].n + i] = RIGHT_NDC;
					kt[h].state[i] = 1;
				}
				else	{
					type[POLE_COUNT*MT_NUM*MAX_MT_LENGTH + h*kt[h].n + i] = SHELL_RIGHT;
					kt[h].state[i] = 0;
				}
			}	
			kt[h].attach[i] = DETACH;
			kt[h].type[i] = kinetochore[i + h*kt[h].n].w;
			kt[h].cos[i] = kinetochore[i + h*kt[h].n].c;
			kt[h].radius[i] = kinetochore[i + h*kt[h].n].rad;
		}
		kt[h].r = &r[POLE_COUNT*MT_NUM*MAX_MT_LENGTH + h*kt[h].n];
	}
}

void KT_bonds(KT* kt, int* harmonicKinCount[KIN_COUNT], int* harmonicKin[KIN_COUNT], float* harmonicKinRadii[KIN_COUNT]){
	int h, i, j;
	int x_num = CORONA_X;
	int y_num = floor(x_num/RATIO);
	for (h = 0; h < KIN_COUNT; h++){
		for (i = 0; i < x_num*y_num; i++){
			for (j = 0; j < x_num*y_num; j++){
				if (i != j){
					if (myDistance(kt[h].r[i], kt[h].r[j]) <= 500){
						harmonicKinCount[h][i]++;
						harmonicKin[h][maxHarmonicKinPerMonomer*i + harmonicKinCount[h][i] - 1] = j;
						harmonicKinRadii[h][maxHarmonicKinPerMonomer*i + harmonicKinCount[h][i] - 1] = myDistance(kt[h].r[i], kt[h].r[j]);
					}		
				}
			}
		}
		for (i = x_num*y_num; i < kt[h].n - 1; i++){
			for (j = x_num*y_num; j < kt[h].n - 1; j++){
				if (i != j){
					if (myDistance(kt[h].r[i], kt[h].r[j]) <= 200 && myDistance(kt[h].r[i], kt[h].r[j]) > 100){
						harmonicKinCount[h][i]++;
						harmonicKin[h][maxHarmonicKinPerMonomer*i + harmonicKinCount[h][i] - 1] = j;
						harmonicKinRadii[h][maxHarmonicKinPerMonomer*i + harmonicKinCount[h][i] - 1] = myDistance(kt[h].r[i], kt[h].r[j]);
					}		
				}
			}
		}
	}
}

void KT_rotation(KT* kt, Config& config){
	int i, j, h;
	int chrom_n = KIN_COUNT/2;
	for (j = 0; j < chrom_n; j++){
		for (h = 0; h < 2; h++){
			for (i = 0; i < kt[h + j*2].n; i++){
				float3 k;
				//X-axis rotation
				k.x = kt[h + j*2].r[i].x;
				k.y = kt[h + j*2].r[i].y*cos(config.x_angle[j]) - kt[h + j*2].r[i].z*sin(config.x_angle[j]);
				k.z = kt[h + j*2].r[i].y*sin(config.x_angle[j]) + kt[h + j*2].r[i].z*cos(config.x_angle[j]);

				kt[h + j*2].r[i] = k;
				//Y-axis rotation
				k.x = kt[h + j*2].r[i].x*cos(config.y_angle[j]) + kt[h + j*2].r[i].z*sin(config.y_angle[j]);
				k.y = kt[h + j*2].r[i].y;
				k.z = -kt[h + j*2].r[i].x*sin(config.y_angle[j]) + kt[h + j*2].r[i].z*cos(config.y_angle[j]);

				kt[h + j*2].r[i] = k;
				//Z-axis rotation
				k.x = kt[h + j*2].r[i].x*cos(config.z_angle[j]) - kt[h + j*2].r[i].y*sin(config.z_angle[j]);
				k.y = kt[h + j*2].r[i].x*sin(config.z_angle[j]) + kt[h + j*2].r[i].y*cos(config.z_angle[j]);
				k.z = kt[h + j*2].r[i].z;

				kt[h + j*2].r[i] = k;
				//3D movement
				kt[h + j*2].r[i].x += config.kinetochore_move_x[j];
				kt[h + j*2].r[i].y += config.kinetochore_move_y[j];
				kt[h + j*2].r[i].z += config.kinetochore_move_z[j];
			}
		}
	}
	for (i = 0; i < KIN_COUNT; i++){
		for (j = 0; j < KIN_COUNT; j++){
			if (i != j){
				float dr = myDistance(kt[i].r[kn - 1], kt[j].r[kn - 1]);
				if (dr < 2*KIN_RADIUS - MT_RADIUS){
					printf("IMPOSSIBLE KINETOCHORE POSITION (KINETOCHORE#%d INTERSECT KINETOCHORE#%d)\n", i, j);
					exit(0);
				} 
			}
		}
	}
}

void KT_rotation_RAND(KT* kt, Config& config){
	int i, j, h;
	int chrom_n = KIN_COUNT/2;
	for (j = 0; j < chrom_n; j++){
		float rand_x_ang = -M_PI + ran2(&rseed)*M_PI*2;
		float rand_y_ang = -M_PI + ran2(&rseed)*M_PI*2;
		float rand_z_ang = -M_PI + ran2(&rseed)*M_PI*2;
		float rand_x = -(0.5*POLE_DISTANCE) + ran2(&rseed)*2*(0.5*POLE_DISTANCE);
		float rand_y = -(0.5*POLE_DISTANCE) + ran2(&rseed)*2*(0.5*POLE_DISTANCE);
		float rand_z = -(0.5*POLE_DISTANCE) + ran2(&rseed)*2*(0.5*POLE_DISTANCE);
		for (h = 0; h < 2; h++){
			for (i = 0; i < kt[h + j*2].n; i++){
				float3 k;
				//X-axis rotation
				k.x = kt[h + j*2].r[i].x;
				k.y = kt[h + j*2].r[i].y*cos(rand_x_ang) - kt[h + j*2].r[i].z*sin(rand_x_ang);
				k.z = kt[h + j*2].r[i].y*sin(rand_x_ang) + kt[h + j*2].r[i].z*cos(rand_x_ang);

				kt[h + j*2].r[i] = k;
				//Y-axis rotation
				k.x = kt[h + j*2].r[i].x*cos(rand_y_ang) + kt[h + j*2].r[i].z*sin(rand_y_ang);
				k.y = kt[h + j*2].r[i].y;
				k.z = -kt[h + j*2].r[i].x*sin(rand_y_ang) + kt[h + j*2].r[i].z*cos(rand_y_ang);

				kt[h + j*2].r[i] = k;
				//Z-axis rotation
				k.x = kt[h + j*2].r[i].x*cos(rand_z_ang) - kt[h + j*2].r[i].y*sin(rand_z_ang);
				k.y = kt[h + j*2].r[i].x*sin(rand_z_ang) + kt[h + j*2].r[i].y*cos(rand_z_ang);
				k.z = kt[h + j*2].r[i].z;

				kt[h + j*2].r[i] = k;
				//3D movement
				kt[h + j*2].r[i].x += rand_x;
				kt[h + j*2].r[i].y += rand_y;
				kt[h + j*2].r[i].z += rand_z;
			}
		}
	}
	/*for (i = 0; i < KIN_COUNT; i++){
		for (j = 0; j < KIN_COUNT; j++){
			if (i != j){
				float dr = myDistance(kt[i].r[kn - 1], kt[j].r[kn - 1]);
				if (dr < 2*KIN_RADIUS - MT_RADIUS){
					printf("IMPOSSIBLE KINETOCHORE POSITION (KINETOCHORE#%d INTERSECT KINETOCHORE#%d)\n", i, j);
					exit(0);
				} 
			}
		}
	}*/
}

void KT_rotation_RING(KT* kt, Config& config){
	int i, j, h;
	int chrom_n = KIN_COUNT/2;
	for (j = 0; j < chrom_n; j++){	
		float c1 = ran2(&rseed);
		float c2 = ran2(&rseed);
		float c3 = ran2(&rseed);
		float c = ran2(&rseed);
		for (h = 0; h < 2; h++){
			for (i = 0; i < kt[h + j*2].n; i++){
				float3 k;
				//Z-axis rotation


				k.x = kt[h + j*2].r[i].x*cos(c1*M_PI*2) - kt[h + j*2].r[i].y*sin(c1*M_PI*2);
				k.y = kt[h + j*2].r[i].x*sin(c1*M_PI*2) + kt[h + j*2].r[i].y*cos(c1*M_PI*2);
				k.z = kt[h + j*2].r[i].z;

				kt[h + j*2].r[i] = k;
				
				//X-axis rotation
				k.x = kt[h + j*2].r[i].x;
				k.y = kt[h + j*2].r[i].y*cos(c2*M_PI*2) - kt[h + j*2].r[i].z*sin(c2*M_PI*2);
				k.z = kt[h + j*2].r[i].y*sin(c2*M_PI*2) + kt[h + j*2].r[i].z*cos(c2*M_PI*2);

				kt[h + j*2].r[i] = k;
				//Y-axis rotation
				k.x = kt[h + j*2].r[i].x*cos(c3*M_PI*2) + kt[h + j*2].r[i].z*sin(c3*M_PI*2);
				k.y = kt[h + j*2].r[i].y;
				k.z = -kt[h + j*2].r[i].x*sin(c3*M_PI*2) + kt[h + j*2].r[i].z*cos(c3*M_PI*2);

				kt[h + j*2].r[i] = k;
				
				//3D movement
				kt[h + j*2].r[i].x += -250 + 500*c;
				kt[h + j*2].r[i].y += RING_RAD*cos(2*j*M_PI/chrom_n);
				kt[h + j*2].r[i].z += RING_RAD*sin(2*j*M_PI/chrom_n);
                //printf("%f\n", RING_RAD);
			}
		}
	}
	for (i = 0; i < KIN_COUNT; i++){
		for (j = 0; j < KIN_COUNT; j++){
			if (i != j){
				float dr = myDistance(kt[i].r[kn - 1], kt[j].r[kn - 1]);
				if (dr < 2*KIN_RADIUS - MT_RADIUS){
					printf("IMPOSSIBLE KINETOCHORE POSITION (KINETOCHORE#%d INTERSECT KINETOCHORE#%d)\n", i, j);
					//exit(0);
				} 
			}
		}
	}
}

vector<kin> generateKinetochoreHollow(){
	float xc = 0.0;
	float yc = 0.0;
	float zc = 0.0;

	kin center;
	center.x = xc;
	center.y = yc;
	center.z = zc;
	center.w = NO_NDC;

	vector <kin> r;
	int chrom_bn = 2*floor((CONTOUR_LENGTH/2 - KIN_RADIUS)/(2*KIN_RADIUS) + 0.5);
	kin ch[chrom_bn];
	kin c1, c2;
	c1.x = c2.x = c1.y = c2.y = c1.z = c2.z = 0.0;
	if (KIN_COUNT != 0){	
		for (int h = 0; h < KIN_COUNT/2; h ++){
			float theta;
			float phi;
			float big_r;
			if (BETA >= 0.925){
				big_r = 5000.0;
			}
			else	if (BETA >= 0 && BETA < 0.925){
				big_r = KIN_RADIUS/(1.0 - BETA);
			}
			else	if (BETA < 0 || BETA > 1){
				printf("INCORRECT BETA PARAMETER\n");
				exit(0);
			}

			int x_num = CORONA_X;
			int y_num = floor(x_num/RATIO);
			float y_delta = KIN_RADIUS/(y_num - 1);
			float alpha = (x_num - 1)*y_delta/(2*M_PI*big_r);
			float gamma = alpha*2*M_PI;
			float d_gamma = gamma/(x_num - 1);
//LEFT KINETOCHORE
			if (KT_SHAPE == KT_RECT){
				for (float y = -KIN_RADIUS*0.5; y <= (y_delta + KIN_RADIUS)*0.5 ; y += y_delta){	
					for (phi = -gamma*0.5; phi <= (gamma + d_gamma)*0.5; phi += d_gamma){
						kin coord;
						coord.x = -(KIN_RADIUS - (big_r - KIN_RADIUS)) - (big_r)*cos(phi);
						coord.y = y;
						coord.z = (big_r)*sin(phi);
						coord.w = NDC_0;

						float3 dr12, dr32;
						dr12.x = coord.x - (-KIN_RADIUS);
						dr12.y = coord.y - 0.0;
						dr12.z = coord.z - 0.0;
						float d_rad = sqrtf(dr12.x*dr12.x + dr12.y*dr12.y + dr12.z*dr12.z);

						dr32.x = KIN_RADIUS - (-KIN_RADIUS);
						dr32.y = 0.0 - 0.0;
						dr32.z = 0.0 - 0.0;

						float r12inv = 1.0f*pow(dr12.x*dr12.x + dr12.y*dr12.y + dr12.z*dr12.z, -0.5);
						float r32inv = 1.0f*pow(dr32.x*dr32.x + dr32.y*dr32.y + dr32.z*dr32.z, -0.5);
						float costheta = (dr12.x*dr32.x + dr12.y*dr32.y + dr12.z*dr32.z)*r12inv*r32inv;
						if (costheta > 1.0f){
							costheta = 1.0f;
						}
						else	if (costheta < -1.0f){
							costheta = -1.0f;
						}
						coord.c = costheta;
						coord.c = acos(coord.c);
						coord.rad = d_rad;

						dr12.x = coord.x - (-KIN_RADIUS);
						dr12.y = coord.y - 0.0;
						dr12.z = coord.z - 0.0;

						dr32.x = -1.5*KIN_RADIUS - (-KIN_RADIUS);
						dr32.y = 2*KIN_RADIUS - 0.0;
						dr32.z = 0.0 - 0.0;

						r12inv = 1.0f*pow(dr12.x*dr12.x + dr12.y*dr12.y + dr12.z*dr12.z, -0.5);
						r32inv = 1.0f*pow(dr32.x*dr32.x + dr32.y*dr32.y + dr32.z*dr32.z, -0.5);
						costheta = (dr12.x*dr32.x + dr12.y*dr32.y + dr12.z*dr32.z)*r12inv*r32inv;
						if (costheta > 1.0f){
							costheta = 1.0f;
						}
						else	if (costheta < -1.0f){
							costheta = -1.0f;
						}
						r.push_back(coord);
						if (coord.x > -(KIN_RADIUS + NO_ATT_ZONE)){
							r.pop_back();
						}
					}
				}
                if (ARMOR == 1){
				    for (theta = 0; theta < 2*M_PI; theta += 0.025*2*M_PI){
					    for (phi = 0; phi < 2*M_PI; phi += 0.025*2*M_PI){
						    kin coord;
						    coord.x = -(2*KIN_RADIUS - 80.0)*sin(theta)*cos(phi);
						    coord.y = (2*KIN_RADIUS - 80.0)*sin(theta)*sin(phi);
						    coord.z = (2*KIN_RADIUS - 80.0)*cos(theta);
						    coord.w = NO_NDC;

						    float3 dr12, dr32;
						    dr12.x = coord.x - (-KIN_RADIUS);
						    dr12.y = coord.y - 0.0;
						    dr12.z = coord.z - 0.0;
						    float d_rad = sqrtf(dr12.x*dr12.x + dr12.y*dr12.y + dr12.z*dr12.z);

						    dr32.x = KIN_RADIUS - (-KIN_RADIUS);
						    dr32.y = 0.0 - 0.0;
						    dr32.z = 0.0 - 0.0;

						    float r12inv = 1.0f*pow(dr12.x*dr12.x + dr12.y*dr12.y + dr12.z*dr12.z, -0.5);
						    float r32inv = 1.0f*pow(dr32.x*dr32.x + dr32.y*dr32.y + dr32.z*dr32.z, -0.5);
						    float costheta = (dr12.x*dr32.x + dr12.y*dr32.y + dr12.z*dr32.z)*r12inv*r32inv;
						    if (costheta > 1.0f){
							    costheta = 1.0f;
						    }
						    else	if (costheta < -1.0f){
							    costheta = -1.0f;
						    }
						    coord.c = costheta;
						    coord.c = acos(coord.c);
						    coord.rad = d_rad;

						    dr12.x = coord.x - (-KIN_RADIUS);
						    dr12.y = coord.y - 0.0;
						    dr12.z = coord.z - 0.0;

						    dr32.x = -1.5*KIN_RADIUS - (-KIN_RADIUS);
						    dr32.y = 2*KIN_RADIUS - 0.0;
						    dr32.z = 0.0 - 0.0;

						    r12inv = 1.0f*pow(dr12.x*dr12.x + dr12.y*dr12.y + dr12.z*dr12.z, -0.5);
						    r32inv = 1.0f*pow(dr32.x*dr32.x + dr32.y*dr32.y + dr32.z*dr32.z, -0.5);
						    costheta = (dr12.x*dr32.x + dr12.y*dr32.y + dr12.z*dr32.z)*r12inv*r32inv;
						    if (costheta > 1.0f){
							    costheta = 1.0f;
						    }
						    else	if (costheta < -1.0f){
							    costheta = -1.0f;
						    }
						    r.push_back(coord);
						    
						    if ((coord.x < -(KIN_RADIUS + NO_ATT_ZONE + 60.0) || coord.x > (KIN_RADIUS + NO_ATT_ZONE + 60.0)) && (coord.y > -KIN_RADIUS*0.5 - 60.0) && (coord.y <= (y_delta + KIN_RADIUS)*0.5 + 60.0)){
							    r.pop_back();
						    }
					    }
				    }
                }
			} 
			c1.x = -KIN_RADIUS;
			c1.y = 0;
			c1.z = 0;
			c1.w = NO_NDC;
			c1.c = 0.0;
			c1.rad = 0.0;
			r.push_back(c1);
			kinIndex1 = r.size() - 1;

			kt[0].centerID = kinIndex1;

			if (CHROM_ARMS == 1){
				int cn;
				for (cn = 0; cn < chrom_bn/2; cn++){
					ch[cn].x = -KIN_RADIUS;
					ch[cn].y = 2*KIN_RADIUS*(cn + 1);
					ch[cn].z = 0;
					ch[cn].w = NO_NDC;
				}
				for (cn = chrom_bn/2; cn < chrom_bn; cn++){
					ch[cn].x = -KIN_RADIUS;
					ch[cn].y = -2*KIN_RADIUS*(cn - chrom_bn/2 + 1);
					ch[cn].z = 0;
					ch[cn].w = NO_NDC;
				}

				float3 dr12, dr32;
				dr12.x = ch[0].x - c1.x;
				dr12.y = ch[0].y - c1.y;
				dr12.z = ch[0].z - c1.z;
				float d_rad = sqrtf(dr12.x*dr12.x + dr12.y*dr12.y + dr12.z*dr12.z);

				dr32.x = ch[chrom_bn/2].x - c1.x;
				dr32.y = ch[chrom_bn/2].y - c1.y;
				dr32.z = ch[chrom_bn/2].z - c1.z;

				float r12inv = 1.0f*pow(dr12.x*dr12.x + dr12.y*dr12.y + dr12.z*dr12.z, -0.5);
				float r32inv = 1.0f*pow(dr32.x*dr32.x + dr32.y*dr32.y + dr32.z*dr32.z, -0.5);
				float costheta = (dr12.x*dr32.x + dr12.y*dr32.y + dr12.z*dr32.z)*r12inv*r32inv;
				if (costheta > 1.0f){
					costheta = 1.0f;
				}
				else	if (costheta < -1.0f){
					costheta = -1.0f;
				}
				ch[0].c = costheta;
				ch[0].c = acos(ch[0].c);
				ch[0].rad = d_rad;
				ch[chrom_bn/2].c = costheta;
				ch[chrom_bn/2].c = acos(ch[chrom_bn/2].c);
				ch[chrom_bn/2].rad = d_rad;

				dr12.x = c1.x - ch[0].x;
				dr12.y = c1.y - ch[0].y;
				dr12.z = c1.z - ch[0].z;

				dr32.x = ch[1].x - ch[0].x;
				dr32.y = ch[1].y - ch[0].y;
				dr32.z = ch[1].z - ch[0].z;

				d_rad = sqrtf(dr32.x*dr32.x + dr32.y*dr32.y + dr32.z*dr32.z);

				r12inv = 1.0f*pow(dr12.x*dr12.x + dr12.y*dr12.y + dr12.z*dr12.z, -0.5);
				r32inv = 1.0f*pow(dr32.x*dr32.x + dr32.y*dr32.y + dr32.z*dr32.z, -0.5);
				costheta = (dr12.x*dr32.x + dr12.y*dr32.y + dr12.z*dr32.z)*r12inv*r32inv;
				if (costheta > 1.0f){
					costheta = 1.0f;
				}
				else	if (costheta < -1.0f){
					costheta = -1.0f;
				}
				ch[1].c = costheta;
				ch[1].c = acos(ch[1].c);
				ch[1].rad = d_rad;

				dr12.x = c1.x - ch[chrom_bn/2].x;
				dr12.y = c1.y - ch[chrom_bn/2].y;
				dr12.z = c1.z - ch[chrom_bn/2].z;

				dr32.x = ch[chrom_bn/2 + 1].x - ch[chrom_bn/2].x;
				dr32.y = ch[chrom_bn/2 + 1].y - ch[chrom_bn/2].y;
				dr32.z = ch[chrom_bn/2 + 1].z - ch[chrom_bn/2].z;

				d_rad = sqrtf(dr32.x*dr32.x + dr32.y*dr32.y + dr32.z*dr32.z);

				r12inv = 1.0f*pow(dr12.x*dr12.x + dr12.y*dr12.y + dr12.z*dr12.z, -0.5);
				r32inv = 1.0f*pow(dr32.x*dr32.x + dr32.y*dr32.y + dr32.z*dr32.z, -0.5);
				costheta = (dr12.x*dr32.x + dr12.y*dr32.y + dr12.z*dr32.z)*r12inv*r32inv;
				if (costheta > 1.0f){
					costheta = 1.0f;
				}
				else	if (costheta < -1.0f){
					costheta = -1.0f;
				}
				ch[chrom_bn/2 + 1].c = costheta;
				ch[chrom_bn/2 + 1].c = acos(ch[chrom_bn/2].c);
				ch[chrom_bn/2 + 1].rad = d_rad;
				for (cn = 2; cn < chrom_bn/2; cn++){
					dr12.x = ch[cn - 2].x - ch[cn - 1].x;
					dr12.y = ch[cn - 2].y - ch[cn - 1].y;
					dr12.z = ch[cn - 2].z - ch[cn - 1].z;

					dr32.x = ch[cn].x - ch[cn - 1].x;
					dr32.y = ch[cn].y - ch[cn - 1].y;
					dr32.z = ch[cn].z - ch[cn - 1].z;

					d_rad = sqrtf(dr32.x*dr32.x + dr32.y*dr32.y + dr32.z*dr32.z);

					r12inv = 1.0f*pow(dr12.x*dr12.x + dr12.y*dr12.y + dr12.z*dr12.z, -0.5);
					r32inv = 1.0f*pow(dr32.x*dr32.x + dr32.y*dr32.y + dr32.z*dr32.z, -0.5);
					costheta = (dr12.x*dr32.x + dr12.y*dr32.y + dr12.z*dr32.z)*r12inv*r32inv;
					if (costheta > 1.0f){
						costheta = 1.0f;
					}
					else	if (costheta < -1.0f){
						costheta = -1.0f;
					}
					ch[cn].c = costheta;
					ch[cn].c = acos(ch[cn].c);
					ch[cn].rad = d_rad;
				}
				for (cn = chrom_bn/2 + 2; cn < chrom_bn; cn++){
					dr12.x = ch[cn - 2].x - ch[cn - 1].x;
					dr12.y = ch[cn - 2].y - ch[cn - 1].y;
					dr12.z = ch[cn - 2].z - ch[cn - 1].z;

					dr32.x = ch[cn].x - ch[cn - 1].x;
					dr32.y = ch[cn].y - ch[cn - 1].y;
					dr32.z = ch[cn].z - ch[cn - 1].z;

					d_rad = sqrtf(dr32.x*dr32.x + dr32.y*dr32.y + dr32.z*dr32.z);

					r12inv = 1.0f*pow(dr12.x*dr12.x + dr12.y*dr12.y + dr12.z*dr12.z, -0.5);
					r32inv = 1.0f*pow(dr32.x*dr32.x + dr32.y*dr32.y + dr32.z*dr32.z, -0.5);
					costheta = (dr12.x*dr32.x + dr12.y*dr32.y + dr12.z*dr32.z)*r12inv*r32inv;
					if (costheta > 1.0f){
						costheta = 1.0f;
					}
					else	if (costheta < -1.0f){
						costheta = -1.0f;
					}
					ch[cn].c = costheta;
					ch[cn].c = acos(ch[cn].c);
					ch[cn].rad = d_rad;
				}
				for (cn = 0; cn < chrom_bn; cn++){
					r.push_back(ch[cn]);
				}
			}
//RIGHT KINETOCHORE
			if (KT_SHAPE == KT_RECT){
				for (float y = -KIN_RADIUS*0.5; y <= (y_delta + KIN_RADIUS)*0.5 ; y += y_delta){	
					for (phi = -gamma*0.5; phi <= (gamma + d_gamma)*0.5; phi += d_gamma){
						kin coord;
						coord.x = (KIN_RADIUS - (big_r - KIN_RADIUS)) + (big_r)*cos(phi);
						coord.y = y;
						coord.z = (big_r)*sin(phi);
						coord.w = NDC_0;

						float3 dr12, dr32;
						dr12.x = coord.x - KIN_RADIUS;
						dr12.y = coord.y - 0.0;
						dr12.z = coord.z - 0.0;
						float d_rad = sqrtf(dr12.x*dr12.x + dr12.y*dr12.y + dr12.z*dr12.z);

						dr32.x = -KIN_RADIUS - KIN_RADIUS;
						dr32.y = 0.0 - 0.0;
						dr32.z = 0.0 - 0.0;

						float r12inv = 1.0f*pow(dr12.x*dr12.x + dr12.y*dr12.y + dr12.z*dr12.z, -0.5);
						float r32inv = 1.0f*pow(dr32.x*dr32.x + dr32.y*dr32.y + dr32.z*dr32.z, -0.5);
						float costheta = (dr12.x*dr32.x + dr12.y*dr32.y + dr12.z*dr32.z)*r12inv*r32inv;
						if (costheta > 1.0f){
							costheta = 1.0f;
						}
						else	if (costheta < -1.0f){
							costheta = -1.0f;
						}
						coord.c = costheta;
						coord.c = acos(coord.c);
						coord.rad = d_rad;

						dr12.x = coord.x - KIN_RADIUS;
						dr12.y = coord.y - 0.0;
						dr12.z = coord.z - 0.0;

						dr32.x = 1.5*KIN_RADIUS - KIN_RADIUS;
						dr32.y = 2*KIN_RADIUS - 0.0;
						dr32.z = 0.0 - 0.0;

						r12inv = 1.0f*pow(dr12.x*dr12.x + dr12.y*dr12.y + dr12.z*dr12.z, -0.5);
						r32inv = 1.0f*pow(dr32.x*dr32.x + dr32.y*dr32.y + dr32.z*dr32.z, -0.5);
						costheta = (dr12.x*dr32.x + dr12.y*dr32.y + dr12.z*dr32.z)*r12inv*r32inv;
						if (costheta > 1.0f){
							costheta = 1.0f;
						}
						else	if (costheta < -1.0f){
							costheta = -1.0f;
						}
						r.push_back(coord);
						if (coord.x < KIN_RADIUS + NO_ATT_ZONE){
							r.pop_back();
						}
					}
				}
                if (ARMOR == 1){
			        for (theta = 0; theta < 2*M_PI; theta += 0.025*2*M_PI){
				        for (phi = 0; phi < 2*M_PI; phi += 0.025*2*M_PI){
					        kin coord;
					        coord.x = -(2*KIN_RADIUS - 80.0)*sin(theta)*cos(phi);
					        coord.y = (2*KIN_RADIUS - 80.0)*sin(theta)*sin(phi);
					        coord.z = (2*KIN_RADIUS - 80.0)*cos(theta);
					        coord.w = NO_NDC;
					        
					        float3 dr12, dr32;
					        dr12.x = coord.x - KIN_RADIUS;
					        dr12.y = coord.y - 0.0;
					        dr12.z = coord.z - 0.0;
					        float d_rad = sqrtf(dr12.x*dr12.x + dr12.y*dr12.y + dr12.z*dr12.z);

					        dr32.x = -KIN_RADIUS - KIN_RADIUS;
					        dr32.y = 0.0 - 0.0;
					        dr32.z = 0.0 - 0.0;

					        float r12inv = 1.0f*pow(dr12.x*dr12.x + dr12.y*dr12.y + dr12.z*dr12.z, -0.5);
					        float r32inv = 1.0f*pow(dr32.x*dr32.x + dr32.y*dr32.y + dr32.z*dr32.z, -0.5);
					        float costheta = (dr12.x*dr32.x + dr12.y*dr32.y + dr12.z*dr32.z)*r12inv*r32inv;
					        if (costheta > 1.0f){
						        costheta = 1.0f;
					        }
					        else	if (costheta < -1.0f){
						        costheta = -1.0f;
					        }
					        coord.c = costheta;
					        coord.c = acos(coord.c);
					        coord.rad = d_rad;

					        dr12.x = coord.x - KIN_RADIUS;
					        dr12.y = coord.y - 0.0;
					        dr12.z = coord.z - 0.0;

					        dr32.x = 1.5*KIN_RADIUS - KIN_RADIUS;
					        dr32.y = 2*KIN_RADIUS - 0.0;
					        dr32.z = 0.0 - 0.0;

					        r12inv = 1.0f*pow(dr12.x*dr12.x + dr12.y*dr12.y + dr12.z*dr12.z, -0.5);
					        r32inv = 1.0f*pow(dr32.x*dr32.x + dr32.y*dr32.y + dr32.z*dr32.z, -0.5);
					        costheta = (dr12.x*dr32.x + dr12.y*dr32.y + dr12.z*dr32.z)*r12inv*r32inv;
					        if (costheta > 1.0f){
						        costheta = 1.0f;
					        }
					        else	if (costheta < -1.0f){
						        costheta = -1.0f;
					        }
					        r.push_back(coord);
					        
					        if ((coord.x < -(KIN_RADIUS + NO_ATT_ZONE + 60.0) || coord.x > (KIN_RADIUS + NO_ATT_ZONE + 60.0)) && (coord.y > -KIN_RADIUS*0.5 - 60.0) && (coord.y <= (y_delta + KIN_RADIUS)*0.5 + 60.0)){
						        r.pop_back();
					        }
				        }
			        }
                }
			}
			c2.z = 0; 
			c2.y = 0;
			c2.x = KIN_RADIUS;
			c2.w = NO_NDC;
			c2.c = 0.0;
			c2.rad = 0.0;
			float kx = c2.x - c1.x;
			float ky = c2.y - c1.y;
			float kz = c2.z - c1.z;

			rk = sqrt(kx*kx + ky*ky + kz*kz);
			r.push_back(c2);
			kinIndex2 = r.size() - 1;	

			kt[1].centerID = kinIndex2;

			if (CHROM_ARMS == 1){
				int cn;
				for (cn = 0; cn < chrom_bn/2; cn++){
					ch[cn].x = KIN_RADIUS;
					ch[cn].y = 2*KIN_RADIUS*(cn + 1);
					ch[cn].z = 0;
					ch[cn].w = NO_NDC;
				}
				for (cn = chrom_bn/2; cn < chrom_bn; cn++){
					ch[cn].x = KIN_RADIUS;
					ch[cn].y = -2*KIN_RADIUS*(cn - chrom_bn/2 + 1);
					ch[cn].z = 0;
					ch[cn].w = NO_NDC;
				}
				
				float3 dr12, dr32;
				dr12.x = ch[0].x - c2.x;
				dr12.y = ch[0].y - c2.y;
				dr12.z = ch[0].z - c2.z;
				float d_rad = sqrtf(dr12.x*dr12.x + dr12.y*dr12.y + dr12.z*dr12.z);

				dr32.x = ch[chrom_bn/2].x - c2.x;
				dr32.y = ch[chrom_bn/2].y - c2.y;
				dr32.z = ch[chrom_bn/2].z - c2.z;

				float r12inv = 1.0f*pow(dr12.x*dr12.x + dr12.y*dr12.y + dr12.z*dr12.z, -0.5);
				float r32inv = 1.0f*pow(dr32.x*dr32.x + dr32.y*dr32.y + dr32.z*dr32.z, -0.5);
				float costheta = (dr12.x*dr32.x + dr12.y*dr32.y + dr12.z*dr32.z)*r12inv*r32inv;
				if (costheta > 1.0f){
					costheta = 1.0f;
				}
				else	if (costheta < -1.0f){
					costheta = -1.0f;
				}
				ch[0].c = costheta;
				ch[0].c = acos(ch[0].c);
				ch[0].rad = d_rad;
				ch[chrom_bn/2].c = costheta;
				ch[chrom_bn/2].c = acos(ch[chrom_bn/2].c);
				ch[chrom_bn/2].rad = d_rad;

				dr12.x = c2.x - ch[0].x;
				dr12.y = c2.y - ch[0].y;
				dr12.z = c2.z - ch[0].z;

				dr32.x = ch[1].x - ch[0].x;
				dr32.y = ch[1].y - ch[0].y;
				dr32.z = ch[1].z - ch[0].z;

				d_rad = sqrtf(dr32.x*dr32.x + dr32.y*dr32.y + dr32.z*dr32.z);

				r12inv = 1.0f*pow(dr12.x*dr12.x + dr12.y*dr12.y + dr12.z*dr12.z, -0.5);
				r32inv = 1.0f*pow(dr32.x*dr32.x + dr32.y*dr32.y + dr32.z*dr32.z, -0.5);
				costheta = (dr12.x*dr32.x + dr12.y*dr32.y + dr12.z*dr32.z)*r12inv*r32inv;
				if (costheta > 1.0f){
					costheta = 1.0f;
				}
				else	if (costheta < -1.0f){
					costheta = -1.0f;
				}
				ch[1].c = costheta;
				ch[1].c = acos(ch[1].c);
				ch[1].rad = d_rad;

				dr12.x = c2.x - ch[chrom_bn/2].x;
				dr12.y = c2.y - ch[chrom_bn/2].y;
				dr12.z = c2.z - ch[chrom_bn/2].z;

				dr32.x = ch[chrom_bn/2 + 1].x - ch[chrom_bn/2].x;
				dr32.y = ch[chrom_bn/2 + 1].y - ch[chrom_bn/2].y;
				dr32.z = ch[chrom_bn/2 + 1].z - ch[chrom_bn/2].z;

				d_rad = sqrtf(dr32.x*dr32.x + dr32.y*dr32.y + dr32.z*dr32.z);

				r12inv = 1.0f*pow(dr12.x*dr12.x + dr12.y*dr12.y + dr12.z*dr12.z, -0.5);
				r32inv = 1.0f*pow(dr32.x*dr32.x + dr32.y*dr32.y + dr32.z*dr32.z, -0.5);
				costheta = (dr12.x*dr32.x + dr12.y*dr32.y + dr12.z*dr32.z)*r12inv*r32inv;
				if (costheta > 1.0f){
					costheta = 1.0f;
				}
				else	if (costheta < -1.0f){
					costheta = -1.0f;
				}
				ch[chrom_bn/2 + 1].c = costheta;
				ch[chrom_bn/2 + 1].c = acos(ch[chrom_bn/2 + 1].c);
				ch[chrom_bn/2 + 1].rad = d_rad;
				for (cn = 2; cn < chrom_bn/2; cn++){
					dr12.x = ch[cn - 2].x - ch[cn - 1].x;
					dr12.y = ch[cn - 2].y - ch[cn - 1].y;
					dr12.z = ch[cn - 2].z - ch[cn - 1].z;

					dr32.x = ch[cn].x - ch[cn - 1].x;
					dr32.y = ch[cn].y - ch[cn - 1].y;
					dr32.z = ch[cn].z - ch[cn - 1].z;

					d_rad = sqrtf(dr32.x*dr32.x + dr32.y*dr32.y + dr32.z*dr32.z);

					r12inv = 1.0f*pow(dr12.x*dr12.x + dr12.y*dr12.y + dr12.z*dr12.z, -0.5);
					r32inv = 1.0f*pow(dr32.x*dr32.x + dr32.y*dr32.y + dr32.z*dr32.z, -0.5);
					costheta = (dr12.x*dr32.x + dr12.y*dr32.y + dr12.z*dr32.z)*r12inv*r32inv;
					if (costheta > 1.0f){
						costheta = 1.0f;
					}
					else	if (costheta < -1.0f){
						costheta = -1.0f;
					}
					ch[cn].c = costheta;
					ch[cn].c = acos(ch[cn].c);
					ch[cn].rad = d_rad;
				}
				for (cn = chrom_bn/2 + 2; cn < chrom_bn; cn++){
					dr12.x = ch[cn - 2].x - ch[cn - 1].x;
					dr12.y = ch[cn - 2].y - ch[cn - 1].y;
					dr12.z = ch[cn - 2].z - ch[cn - 1].z;

					dr32.x = ch[cn].x - ch[cn - 1].x;
					dr32.y = ch[cn].y - ch[cn - 1].y;
					dr32.z = ch[cn].z - ch[cn - 1].z;

					d_rad = sqrtf(dr32.x*dr32.x + dr32.y*dr32.y + dr32.z*dr32.z);

					r12inv = 1.0f*pow(dr12.x*dr12.x + dr12.y*dr12.y + dr12.z*dr12.z, -0.5);
					r32inv = 1.0f*pow(dr32.x*dr32.x + dr32.y*dr32.y + dr32.z*dr32.z, -0.5);
					costheta = (dr12.x*dr32.x + dr12.y*dr32.y + dr12.z*dr32.z)*r12inv*r32inv;
					if (costheta > 1.0f){
						costheta = 1.0f;
					}
					else	if (costheta < -1.0f){
						costheta = -1.0f;
					}
					ch[cn].c = costheta;
					ch[cn].c = acos(ch[cn].c);
					ch[cn].rad = d_rad;
				}
				for (cn = 0; cn < chrom_bn; cn++){
					r.push_back(ch[cn]);
				}
			}
			//r.push_back(center);	
		}
	}
	return r;
}
