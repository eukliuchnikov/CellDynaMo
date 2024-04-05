/*
 * implicit components update module
 *
 *  Created on: Mar 15, 2018
 *	Edited on:	Mar 30, 2019
 *      Author: kliuchnikov
 */

#include "implicit.h"

void membrane_impl(int* cell_type){
	int i, j, k, h;
	int3 cell_index;
	for (i = 1; i < num_x_cell + 1; i++){
		for (j = 1; j < num_y_cell + 1; j++){
			for (k = 1; k < num_z_cell + 1; k++){
				int3 cell;
				cell.x = i;
				cell.y = j;
				cell.z = k;
				float3 r;
				r = get_coord(cell);
				h = COORintoARRAY(i, j, k);
				if (MEM_SHAPE == MEM_SPHERE || MEM_SHAPE == MEM_ELIPSOID){
					if (pow((r.x - 0.5*float(XBOX)), 2)/pow(ELLIPSE_A, 2) + pow((r.y - 0.5*float(YBOX)), 2)/pow(ELLIPSE_B, 2) + pow((r.z - 0.5*float(ZBOX)), 2)/pow(ELLIPSE_C, 2) > 1){
						cell_type[h] = OUTER_SPACE;	//no diffusion (static)
					}
					else{
						cell_type[h] = INNER_SPACE;	//diffusion
					}
				}
				else	if (MEM_SHAPE == MEM_CUBE || MEM_SHAPE == MEM_CUBOID){
					cell_type[h] = INNER_SPACE;
				}				
			}
		}
	}
}

void mt_dynamic(int* cell_type, float3* r, int* mt_pol_conc, int* mt_depol_conc, int* mt_det_conc, int* mt_att_conc, int* type){
	int k, m, i, l, j, h;
	for (i = 1; i < num_x_cell + 1; i++){
		for (j = 1; j < num_y_cell + 1; j++){
			for (k = 1; k < num_z_cell + 1; k++){
				l = COORintoARRAY(i, j, k);
				mt_pol_conc[l] = 0;
				mt_depol_conc[l] = 0;
				mt_det_conc[l] = 0;
				mt_att_conc[l] = 0;
			}
		}
	}
	int x_num = CORONA_X;
	int y_num = floor(x_num/RATIO);
	for (k = 0; k < POLE_COUNT; k++){
		for (m = 0; m < MT_NUM; m++){
			int occup = 0;
			i = k*MT_NUM*MAX_MT_LENGTH + m*MAX_MT_LENGTH + MAX_MT_LENGTH - 1;
			float3 ri = r[i];
			int3 mtindex = get_cell(COORtransfer1(ri));
			l = COORintoARRAY(mtindex.x, mtindex.y, mtindex.z);
			if (type[i] == PLUS_DET && pol[k][m].NDCn == 0){
				for (j = 0; j < KIN_COUNT; j++){
					for (h = 0; h < x_num*y_num; h++){
						if (kt[j].cell[h] == l){
							float3 rj = r[POLE_COUNT*MT_NUM*MAX_MT_LENGTH + j*kt[j].n + h];
							if (myDistance(ri, rj) < NDC_LENGTH && occup == 0){
								mt_det_conc[l]++;
								occup = 1;
							}
						}
					}
				}
				if (occup == 0){
					type[i] = PLUS_DET_INVALID;
				}
			}
			else	if (type[i] == PLUS_DET_INVALID){
				for (j = 0; j < KIN_COUNT; j++){
					for (h = 0; h < x_num*y_num; h++){
						if (kt[j].cell[h] == l){
							float3 rj = r[POLE_COUNT*MT_NUM*MAX_MT_LENGTH + j*kt[j].n + h];
							if (myDistance(ri, rj) < NDC_LENGTH){
								type[i] = PLUS_DET;
							}
						}
					}
				}
			}
			else	if (type[i] == PLUS_ATT){
				mt_att_conc[l]++;
				float minR = 1000.0;
                float maxR = 0.0;
				for (j = 0; j < KIN_COUNT; j++){
					float3 rkin = r[POLE_COUNT*MT_NUM*MAX_MT_LENGTH + (j + 1)*kt[j].n - 1];

					for (h = 0; h < x_num*y_num; h++){
						float3 rj = r[POLE_COUNT*MT_NUM*MAX_MT_LENGTH + j*kt[j].n + h];
                        if (myDistance(ri, rj) < minR){
                        	float drc = myDistance(rkin, rj);
	                        minR = myDistance(ri, rj);
	                        maxR = drc;
                        } 
                    }
                }
                if (pol[k][m].NDCn < NDC_PER_MT){
					mt_det_conc[l]++;
                }
            }
			if (pol[k][m].state == MT_STATE_POL){
				mt_pol_conc[l]++;
			}
			else	if (pol[k][m].state == MT_STATE_DEPOL){
				mt_depol_conc[l]++;
			}
		}
	}
}

void pole_initial(int* cell_type, float3* r){
	int i, j, k, h, l;
	
	float3 max, min;
	max.x = POLE_DISTANCE + POLE_RADIUS;
	max.y = POLE_RADIUS;
	max.z = POLE_RADIUS;
	min.x = -max.x;
	min.y = -max.y;
	min.z = -max.z;

	float3 max_impl, min_impl;
	max_impl = COORtransfer1(max);
	min_impl = COORtransfer1(min);
	
	int3 max_impl_index = get_cell(max_impl);	
	int3 min_impl_index = get_cell(min_impl);

	int imin, imax, jmin, jmax, kmin, kmax;
	imin = min_impl_index.x;
	jmin = min_impl_index.y;
	kmin = min_impl_index.z;
	imax = max_impl_index.x;
	jmax = max_impl_index.y;
	kmax = max_impl_index.z;
	for (i = imin; i < imax + 1; i++){
		for (j = jmin; j < jmax + 1; j++){
			for (k = kmin; k < kmax + 1; k++){
				l = COORintoARRAY(i, j, k);
				int3 ijk;
                ijk.x = i;
                ijk.y = j;
                ijk.z = k;
                float3 xyz;
                xyz = COORtransfer(get_coord(ijk)); 
                if ((xyz.x - POLE_DISTANCE)*(xyz.x - POLE_DISTANCE) + (xyz.y)*(xyz.y) + (xyz.z)*(xyz.z) <= POLE_RADIUS*POLE_RADIUS || (xyz.x + POLE_DISTANCE)*(xyz.x + POLE_DISTANCE) + (xyz.y)*(xyz.y) + (xyz.z)*(xyz.z) <= POLE_RADIUS*POLE_RADIUS){
					cell_type[l] = OUTER_SPACE;	// no diffusion (static) 
				}
				else 	if ((xyz.x - POLE_DISTANCE)*(xyz.x - POLE_DISTANCE) + (xyz.y)*(xyz.y) + (xyz.z)*(xyz.z) <= 2000*2000 || (xyz.x + POLE_DISTANCE)*(xyz.x + POLE_DISTANCE) + (xyz.y)*(xyz.y) + (xyz.z)*(xyz.z) <= 2000*2000){
					cell_type[l] = AUA_SPACE;
				}
				else{
					cell_type[l] = INNER_SPACE;
				}
			}
		}
	}
}

void NDC80_dynamic(int* cell_type, float3* r, int* NDC80_0_conc, int* NDC80_1_conc, int* NDC80_2_conc, int* NDC80_3_conc, int* NDC80_4_conc, int* NDC80_5_conc, int* NDC80_6_conc, int* NDC80_7_conc, int* NDC80_0_ATT_conc, int* NDC80_1_ATT_conc, int* NDC80_2_ATT_conc, int* NDC80_3_ATT_conc, int* NDC80_4_ATT_conc, int* NDC80_5_ATT_conc, int* NDC80_6_ATT_conc, int* NDC80_7_ATT_conc){
	int i, j, k, h, m, p, l;

	for (i = 1; i < num_x_cell + 1; i++){
		for (j = 1; j < num_y_cell + 1; j++){
			for (k = 1; k < num_z_cell + 1; k++){
				p = COORintoARRAY(i, j, k);
				NDC80_0_conc[p] = 0;
				NDC80_1_conc[p] = 0;
				NDC80_2_conc[p] = 0;
				NDC80_3_conc[p] = 0;
				NDC80_4_conc[p] = 0;
				NDC80_5_conc[p] = 0;
				NDC80_6_conc[p] = 0;
				NDC80_7_conc[p] = 0;
				NDC80_0_ATT_conc[p] = 0;
				NDC80_1_ATT_conc[p] = 0;
				NDC80_2_ATT_conc[p] = 0;
				NDC80_3_ATT_conc[p] = 0;
				NDC80_4_ATT_conc[p] = 0;
				NDC80_5_ATT_conc[p] = 0;
				NDC80_6_ATT_conc[p] = 0;
				NDC80_7_ATT_conc[p] = 0;

				if (cell_type[p] == DIFFUSION_SPACE){
					cell_type[p] == INNER_SPACE;
				}
				if (cell_type[p] == DIFFUSION_MIX_SPACE){
					cell_type[p] == AUA_SPACE;
				}
			}
		}
	}
	int x_num = CORONA_X;
	int y_num = floor(x_num/RATIO);
	for (h = 0; h < KIN_COUNT; h++){	
		for (m = 0; m < x_num*y_num; m++){
			float3 rj = r[POLE_COUNT*MT_NUM*MAX_MT_LENGTH + h*kt[h].n + m];
			int3 kinindex = get_cell(COORtransfer1(rj));
			p = COORintoARRAY(kinindex.x, kinindex.y, kinindex.z);
			kt[h].cell[m] = p;
			if (kt[h].type[m] == NDC_0){
				NDC80_0_conc[p]++;
				if (kt[h].attach[m] == ATTACH){
					NDC80_0_conc[p]--;
					NDC80_0_ATT_conc[p]++;
				}
			}
			else	if (kt[h].type[m] == NDC_1){
				NDC80_1_conc[p]++;
				if (kt[h].attach[m] == ATTACH){
					NDC80_1_conc[p]--;
					NDC80_1_ATT_conc[p]++;
				}
			}
			else	if (kt[h].type[m] == NDC_2){
				NDC80_2_conc[p]++;
				if (kt[h].attach[m] == ATTACH){
					NDC80_2_conc[p]--;
					NDC80_2_ATT_conc[p]++;
				}
			}
			else	if (kt[h].type[m] == NDC_3){
				NDC80_3_conc[p]++;
				if (kt[h].attach[m] == ATTACH){
					NDC80_3_conc[p]--;
					NDC80_3_ATT_conc[p]++;
				}
			}
			else	if (kt[h].type[m] == NDC_4){
				NDC80_4_conc[p]++;
				if (kt[h].attach[m] == ATTACH){
					NDC80_4_conc[p]--;
					NDC80_4_ATT_conc[p]++;
				}
			}
			else	if (kt[h].type[m] == NDC_5){
				NDC80_5_conc[p]++;
				if (kt[h].attach[m] == ATTACH){
					NDC80_5_conc[p]--;
					NDC80_5_ATT_conc[p]++;
				}
			}
			else	if (kt[h].type[m] == NDC_6){
				NDC80_6_conc[p]++;
				if (kt[h].attach[m] == ATTACH){
					NDC80_6_conc[p]--;
					NDC80_6_ATT_conc[p]++;
				}
			}
			else	if (kt[h].type[m] == NDC_7){
				NDC80_7_conc[p]++;
				if (kt[h].attach[m] == ATTACH){
					NDC80_7_conc[p]--;
					NDC80_7_ATT_conc[p]++;
				}
			}
		}
		if (h % 2 == 0){
			int n_chrom = 2*floor((CONTOUR_LENGTH/2 - KIN_RADIUS)/(2*KIN_RADIUS) + 0.5);
			float3 kt1_center = r[POLE_COUNT*MT_NUM*MAX_MT_LENGTH + (h + 1)*kn - n_chrom - 1];
			float3 kt2_center = r[POLE_COUNT*MT_NUM*MAX_MT_LENGTH + (h + 2)*kn - n_chrom - 1];
			float3 kin_center;
			kin_center.x = (kt1_center.x + kt2_center.x)/2;
			kin_center.y = (kt1_center.y + kt2_center.y)/2;
			kin_center.z = (kt1_center.z + kt2_center.z)/2;
			for (i = 1; i < num_x_cell + 1; i++){
				for (j = 1; j < num_y_cell + 1; j++){
					for (k = 1; k < num_z_cell + 1; k++){
						int3 cell;
						cell.x = i;
						cell.y = j;
						cell.z = k;
						float3 r0;
						r0 = COORtransfer(get_coord(cell));
						l = COORintoARRAY(i, j, k);
						float Rad_allow = myDistance(kin_center, kt1_center) + 2*sv_size + KIN_RADIUS;
						if ((r0.x - kin_center.x)*(r0.x - kin_center.x) + (r0.y - kin_center.y)*(r0.y - kin_center.y) + (r0.z - kin_center.z)*(r0.z - kin_center.z) < Rad_allow*Rad_allow && cell_type[l] == INNER_SPACE){
							cell_type[l] = DIFFUSION_SPACE;
						}
						else	if ((r0.x - kin_center.x)*(r0.x - kin_center.x) + (r0.y - kin_center.y)*(r0.y - kin_center.y) + (r0.z - kin_center.z)*(r0.z - kin_center.z) < Rad_allow*Rad_allow && cell_type[l] == AUA_SPACE){
							cell_type[l] = DIFFUSION_MIX_SPACE;
						}
					}
				}
			}
		}
	}
}

void count_initial(float3* r0, int* cell_type, int* aurora_B_conc, int* aurora_A_conc, int N, int* number_AB){
	int i, j, k, h, n, p, b;
	int3 cell_index;
	p = 0;
	b = 0;
	for (i = 1; i < num_x_cell + 1; i++){
		for (j = 1; j < num_y_cell + 1; j++){
			for (k = 1; k < num_z_cell + 1; k++){
				int3 cell;
				cell.x = i;
				cell.y = j;
				cell.z = k;
				float3 r;
				r = get_coord(cell);
				h = COORintoARRAY(i, j, k);
				if ((pow((r.x - 0.5*float(XBOX)), 2)/pow(400, 2) <= 1) && (pow((r.y - 0.5*float(YBOX)), 2)/pow(400, 2) <= 1) && (pow((r.z - 0.5*float(ZBOX)), 2)/pow(400, 2) <= 1)){
					p++;
				}
				if ((pow((r.x - 0.5*float(XBOX)), 2)/pow(2000, 2) <= 1) && (pow((r.y - 0.5*float(YBOX)), 2)/pow(2000, 2) <= 1) && (pow((r.z - 0.5*float(ZBOX)), 2)/pow(2000, 2) <= 1)){
					b++;
					if ((pow((r.x - 0.5*float(XBOX)), 2)/pow(POLE_RADIUS, 2) > 1) && (pow((r.y - 0.5*float(YBOX)), 2)/pow(POLE_RADIUS, 2) > 1) && (pow((r.z - 0.5*float(ZBOX)), 2)/pow(POLE_RADIUS, 2) > 1)){
						b--;
					}
				}		
			}
		}
	}
	float subV = sv_size*sv_size*sv_size;
	float ApV = subV*AU_CONC*6.02*pow(10, -7);	
	int ABN = floor(ApV*p);
    int AAN;
    if (AURORA_A == 1){
	    AAN = floor(ApV*b);
    }
    else if (AURORA_A == 0){
	    AAN = 0;
    }
	//float volume = 4/3*M_PI*400*400*400;
	//int ABN = floor(AU_CONC*0.6*volume);

	for (int s = 0; s < POLE_COUNT; s++){
		float3 rpole = r0[N - s - 1];
		int n_aa = AAN;
		int number_AA = 0;		
		for (i = 1; i < num_x_cell + 1; i++){
			for (j = 1; j < num_y_cell + 1; j++){
				for (k = 1; k < num_z_cell + 1; k++){
					int3 cell;
					cell.x = i;
					cell.y = j;
					cell.z = k;
					float3 r;
					r = COORtransfer(get_coord(cell));
					h = COORintoARRAY(i, j, k);

					float rijk = sqrt(pow((r.x - rpole.x), 2) + pow((r.y - rpole.y), 2) + pow((r.z - rpole.z), 2));
					int portionAA = floor(10*AAN*exp(-rijk*rijk/(2*400*400))/(400*sqrt(2*M_PI)));
					if (cell_type[h] != OUTER_SPACE){			
						aurora_A_conc[h] += portionAA;
						number_AA += portionAA;
						//n_aa -= portionAA;
					}
				}
			}
		}
		printf("%d molecules of Aurora A kinase were placed around centrosome #%d\n", number_AA, s);
		printf("******************************************************************\n");
	}
    number_AB[0] = 0;
	for (int m = 0; m < KIN_COUNT; m += 2){
		int n_chrom = 2*floor((CONTOUR_LENGTH/2 - KIN_RADIUS)/(2*KIN_RADIUS) + 0.5);
		float3 rcom1 = r0[POLE_COUNT*MT_NUM*MAX_MT_LENGTH + (m + 1)*kt[0].n - n_chrom - 1];
		float3 rcom2 = r0[POLE_COUNT*MT_NUM*MAX_MT_LENGTH + (m + 2)*kt[0].n - n_chrom - 1];
        
        float3 kincenter;
        kincenter.x = (rcom1.x + rcom2.x)/2;
        kincenter.y = (rcom1.y + rcom2.y)/2;
        kincenter.z = (rcom1.z + rcom2.z)/2;

		int n_ab = ABN;
		int number = m + 1;
        int counter = number_AB[0];
		for (i = 1; i < num_x_cell + 1; i++){
			for (j = 1; j < num_y_cell + 1; j++){
				for (k = 1; k < num_z_cell + 1; k++){
					int3 cell;
					cell.x = i;
					cell.y = j;
					cell.z = k;
					float3 r;
					r = COORtransfer(get_coord(cell));
					h = COORintoARRAY(i, j, k);

					float rijk = sqrt(pow((r.x - kincenter.x), 2) + pow((r.y - kincenter.y), 2) + pow((r.z - kincenter.z), 2));
					//int portionAB = floor(ABN*exp(-rijk*rijk/(2*(2*KIN_RADIUS + sv_size)*(2*KIN_RADIUS + sv_size)))/((2*KIN_RADIUS + sv_size)*sqrt(2*M_PI)));
                    if (m == 0){
					    int portionAB = floor(30*n_ab*exp(-rijk*rijk/(2*350*350))/(350*sqrt(2*M_PI)));
					    if (cell_type[h] != OUTER_SPACE){
						    aurora_B_conc[h] += portionAB;
						    counter += portionAB;
					    }
                    }
                    else    {
                        
                        int portionAB = floor(30*n_ab*exp(-rijk*rijk/(2*350*350))/(350*sqrt(2*M_PI)));
					    if (cell_type[h] != OUTER_SPACE){
						    aurora_B_conc[h] += portionAB;
						    counter -= portionAB;
                        }
                    } 
				}
			}
		}
        if (m == 0){
            number_AB[0] = counter;
        }
        else if (counter != 0){
            int3 null = get_cell(COORtransfer1(kincenter));
			n = COORintoARRAY(null.x, null.y, null.z);
            aurora_B_conc[n] += counter;
            counter = 0;
        }
		printf("%d molecules of Aurora B kinase were placed in chromosome #%d\n", number_AB[0], number);
		printf("******************************************************************\n");
		/*printf("It remains to place %d particles of Aurora-B in the center of chromosome #%d\n", n_ab, number);
		printf("%d particles of Aurora-B were added in the center of chromosome #%d\n", n_ab, number);
		int3 null = get_cell(COORtransfer1(kincenter));
		n = COORintoARRAY(null.x, null.y, null.z);
		aurora_B_conc[n] += n_ab;
		n_ab -= n_ab;
		printf("******************************************************************\n");*/
	}
}

void count_dyn(float3* r0, int* cell_type, int* aurora_B_conc, int N, int* number_AB){
	int i, j, k, h, n, p, b;
	int3 cell_index;
	p = 0;
	b = 0;
	for (i = 1; i < num_x_cell + 1; i++){
		for (j = 1; j < num_y_cell + 1; j++){
			for (k = 1; k < num_z_cell + 1; k++){
				int3 cell;
				cell.x = i;
				cell.y = j;
				cell.z = k;
				float3 r;
				r = get_coord(cell);
				h = COORintoARRAY(i, j, k);
				aurora_B_conc[h] = 0;
				if ((pow((r.x - 0.5*float(XBOX)), 2)/pow(400, 2) <= 1) && (pow((r.y - 0.5*float(YBOX)), 2)/pow(400, 2) <= 1) && (pow((r.z - 0.5*float(ZBOX)), 2)/pow(400, 2) <= 1)){
					p++;
				}		
			}
		}
	}
	float subV = sv_size*sv_size*sv_size;
	float ApV = subV*AU_CONC*6.02*pow(10, -7);	
	int ABN = floor(ApV*p);
	for (int m = 0; m < KIN_COUNT; m += 2){
		int n_chrom = 2*floor((CONTOUR_LENGTH/2 - KIN_RADIUS)/(2*KIN_RADIUS) + 0.5);
		float3 rcom1 = r0[POLE_COUNT*MT_NUM*MAX_MT_LENGTH + (m + 1)*kt[0].n - n_chrom - 1];
		float3 rcom2 = r0[POLE_COUNT*MT_NUM*MAX_MT_LENGTH + (m + 2)*kt[0].n - n_chrom - 1];
        
        float3 kincenter;
        kincenter.x = (rcom1.x + rcom2.x)/2;
        kincenter.y = (rcom1.y + rcom2.y)/2;
        kincenter.z = (rcom1.z + rcom2.z)/2;

		int n_ab = ABN;
		int number = m + 1;
        int counter = number_AB;

		for (i = 1; i < num_x_cell + 1; i++){
			for (j = 1; j < num_y_cell + 1; j++){
				for (k = 1; k < num_z_cell + 1; k++){
					int3 cell;
					cell.x = i;
					cell.y = j;
					cell.z = k;
					float3 r;
					r = COORtransfer(get_coord(cell));
					h = COORintoARRAY(i, j, k);

					float rijk = sqrt(pow((r.x - kincenter.x), 2) + pow((r.y - kincenter.y), 2) + pow((r.z - kincenter.z), 2));
					//int portionAB = floor(ABN*exp(-rijk*rijk/(2*(2*KIN_RADIUS + sv_size)*(2*KIN_RADIUS + sv_size)))/((2*KIN_RADIUS + sv_size)*sqrt(2*M_PI)));
					int portionAB = floor(30*n_ab*exp(-rijk*rijk/(2*350*350))/(350*sqrt(2*M_PI)));
				    if (cell_type[h] != OUTER_SPACE){
                        /*if (rijk > KIN_RADIUS && portionAB > 0){
							float c = ran2(&rseed);
							int shift = int(10.0*c - 5);
							portionAB += shift;
							if (portionAB < 0){
								portionAB = 0;
							}
						}*/				    
                        aurora_B_conc[h] += portionAB;
					    counter -= portionAB;
                        
                    }
				}
			}
		}
        printf("COUNTER1: %d\n", counter);
		if (counter != 0){
            int3 null = get_cell(COORtransfer1(kincenter));
			n = COORintoARRAY(null.x, null.y, null.z);
            aurora_B_conc[n] += counter;
            counter = 0;
        }
        printf("COUNTER2: %d\n", counter);
	}
}
