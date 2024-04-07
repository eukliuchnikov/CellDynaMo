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

void mt_dynamic(int* cell_type, float3* r, int* mt_pol_conc, int* mt_depol_conc, int* mt_det_conc, int* mt_att_conc, int* type, 
	int* connect_1, int* connect_2, float* len1, float* len2, float3* center, int* kinesin_att_conc, int* kinesin_det_conc, float* h_att_l){

	int k, m, i, l, j, h;
	for (i = 1; i < num_x_cell + 1; i++){
		for (j = 1; j < num_y_cell + 1; j++){
			for (k = 1; k < num_z_cell + 1; k++){
				l = COORintoARRAY(i, j, k);
				mt_pol_conc[l] = 0;
				mt_depol_conc[l] = 0;
				mt_det_conc[l] = 0;
				mt_att_conc[l] = 0;
				kinesin_att_conc[l] = 0;
				kinesin_det_conc[l] = 0;
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
			if (mtindex.x > num_x_cell || mtindex.y > num_y_cell || mtindex.z > num_z_cell || mtindex.x <= 0 || mtindex.y <= 0 || mtindex.z <= 0){
				printf("mt_dyn1: %d\n", i);
			}
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
	int n_chrom = 2*floor((CONTOUR_LENGTH/2 - KIN_RADIUS)/(2*KIN_RADIUS) + 0.5);
	for (k = 0; k < POLE_COUNT; k++){
		for (m = 0; m < MT_NUM; m++){
			if (pol[k][m].kinesin == 0){
				for (int mi = 0; mi < MAX_MT_LENGTH - 1; mi++){
					int i = k*MT_NUM*MAX_MT_LENGTH + m*MAX_MT_LENGTH + mi;
					int j = i + 1;

					float3 ri1 = r[i];  //a0
				    float3 rj1 = r[j];  //a1
		            
		            float3 intersect1, intersect2;
		            float min_dist = 10000.0;
		            float dx, dy, dz;

		            float3 p_vector1, p1_norm;
		            p_vector1.x = rj1.x - ri1.x;  //A
		            p_vector1.y = rj1.y - ri1.y;
		            p_vector1.z = rj1.z - ri1.z;

		            float dp1 = sqrtf(p_vector1.x*p_vector1.x + p_vector1.y*p_vector1.y + p_vector1.z*p_vector1.z); //magA
		            p1_norm.x = p_vector1.x/dp1; //_A
		            p1_norm.y = p_vector1.y/dp1;
		            p1_norm.z = p_vector1.z/dp1;

		            for (int p = 0; p < KIN_COUNT; p++){
		                for (int ch = 2; ch <= (n_chrom/2); ch++){
		                    int js = POLE_COUNT*MT_NUM*MAX_MT_LENGTH + (p + 1)*kn - ch;
		                    float3 ri2 = r[js]; //b0
				            float3 rj2 = r[js + 1];  //b1
		                    float3 p_vector2, p2_norm;
		                    p_vector2.x = rj2.x - ri2.x; //B
		                    p_vector2.y = rj2.y - ri2.y;
		                    p_vector2.z = rj2.z - ri2.z;

		                    float dp2 = sqrtf(p_vector2.x*p_vector2.x + p_vector2.y*p_vector2.y + p_vector2.z*p_vector2.z); //magB
		                    p2_norm.x = p_vector2.x/dp2; //_B
		                    p2_norm.y = p_vector2.y/dp2;
		                    p2_norm.z = p_vector2.z/dp2;

		                    float3 cross;
		                    cross.x = p1_norm.y*p2_norm.z - p2_norm.y*p1_norm.z; //cross
		                    cross.y = p1_norm.z*p2_norm.x - p2_norm.z*p1_norm.x;
		                    cross.z = p1_norm.x*p2_norm.y - p2_norm.x*p1_norm.y;

		                    float d_cross = sqrtf(cross.x*cross.x + cross.y*cross.y + cross.z*cross.z); //denom
		                    d_cross = d_cross*d_cross;
		                    if (d_cross == 0){
		                        float d0 = p1_norm.x*(ri2.x - ri1.x) + p1_norm.y*(ri2.y - ri1.y) + p1_norm.z*(ri2.z - ri1.z);
		                        float d1 = p1_norm.x*(rj2.x - ri1.x) + p1_norm.y*(rj2.y - ri1.y) + p1_norm.z*(rj2.z - ri1.z);
		                        if (d0 <= 0 && 0 >= d1){
		                            if (sqrtf(d0*d0) < sqrtf(d1*d1)){
		                                intersect1 = ri1;
		                                intersect2 = ri2;
		                                dx = intersect2.x - intersect1.x;
		                                dy = intersect2.y - intersect1.y;
		                                dz = intersect2.z - intersect1.z;
		                                min_dist = sqrtf(dx*dx + dy*dy + dz*dz);
		                            }
		                            else    {
		                                intersect1 = ri1;
		                                intersect2 = rj2;
		                                
		                                dx = intersect2.x - intersect1.x;
		                                dy = intersect2.y - intersect1.y;
		                                dz = intersect2.z - intersect1.z;
		                                min_dist = sqrtf(dx*dx + dy*dy + dz*dz);
		                            }
		                        }
		                        else if (d0 >= dp1 && dp1 <= d1){
		                            if (sqrtf(d0*d0) < sqrtf(d1*d1)){
		                                intersect1 = rj1;
		                                intersect2 = ri2;
		                                dx = intersect2.x - intersect1.x;
		                                dy = intersect2.y - intersect1.y;
		                                dz = intersect2.z - intersect1.z;
		                                min_dist = sqrtf(dx*dx + dy*dy + dz*dz);
		                            }
		                            else    {
		                                intersect1 = rj1;
		                                intersect2 = rj2;
		                                dx = intersect2.x - intersect1.x;
		                                dy = intersect2.y - intersect1.y;
		                                dz = intersect2.z - intersect1.z;
		                                min_dist = sqrtf(dx*dx + dy*dy + dz*dz);
		                            }
		                        }
		                    }
		                    else {
		                        float3 t;
		                        t.x = ri2.x - ri1.x;
		                        t.y = ri2.y - ri1.y;
		                        t.z = ri2.z - ri1.z;

		                        float det0 = t.x*(p2_norm.y*cross.z - cross.y*p2_norm.z) + t.y*(p2_norm.z*cross.x - cross.z*p2_norm.x) + t.z*(p2_norm.x*cross.y - cross.x*p2_norm.y);
		                        float det1 = t.x*(p1_norm.y*cross.z - cross.y*p1_norm.z) + t.y*(p1_norm.z*cross.x - cross.z*p1_norm.x) + t.z*(p1_norm.x*cross.y - cross.x*p1_norm.y);
		                        float t0 = det0/d_cross;
		                        float t1 = det1/d_cross;

		                        float3 p1, p2;
		                        p1.x = ri1.x + t0*p1_norm.x;
		                        p1.y = ri1.y + t0*p1_norm.y;
		                        p1.z = ri1.z + t0*p1_norm.z;

		                        p2.x = ri2.x + t1*p2_norm.x;
		                        p2.y = ri2.y + t1*p2_norm.y;
		                        p2.z = ri2.z + t1*p2_norm.z;

		                        if (t0 < 0){
		                            p1.x = ri1.x;
		                            p1.y = ri1.y;
		                            p1.z = ri1.z;
		                        }
		                        else if (t0 > dp1){
		                            p1.x = rj1.x;
		                            p1.y = rj1.y;
		                            p1.z = rj1.z;
		                        }
		                        if (t1 < 0){
		                            p2.x = ri2.x;
		                            p2.y = ri2.y;
		                            p2.z = ri2.z;
		                        }
		                        else if (t1 > dp2){
		                            p2.x = rj2.x;
		                            p2.y = rj2.y;
		                            p2.z = rj2.z;
		                        }
		                        float dot;
		                        if (t0 < 0 || t0 > dp1){
		                            dot = p2_norm.x*(p1.x - ri2.x) + p2_norm.y*(p1.y - ri2.y) + p2_norm.z*(p1.z - ri2.z);
		                            if (dot < 0){
		                                dot = 0.0;
		                            }
		                            else if (dot > dp2){
		                                dot = dp2;
		                            }
		                            p2.x = ri2.x + dot*p2_norm.x;
		                            p2.y = ri2.y + dot*p2_norm.y;
		                            p2.z = ri2.z + dot*p2_norm.z;
		                        }

		                        if (t1 < 0 || t1 > dp2){
		                            dot = p1_norm.x*(p2.x - ri1.x) + p1_norm.y*(p2.y - ri1.y) + p1_norm.z*(p2.z - ri1.z);
		                            if (dot < 0){
		                                dot = 0.0;
		                            }
		                            else if (dot > dp1){
		                                dot = dp1;
		                            }
		                            p1.x = ri1.x + dot*p1_norm.x;
		                            p1.y = ri1.y + dot*p1_norm.y;
		                            p1.z = ri1.z + dot*p1_norm.z;
		                        }

		                        intersect1 = p1;
		                        intersect2 = p2;
		                        dx = p2.x - p1.x;
		                        dy = p2.y - p1.y;
		                        dz = p2.z - p1.z;
		                        min_dist = sqrtf(dx*dx + dy*dy + dz*dz);
		                    }

		                    if ((min_dist < MT_RADIUS + KIN_RADIUS + 50.0) && (myDistance(intersect1, r[k*MT_NUM*MAX_MT_LENGTH + m*MAX_MT_LENGTH + MAX_MT_LENGTH - 1]) > 100.0*MT_RADIUS)){
		                    	int mt_dom = POLE_COUNT*MT_NUM*MAX_MT_LENGTH + KIN_COUNT*kn + 2*k*MT_NUM + m*2;
		                    	int ch_dom = POLE_COUNT*MT_NUM*MAX_MT_LENGTH + KIN_COUNT*kn + 2*k*MT_NUM + m*2 + 1;
		                    	h_att_l[mt_dom] = mi*(pol[k][m].length/(MAX_MT_LENGTH - 1)) + myDistance(intersect1, r[i]);
		                    	//printf("h_att_l: %f\n", h_att_l[mt_dom]);

		                    	center[ch_dom].x = intersect2.x - KIN_RADIUS*dx/min_dist;
		                    	center[ch_dom].y = intersect2.y - KIN_RADIUS*dy/min_dist;
		                    	center[ch_dom].z = intersect2.z - KIN_RADIUS*dz/min_dist;

		                    	center[mt_dom].x = intersect1.x + MT_RADIUS*dx/min_dist;
		                    	center[mt_dom].y = intersect1.y + MT_RADIUS*dy/min_dist;
		                    	center[mt_dom].z = intersect1.z + MT_RADIUS*dz/min_dist;
		                    	//printf("INTERSECT: %f\t%f\n", center[mt_dom].x, center[ch_dom].x);

		                    	connect_1[mt_dom] = i;
		                    	connect_2[mt_dom] = j;

		                    	connect_1[ch_dom] = js;
		                    	connect_2[ch_dom] = js + 1;
		                    	//printf("DISTANCE: %f\n%f\t%f\t%f\n%f\t%f\t%f\n", myDistance(r[j], center[mt_dom]), r[j].x, r[j].y, r[j].z, center[mt_dom].x, center[mt_dom].y, center[mt_dom].z);
		                    	len1[mt_dom] = myDistance(r[i], center[mt_dom]);
		                    	len2[mt_dom] = myDistance(r[j], center[mt_dom]);

		                    	len1[ch_dom] = myDistance(r[js], center[ch_dom]);
		                    	len2[ch_dom] = myDistance(r[js + 1], center[ch_dom]);

		                    	int3 kines_index = get_cell(COORtransfer1(center[mt_dom]));
                                if (kines_index.x > num_x_cell || kines_index.y > num_y_cell || kines_index.z > num_z_cell || kines_index.x <= 0 || kines_index.y <= 0 || kines_index.z <= 0){
				                    printf("mt_dyn_kin1: %d\n", i);
			                    }
								l = COORintoARRAY(kines_index.x, kines_index.y, kines_index.z);

		                    	kinesin_det_conc[l] = 1;
		                    }
						}
		                for (int ch = n_chrom/2 + 2; ch <= (n_chrom); ch++){
		                    int js = POLE_COUNT*MT_NUM*MAX_MT_LENGTH + (p + 1)*kn - ch;
		                    float3 ri2 = r[js]; //b0
				            float3 rj2 = r[js + 1];  //b1
		                    float3 p_vector2, p2_norm;
		                    p_vector2.x = rj2.x - ri2.x; //B
		                    p_vector2.y = rj2.y - ri2.y;
		                    p_vector2.z = rj2.z - ri2.z;

		                    float dp2 = sqrtf(p_vector2.x*p_vector2.x + p_vector2.y*p_vector2.y + p_vector2.z*p_vector2.z); //magB
		                    p2_norm.x = p_vector2.x/dp2; //_B
		                    p2_norm.y = p_vector2.y/dp2;
		                    p2_norm.z = p_vector2.z/dp2;

		                    float3 cross;
		                    cross.x = p1_norm.y*p2_norm.z - p2_norm.y*p1_norm.z; //cross
		                    cross.y = p1_norm.z*p2_norm.x - p2_norm.z*p1_norm.x;
		                    cross.z = p1_norm.x*p2_norm.y - p2_norm.x*p1_norm.y;

		                    float d_cross = sqrtf(cross.x*cross.x + cross.y*cross.y + cross.z*cross.z); //denom
		                    d_cross = d_cross*d_cross;
		                    if (d_cross == 0){
		                        float d0 = p1_norm.x*(ri2.x - ri1.x) + p1_norm.y*(ri2.y - ri1.y) + p1_norm.z*(ri2.z - ri1.z);
		                        float d1 = p1_norm.x*(rj2.x - ri1.x) + p1_norm.y*(rj2.y - ri1.y) + p1_norm.z*(rj2.z - ri1.z);
		                        if (d0 <= 0 && 0 >= d1){
		                            if (sqrtf(d0*d0) < sqrtf(d1*d1)){
		                                intersect1 = ri1;
		                                intersect2 = ri2;
		                                dx = intersect2.x - intersect1.x;
		                                dy = intersect2.y - intersect1.y;
		                                dz = intersect2.z - intersect1.z;
		                                min_dist = sqrtf(dx*dx + dy*dy + dz*dz);
		                            }
		                            else    {
		                                intersect1 = ri1;
		                                intersect2 = rj2;
		                                
		                                dx = intersect2.x - intersect1.x;
		                                dy = intersect2.y - intersect1.y;
		                                dz = intersect2.z - intersect1.z;
		                                min_dist = sqrtf(dx*dx + dy*dy + dz*dz);
		                            }
		                        }
		                        else if (d0 >= dp1 && dp1 <= d1){
		                            if (sqrtf(d0*d0) < sqrtf(d1*d1)){
		                                intersect1 = rj1;
		                                intersect2 = ri2;
		                                dx = intersect2.x - intersect1.x;
		                                dy = intersect2.y - intersect1.y;
		                                dz = intersect2.z - intersect1.z;
		                                min_dist = sqrtf(dx*dx + dy*dy + dz*dz);
		                            }
		                            else    {
		                                intersect1 = rj1;
		                                intersect2 = rj2;
		                                dx = intersect2.x - intersect1.x;
		                                dy = intersect2.y - intersect1.y;
		                                dz = intersect2.z - intersect1.z;
		                                min_dist = sqrtf(dx*dx + dy*dy + dz*dz);
		                            }
		                        }
		                    }
		                    else {
		                        float3 t;
		                        t.x = ri2.x - ri1.x;
		                        t.y = ri2.y - ri1.y;
		                        t.z = ri2.z - ri1.z;

		                        float det0 = t.x*(p2_norm.y*cross.z - cross.y*p2_norm.z) + t.y*(p2_norm.z*cross.x - cross.z*p2_norm.x) + t.z*(p2_norm.x*cross.y - cross.x*p2_norm.y);
		                        float det1 = t.x*(p1_norm.y*cross.z - cross.y*p1_norm.z) + t.y*(p1_norm.z*cross.x - cross.z*p1_norm.x) + t.z*(p1_norm.x*cross.y - cross.x*p1_norm.y);
		                        float t0 = det0/d_cross;
		                        float t1 = det1/d_cross;

		                        float3 p1, p2;
		                        p1.x = ri1.x + t0*p1_norm.x;
		                        p1.y = ri1.y + t0*p1_norm.y;
		                        p1.z = ri1.z + t0*p1_norm.z;

		                        p2.x = ri2.x + t1*p2_norm.x;
		                        p2.y = ri2.y + t1*p2_norm.y;
		                        p2.z = ri2.z + t1*p2_norm.z;

		                        if (t0 < 0){
		                            p1.x = ri1.x;
		                            p1.y = ri1.y;
		                            p1.z = ri1.z;
		                        }
		                        else if (t0 > dp1){
		                            p1.x = rj1.x;
		                            p1.y = rj1.y;
		                            p1.z = rj1.z;
		                        }
		                        if (t1 < 0){
		                            p2.x = ri2.x;
		                            p2.y = ri2.y;
		                            p2.z = ri2.z;
		                        }
		                        else if (t1 > dp2){
		                            p2.x = rj2.x;
		                            p2.y = rj2.y;
		                            p2.z = rj2.z;
		                        }
		                        float dot;
		                        if (t0 < 0 || t0 > dp1){
		                            dot = p2_norm.x*(p1.x - ri2.x) + p2_norm.y*(p1.y - ri2.y) + p2_norm.z*(p1.z - ri2.z);
		                            if (dot < 0){
		                                dot = 0.0;
		                            }
		                            else if (dot > dp2){
		                                dot = dp2;
		                            }
		                            p2.x = ri2.x + dot*p2_norm.x;
		                            p2.y = ri2.y + dot*p2_norm.y;
		                            p2.z = ri2.z + dot*p2_norm.z;
		                        }

		                        if (t1 < 0 || t1 > dp2){
		                            dot = p1_norm.x*(p2.x - ri1.x) + p1_norm.y*(p2.y - ri1.y) + p1_norm.z*(p2.z - ri1.z);
		                            if (dot < 0){
		                                dot = 0.0;
		                            }
		                            else if (dot > dp1){
		                                dot = dp1;
		                            }
		                            p1.x = ri1.x + dot*p1_norm.x;
		                            p1.y = ri1.y + dot*p1_norm.y;
		                            p1.z = ri1.z + dot*p1_norm.z;
		                        }

		                        intersect1 = p1;
		                        intersect2 = p2;
		                        dx = p2.x - p1.x;
		                        dy = p2.y - p1.y;
		                        dz = p2.z - p1.z;
		                        min_dist = sqrtf(dx*dx + dy*dy + dz*dz);
		                    }
		                    
		                    if ((min_dist < MT_RADIUS + KIN_RADIUS + 50.0) && (myDistance(intersect1, r[k*MT_NUM*MAX_MT_LENGTH + m*MAX_MT_LENGTH + MAX_MT_LENGTH - 1]) > 100.0*MT_RADIUS)){
		                    	int mt_dom = POLE_COUNT*MT_NUM*MAX_MT_LENGTH + KIN_COUNT*kn + 2*k*MT_NUM + m*2;
		                    	int ch_dom = POLE_COUNT*MT_NUM*MAX_MT_LENGTH + KIN_COUNT*kn + 2*k*MT_NUM + m*2 + 1;
		                    	h_att_l[mt_dom] = mi*(pol[k][m].length/(MAX_MT_LENGTH - 1)) + myDistance(intersect1, r[i]);
		                    	printf("h_att_l: %f\n", h_att_l[mt_dom]);

		                    	center[ch_dom].x = intersect2.x - KIN_RADIUS*dx/min_dist;
		                    	center[ch_dom].y = intersect2.y - KIN_RADIUS*dy/min_dist;
		                    	center[ch_dom].z = intersect2.z - KIN_RADIUS*dz/min_dist;

		                    	center[mt_dom].x = intersect1.x + MT_RADIUS*dx/min_dist;
		                    	center[mt_dom].y = intersect1.y + MT_RADIUS*dy/min_dist;
		                    	center[mt_dom].z = intersect1.z + MT_RADIUS*dz/min_dist;
		                    	//printf("INTERSECT: %f\t%f\n", center[mt_dom].x, center[ch_dom].x);

		                    	connect_1[mt_dom] = i;
		                    	connect_2[mt_dom] = j;

		                    	connect_1[ch_dom] = js;
		                    	connect_2[ch_dom] = js + 1;
		                    	//printf("DISTANCE: %f\n%f\t%f\t%f\n%f\t%f\t%f\n", myDistance(r[j], center[mt_dom]), r[j].x, r[j].y, r[j].z, center[mt_dom].x, center[mt_dom].y, center[mt_dom].z);
		                    	len1[mt_dom] = myDistance(r[i], center[mt_dom]);
		                    	len2[mt_dom] = myDistance(r[j], center[mt_dom]);

		                    	len1[ch_dom] = myDistance(r[js], center[ch_dom]);
		                    	len2[ch_dom] = myDistance(r[js + 1], center[ch_dom]);

		                    	int3 kines_index = get_cell(COORtransfer1(center[mt_dom]));
                                if (kines_index.x > num_x_cell || kines_index.y > num_y_cell || kines_index.z > num_z_cell || kines_index.x <= 0 || kines_index.y <= 0 || kines_index.z <= 0){
				                    printf("mt_dyn_kin2: %d\n", i);
			                    }
								l = COORintoARRAY(kines_index.x, kines_index.y, kines_index.z);

		                    	kinesin_det_conc[l] = 1;
		                    }
		                }
		            }
		        }
		    }
		    i = POLE_COUNT*MT_NUM*MAX_MT_LENGTH + KIN_COUNT*kn + 2*k*MT_NUM + m*2;
		    if (type[i] == LIG_MT_ACTIVE){
		    	int3 kines_index = get_cell(COORtransfer1(r[i]));
                if (kines_index.x > num_x_cell || kines_index.y > num_y_cell || kines_index.z > num_z_cell || kines_index.x <= 0 || kines_index.y <= 0 || kines_index.z <= 0){
                    printf("mt_dyn_kin3: %d\n", i);
                }
				l = COORintoARRAY(kines_index.x, kines_index.y, kines_index.z);
            	kinesin_att_conc[l] = 1;
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
            if (kinindex.x > num_x_cell || kinindex.y > num_y_cell || kinindex.z > num_z_cell || kinindex.x <= 0 || kinindex.y <= 0 || kinindex.z <= 0){
                printf("NDC: %d\n", i);
            }
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

	//int AAN = floor(ApV*b);
	int AAN = 0;
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
					int portionAA = floor(100*AAN*exp(-rijk*rijk/(2*400*400))/(400*sqrt(2*M_PI)));
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
	for (int m = 0; m < KIN_COUNT; m++){
		int n_chrom = 2*floor((CONTOUR_LENGTH/2 - KIN_RADIUS)/(2*KIN_RADIUS) + 0.5);
		float3 rcom = r0[POLE_COUNT*MT_NUM*MAX_MT_LENGTH + (m + 1)*kt[0].n - n_chrom - 1];

		int n_ab = floor(ABN/2);
		int number = m + 1;
		number_AB[0] = 0;
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

					float rijk = sqrt(pow((r.x - rcom.x), 2) + pow((r.y - rcom.y), 2) + pow((r.z - rcom.z), 2));
					//int portionAB = floor(ABN*exp(-rijk*rijk/(2*(2*KIN_RADIUS + sv_size)*(2*KIN_RADIUS + sv_size)))/((2*KIN_RADIUS + sv_size)*sqrt(2*M_PI)));
					int portionAB = floor(100*n_ab*exp(-rijk*rijk/(2*200*200))/(200*sqrt(2*M_PI)));
					if (cell_type[h] != OUTER_SPACE){
						aurora_B_conc[h] += portionAB;
						number_AB[0] += portionAB;
					}
				}
			}
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
	for (int m = 0; m < KIN_COUNT; m++){
		int n_chrom = 2*floor((CONTOUR_LENGTH/2 - KIN_RADIUS)/(2*KIN_RADIUS) + 0.5);
		float3 rcom = r0[POLE_COUNT*MT_NUM*MAX_MT_LENGTH + (m + 1)*kt[0].n - n_chrom - 1];

		int n_ab = number_AB;
		n_ab = n_ab/2;
		int number = m + 1;
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

					float rijk = sqrt(pow((r.x - rcom.x), 2) + pow((r.y - rcom.y), 2) + pow((r.z - rcom.z), 2));
					//int portionAB = floor(ABN*exp(-rijk*rijk/(2*(2*KIN_RADIUS + sv_size)*(2*KIN_RADIUS + sv_size)))/((2*KIN_RADIUS + sv_size)*sqrt(2*M_PI)));
					int portionAB = floor(100*n_ab*exp(-rijk*rijk/(2*200*200))/(200*sqrt(2*M_PI)));
					if (n_ab > 0 && portionAB < n_ab && cell_type[h] != OUTER_SPACE){
						if (rijk > KIN_RADIUS && portionAB > 0){
							float c = ran2(&rseed);
							int shift = int(10.0*c - 5);
							portionAB += shift;
							if (portionAB < 0){
								portionAB = 0;
							}
						}
						aurora_B_conc[h] += portionAB;
						//number_AB += portionAB;
						n_ab -= portionAB;
					}
					//n_ab -= portionAB;
				}
			}
		}
		if (n_ab > 0){
			int3 null = get_cell(COORtransfer1(rcom));
            if (null.x > num_x_cell || null.y > num_y_cell || null.z > num_z_cell || null.x <= 0 || null.y <= 0 || null.z <= 0){
                printf("null: %d\n", i);
            }
			n = COORintoARRAY(null.x, null.y, null.z);
			aurora_B_conc[n] += n_ab;
			n_ab -= n_ab;
		}
	}
}
