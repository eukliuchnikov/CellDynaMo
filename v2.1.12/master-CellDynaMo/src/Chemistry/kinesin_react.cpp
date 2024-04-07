#include "kinesin_react.h"

void kinesin_attach(float3* r, int* type, int l, float3* center, int* connect_1){
	int h, m;
	int n_chrom = 2*floor((CONTOUR_LENGTH/2 - KIN_RADIUS)/(2*KIN_RADIUS) + 0.5);
	for (h = 0; h < POLE_COUNT; h++){	
		for (m = 0; m < MT_NUM; m++){
			int i = POLE_COUNT*MT_NUM*MAX_MT_LENGTH + KIN_COUNT*kn + 2*h*MT_NUM + m*2;
			if (type[i] == LIG_MT_INACTIVE && connect_1[i] != 0){
				type[i] = LIG_MT_ACTIVE;
				type[i + 1] = LIG_CH_ACTIVE;
				pol[h][m].kinesin = 1;
				r[i] = center[i];
				r[i + 1] = center[i + 1];
				//printf("ATTACHED: %f\t%f\n", r[i].x, r[i + 1].x);
			}
		}
	}
}

void kinesin_detach(float3* r, int* type, int l, float3* center, int* connect_1, int* connect_2, float* len1, float* len2, float* att_l){
	int h, m;
	int countMT = 0;
	int rseed_kin = rseed + 123;
	for (h = 0; h < POLE_COUNT; h++){	
		for (m = 0; m < MT_NUM; m++){
			int i = POLE_COUNT*MT_NUM*MAX_MT_LENGTH + KIN_COUNT*kn + 2*h*MT_NUM + m*2;
			float3 ri = r[i];
			int3 mtindex = get_cell(COORtransfer1(ri));
			if (mtindex.x > num_x_cell || mtindex.y > num_y_cell || mtindex.z > num_z_cell || mtindex.x <= 0 || mtindex.y <= 0 || mtindex.z <= 0){
				printf("kin_react1: %d\n", i);
			}
			int s = COORintoARRAY(mtindex.x, mtindex.y, mtindex.z);
			if (s == l && type[i] == LIG_MT_ACTIVE && pol[h][m].kinesin > 0){
				countMT++;
				pol[h][m].rand = countMT;
			}
		}
	}
	float c = countMT*ran2(&rseed_kin);
	for (h = 0; h < POLE_COUNT; h++){	
		for (m = 0; m < MT_NUM; m++){
			int i = POLE_COUNT*MT_NUM*MAX_MT_LENGTH + KIN_COUNT*kn + 2*h*MT_NUM + m*2;
			float3 ri = r[i];
			int3 mtindex = get_cell(COORtransfer1(ri));
			if (type[i] == LIG_MT_INACTIVE && pol[h][m].kinesin > 0){
				pol[h][m].kinesin = 0;
			}
            if (mtindex.x > num_x_cell || mtindex.y > num_y_cell || mtindex.z > num_z_cell || mtindex.x <= 0 || mtindex.y <= 0 || mtindex.z <= 0){
				printf("kin_react2: %d\n", i);
			}
			int s = COORintoARRAY(mtindex.x, mtindex.y, mtindex.z);
			if (s == l && type[i] == LIG_MT_ACTIVE && pol[h][m].kinesin > 0){
				if (pol[h][m].rand - 1 < c && c <= pol[h][m].rand){
					pol[h][m].kinesin = 0;
					type[i] = LIG_MT_INACTIVE;
					type[i + 1] = LIG_CH_INACTIVE;
					connect_1[i] = 0;
					connect_2[i] = 0;
					len1[i] = 0.0;
					len2[i] = 0.0;
					center[i].x = POLE_DISTANCE;
					center[i].y = 0.0;
					center[i].z = 0.0;
					r[i] = center[i];
					att_l[i] = 0.0;

					connect_1[i + 1] = 0;
					connect_2[i + 1] = 0;
					len1[i + 1] = 0.0;
					len2[i + 1] = 0.0;
					center[i + 1].x = POLE_DISTANCE;
					center[i + 1].y = 0.0;
					center[i + 1].z = 0.0;
					r[i + 1] = center[i + 1];
					att_l[i + 1] = 0.0;
				}
			}
		}
	}
}	

void kinesin_step(float3* r, int* type, float* length, int l, float3* center, int* connect_1, int* connect_2, float* len1, float* len2, float* f_out, float* att_l){
	int h, m, mi;
	int countMT = 0;
	int rseed_kin = rseed + 123;
	for (h = 0; h < POLE_COUNT; h++){	
		for (m = 0; m < MT_NUM; m++){
			int i = POLE_COUNT*MT_NUM*MAX_MT_LENGTH + KIN_COUNT*kn + 2*h*MT_NUM + m*2;
			float3 ri = r[i];
			int3 mtindex = get_cell(COORtransfer1(ri));
			if (mtindex.x > num_x_cell || mtindex.y > num_y_cell || mtindex.z > num_z_cell || mtindex.x <= 0 || mtindex.y <= 0 || mtindex.z <= 0){
				printf("kin_step1: %d\n", i);
			}	
			int s = COORintoARRAY(mtindex.x, mtindex.y, mtindex.z);
			if (s == l && type[i] == LIG_MT_ACTIVE){
				countMT++;
				pol[h][m].rand = countMT;
			}
		}
	}
	float c = countMT*ran2(&rseed_kin);
	for (h = 0; h < POLE_COUNT; h++){	
		for (m = 0; m < MT_NUM; m++){
			int i = POLE_COUNT*MT_NUM*MAX_MT_LENGTH + KIN_COUNT*kn + 2*h*MT_NUM + m*2;
			float3 ri = r[i];
			int3 mtindex = get_cell(COORtransfer1(ri));
			if (mtindex.x > num_x_cell || mtindex.y > num_y_cell || mtindex.z > num_z_cell || mtindex.x <= 0 || mtindex.y <= 0 || mtindex.z <= 0){
				printf("kin_step1: %d\n", i);
			}
			int s = COORintoARRAY(mtindex.x, mtindex.y, mtindex.z);
			if (s == l && type[i] == LIG_MT_ACTIVE){
				if (pol[h][m].rand - 1 < c && c <= pol[h][m].rand){
					int m1 = connect_1[i];
					float3 r_m1 = r[m1];
					int m2 = connect_2[i];            
	                float3 r_m12 = r[m2];
	                float dr2 = MT_RADIUS;
	                
	                float dx_diff = r_m12.x - r_m1.x;
	                float dy_diff = r_m12.y - r_m1.y;
	                float dz_diff = r_m12.z - r_m1.z;

	                float diff = sqrtf(dx_diff*dx_diff + dy_diff*dy_diff + dz_diff*dz_diff);
	        		
	        		float f0 = 15.0;
	        		float step = 2.0*MT_RADIUS*(1 - f_out[i]/f0);
	        		if (step < 0){
	        			step = 0.0;
	        		}
	        		//printf("%f\t%f\n", step, f_out[i]);
	                float3 r_m2;
	                if (m1 > m2){
		                r[i].x -= step*dx_diff/diff;
		                r[i].y -= step*dy_diff/diff;
		                r[i].z -= step*dz_diff/diff;
		            }
		            else	{
		            	r[i].x += step*dx_diff/diff;
		                r[i].y += step*dy_diff/diff;
		                r[i].z += step*dz_diff/diff;
		            }
		            att_l[i] += step;
	                if (connect_1[i] != 0 && type[i] == LIG_MT_ACTIVE){
	                	int min_ID = floor(att_l[i]/(length[m1]/(MAX_MT_LENGTH - 1)));

	                	connect_1[i] = h*MT_NUM*MAX_MT_LENGTH + m*MAX_MT_LENGTH + min_ID + 1;
	                	connect_2[i] = h*MT_NUM*MAX_MT_LENGTH + m*MAX_MT_LENGTH + min_ID;

	                	len2[i] = att_l[i] - min_ID*(length[m1]/(MAX_MT_LENGTH - 1));
	                	len1[i] = length[m1]/(MAX_MT_LENGTH - 1) - len2[i];
        			}

		            if (myDistance(pol[h][m].r[MAX_MT_LENGTH - 1], r[i]) < MT_RADIUS || f_out[i] > 20.0){
			            pol[h][m].kinesin = 0;
						type[i] = LIG_MT_INACTIVE;
						type[i + 1] = LIG_CH_INACTIVE;
						connect_1[i] = 0;
						connect_2[i] = 0;
						len1[i] = 0.0;
						len2[i] = 0.0;
						center[i].x = POLE_DISTANCE;
						center[i].y = 0.0;
						center[i].z = 0.0;
						r[i] = center[i];
						att_l[i] = 0.0;

						connect_1[i + 1] = 0;
						connect_2[i + 1] = 0;
						len1[i + 1] = 0.0;
						len2[i + 1] = 0.0;
						center[i + 1].x = POLE_DISTANCE;
						center[i + 1].y = 0.0;
						center[i + 1].z = 0.0;
						r[i + 1] = center[i + 1];
						att_l[i + 1] = 0.0;
		            }
				}
			}
		}
	}
}
