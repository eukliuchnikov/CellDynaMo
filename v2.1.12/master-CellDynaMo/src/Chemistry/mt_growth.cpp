#include "mt_growth.h"

void GROWTHfunc(int l, float3* r, int* type, int & minimiz, float* length, int* connect_1, int* connect_2, float* len1, float* len2, float* h_att_l){
	int m, j, i, q, mi;
	int p;
	int countMT = 0;
	int rseed_kin = rseed + 123;
	for (j = 0; j < POLE_COUNT; j++){
		for (m = 0; m < MT_NUM; m++){
			i = j*MT_NUM*MAX_MT_LENGTH + m*MAX_MT_LENGTH + MAX_MT_LENGTH - 1;
			float3 ri = r[i];
			int3 mtindex = get_cell(COORtransfer1(ri));
			if (mtindex.x > num_x_cell || mtindex.y > num_y_cell || mtindex.z > num_z_cell || mtindex.x <= 0 || mtindex.y <= 0 || mtindex.z <= 0){
				printf("growth1: %d\n", i);
			}
			q = COORintoARRAY(mtindex.x, mtindex.y, mtindex.z);
			if (l == q && pol[j][m].state == MT_STATE_POL && type[i] != PLUS_ATT){
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
				printf("growth2: %d\n", i);
			}
			q = COORintoARRAY(mtindex.x, mtindex.y, mtindex.z);
			if (l == q && pol[j][m].state == MT_STATE_POL && type[i] != PLUS_ATT){
				if (pol[j][m].rand - 1 < c && c <= pol[j][m].rand){
					//printf("mt index: %d\n", i);
					for (mi = 1; mi < MAX_MT_LENGTH; mi ++){
						if (mi == 1){
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

							pol[j][m].r[MAX_MT_LENGTH - mi].x += dx0;
							pol[j][m].r[MAX_MT_LENGTH - mi].y += dy0;
							pol[j][m].r[MAX_MT_LENGTH - mi].z += dz0;
						}
						else	{
							float dx = pol[j][m].r[MAX_MT_LENGTH - mi + 1].x - pol[j][m].r[MAX_MT_LENGTH - mi].x;
							float dy = pol[j][m].r[MAX_MT_LENGTH - mi + 1].y - pol[j][m].r[MAX_MT_LENGTH - mi].y;
							float dz = pol[j][m].r[MAX_MT_LENGTH - mi + 1].z - pol[j][m].r[MAX_MT_LENGTH - mi].z;
							float dr = sqrt(dx*dx + dy*dy + dz*dz);
							dx = 2.0*MT_RADIUS*dx/dr;
							dy = 2.0*MT_RADIUS*dy/dr;
							dz = 2.0*MT_RADIUS*dz/dr;

							float dx0 = (MAX_MT_LENGTH - mi)*dx/(MAX_MT_LENGTH - 1);
							float dy0 = (MAX_MT_LENGTH - mi)*dy/(MAX_MT_LENGTH - 1);
							float dz0 = (MAX_MT_LENGTH - mi)*dz/(MAX_MT_LENGTH - 1);
							float dr0 = sqrt(dx0*dx0 + dy0*dy0 + dz0*dz0);

							pol[j][m].r[MAX_MT_LENGTH - mi].x += dx0;
							pol[j][m].r[MAX_MT_LENGTH - mi].y += dy0;
							pol[j][m].r[MAX_MT_LENGTH - mi].z += dz0;
						}
					}
					pol[j][m].length += 2.0*MT_RADIUS;
					length[i] = pol[j][m].length;
					for (mi = 1; mi < MAX_MT_LENGTH; mi ++){
						length[i - mi] = pol[j][m].length;
					}
					int mt_dom = POLE_COUNT*MT_NUM*MAX_MT_LENGTH + KIN_COUNT*kn + 2*j*MT_NUM + m*2;
					
					float3 r_new = pol[j][m].r[MAX_MT_LENGTH - 1];
					if (MEM_SHAPE == MEM_SPHERE || MEM_SHAPE == MEM_ELIPSOID){
						float re2 = pow(r_new.x, 2)*pow(ELLIPSE_A, -2) + pow(r_new.y, 2)*pow(ELLIPSE_B, -2) + pow(r_new.z, 2)*pow(ELLIPSE_C, -2); 
						if (re2 > 1.0){
							pol[j][m].state = MT_STATE_DEPOL;
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
							for (mi = 1; mi < MAX_MT_LENGTH; mi++){
								length[i - mi] = pol[j][m].length;
							}
							int mt_dom = POLE_COUNT*MT_NUM*MAX_MT_LENGTH + KIN_COUNT*kn + 2*j*MT_NUM + m*2;
						}
					}
					else	if (MEM_SHAPE == MEM_CUBE || MEM_SHAPE == MEM_CUBOID){
						if (r_new.x > ELLIPSE_A/2 || r_new.x < -ELLIPSE_A/2 || r_new.y > ELLIPSE_B/2 || r_new.y < -ELLIPSE_B/2 || r_new.z > ELLIPSE_C/2 || r_new.z < -ELLIPSE_C/2){
							pol[j][m].state = MT_STATE_DEPOL;
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
						}
					}
					for (p = 0; p < KIN_COUNT; p++){
						int n_chrom = 2*floor((CONTOUR_LENGTH/2 - KIN_RADIUS)/(2*KIN_RADIUS) + 0.5);
    					int k = POLE_COUNT*MT_NUM*MAX_MT_LENGTH + (p + 1)*kn - n_chrom - 1;
    					/*for (int chr = 1; chr < n_chrom; chr++){
    						int c_id = POLE_COUNT*MT_NUM*MAX_MT_LENGTH + (p + 1)*kn - chr;
    						float3 rj = r[c_id];
							float dxK = rj.x - pol[j][m].r[MAX_MT_LENGTH - 1].x;
							float dyK = rj.y - pol[j][m].r[MAX_MT_LENGTH - 1].y;
							float dzK = rj.z - pol[j][m].r[MAX_MT_LENGTH - 1].z;
							float drK = sqrtf(dxK*dxK + dyK*dyK + dzK*dzK);

							if (drK < KIN_RADIUS + 2*MT_RADIUS){
							//pol[j][m].state = MT_STATE_DEPOL;
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
									length[i + mi] = pol[j][m].length;
								}
								minimiz = 1;	
							}
    					}*/
    					float3 rj = r[k];
						float dxK = rj.x - pol[j][m].r[MAX_MT_LENGTH - 1].x;
						float dyK = rj.y - pol[j][m].r[MAX_MT_LENGTH - 1].y;
						float dzK = rj.z - pol[j][m].r[MAX_MT_LENGTH - 1].z;
						float drK = sqrtf(dxK*dxK + dyK*dyK + dzK*dzK);

						float3 kt1;
						float3 kt2;
						if (p % 2 == 0){
							kt1 = r[POLE_COUNT*MT_NUM*MAX_MT_LENGTH + (p + 1)*kn - n_chrom - 1];
							kt2 = r[POLE_COUNT*MT_NUM*MAX_MT_LENGTH + (p + 2)*kn - n_chrom - 1];
						}
						else	{
							kt1 = r[POLE_COUNT*MT_NUM*MAX_MT_LENGTH + p*kn - n_chrom - 1];
							kt2 = r[POLE_COUNT*MT_NUM*MAX_MT_LENGTH + (p + 1)*kn - n_chrom - 1];
						}
						float3 normal;
						normal.x = kt2.x - kt1.x;
						normal.y = kt2.y - kt1.y;
						normal.z = kt2.z - kt1.z; 
						float dn = sqrt(normal.x*normal.x + normal.y*normal.y + normal.z*normal.z);
						float3 delta;
						delta.x = NO_ATT_ZONE*normal.x/dn;
						delta.y = NO_ATT_ZONE*normal.y/dn;
						delta.z = NO_ATT_ZONE*normal.z/dn;
						float3 right, left;
						right.x = kt2.x + delta.x;
						right.y = kt2.y + delta.y;
						right.z = kt2.z + delta.z;

						left.x = kt1.x - delta.x;
						left.y = kt1.y - delta.y;
						left.z = kt1.z - delta.z;
						float big_r;
						if (BETA >= 0.925){
							big_r = 5000.0;
						}
						else	if (BETA >= 0 && BETA < 0.925){
							big_r = KIN_RADIUS/(1.0 - BETA);
						}
						float high = big_r - KIN_RADIUS + NO_ATT_ZONE;
						float limit = sqrtf(pow(NO_ATT_ZONE, 2) + pow(big_r, 2) - pow(high, 2));
						if (drK >= KIN_RADIUS && drK < limit + 4*MT_RADIUS){
							minimiz = 1;
						}
						if (drK < KIN_RADIUS || drK < limit && ((normal.x*(pol[j][m].r[MAX_MT_LENGTH - 1].x - right.x) + normal.y*(pol[j][m].r[MAX_MT_LENGTH - 1].y - right.y) + normal.z*(pol[j][m].r[MAX_MT_LENGTH - 1].z - right.z) < 0) && (normal.x*(pol[j][m].r[MAX_MT_LENGTH - 1].x - left.x) + normal.y*(pol[j][m].r[MAX_MT_LENGTH - 1].y - left.y) + normal.z*(pol[j][m].r[MAX_MT_LENGTH - 1].z - left.z) > 0))){
							pol[j][m].state = MT_STATE_DEPOL;
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
							minimiz = 1;	
						}
						if (p % 2 == 0){
							int k_next = k + kn;
							float3 rjn = r[k_next];
							float3 r_center;
							r_center.x = (rj.x + rjn.x)/2;
							r_center.y = (rj.y + rjn.y)/2;
							r_center.z = (rj.z + rjn.z)/2;

							if (myDistance(r_center, pol[j][m].r[MAX_MT_LENGTH - 1]) < 2*KIN_RADIUS){
								pol[j][m].state = MT_STATE_DEPOL;
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
								minimiz = 1;	
							}
						}
						/*int x_num = CORONA_X;
						int y_num = floor(x_num/RATIO);
						for (int cr = x_num*y_num; cr < kn; cr++){
							int c_id = POLE_COUNT*MT_NUM*MAX_MT_LENGTH + p*kn + cr;
							float3 rc = r[c_id];
							float dxc = rc.x - pol[j][m].r[MAX_MT_LENGTH - 1].x;
							float dyc = rc.y - pol[j][m].r[MAX_MT_LENGTH - 1].y;
							float dzc = rc.z - pol[j][m].r[MAX_MT_LENGTH - 1].z;
							float drc = sqrtf(dxc*dxc + dyc*dyc + dzc*dzc);
							if (drc < 60.0){
								pol[j][m].state = MT_STATE_DEPOL;
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
									length[i + mi] = pol[j][m].length;
								}
								minimiz = 1;	
							}
						}*/
					}
					if (connect_1[mt_dom] != 0 && type[mt_dom] == LIG_MT_ACTIVE){
	                	//float min_dist = 10000;


	                	//int min_ID = connect_1[mt_dom];

	                	int min_ID = floor(h_att_l[mt_dom]/(length[i]/(MAX_MT_LENGTH - 1)));


	                	/*for (mi = 1; mi < MAX_MT_LENGTH; mi ++){
	                		if (myDistance(pol[j][m].r[MAX_MT_LENGTH - mi], r[mt_dom]) < min_dist){
	                			min_dist = myDistance(pol[j][m].r[MAX_MT_LENGTH - mi], r[mt_dom]);
	                			min_ID = j*MT_NUM*MAX_MT_LENGTH + m*MAX_MT_LENGTH + MAX_MT_LENGTH - mi;
	                		}
	                	}*/
	                	connect_1[mt_dom] = j*MT_NUM*MAX_MT_LENGTH + m*MAX_MT_LENGTH + min_ID + 1;
	                	connect_2[mt_dom] = j*MT_NUM*MAX_MT_LENGTH + m*MAX_MT_LENGTH + min_ID;

	                	len2[mt_dom] = h_att_l[mt_dom] - min_ID*(length[i]/(MAX_MT_LENGTH - 1));
	                	len1[mt_dom] = length[i]/(MAX_MT_LENGTH - 1) - len2[mt_dom];
	                	/*len1[mt_dom] = myDistance(r[connect_1[mt_dom]], r[mt_dom]);
	                	if (connect_1[mt_dom] != i && myDistance(r[connect_1[mt_dom] - 1], r[mt_dom]) >= myDistance(r[connect_1[mt_dom] + 1], r[mt_dom])){
	                		connect_2[mt_dom] = min_ID + 1;
	                	}
	                	else	if (myDistance(r[connect_1[mt_dom] - 1], r[mt_dom]) < myDistance(r[connect_1[mt_dom] + 1], r[mt_dom])){
	                		connect_2[mt_dom] = min_ID - 1;
	                	}
	                	len2[mt_dom] = length[i]/(MAX_MT_LENGTH - 1) - len1[mt_dom];*/
        			}
				}
			}
		}
	}
	minimiz = 1;
}
