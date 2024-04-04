#include "mt_growth.h"

void GROWTHfunc(int l, float3* r, int* type, int & minimiz, float* length){
	int m, j, i, q;
	int p;
	int countMT = 0;
	int rseed_kin = rseed + 123;
	for (j = 0; j < POLE_COUNT; j++){
		for (m = 0; m < MT_NUM; m++){
			i = j*MT_NUM*MAX_MT_LENGTH + m*MAX_MT_LENGTH + MAX_MT_LENGTH - 1;
			float3 ri = r[i];
			int3 mtindex = get_cell(COORtransfer1(ri));
			if (mtindex.x > num_x_cell || mtindex.y > num_y_cell || mtindex.z > num_z_cell || mtindex.x <= 0 || mtindex.y <= 0 && mtindex.z <= 0){
				printf("growth: %d\n", i);
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
			if (mtindex.x > num_x_cell || mtindex.y > num_y_cell || mtindex.z > num_z_cell || mtindex.x <= 0 || mtindex.y <= 0 && mtindex.z <= 0){
				printf("growth: %d\n", i);
			}
			q = COORintoARRAY(mtindex.x, mtindex.y, mtindex.z);
			if (l == q && pol[j][m].state == MT_STATE_POL && type[i] != PLUS_ATT){
				if (pol[j][m].rand - 1 < c && c <= pol[j][m].rand){
					//printf("mt index: %d\n", i);
					float dx = pol[j][m].r[2].x - pol[j][m].r[1].x;
					float dy = pol[j][m].r[2].y - pol[j][m].r[1].y;
					float dz = pol[j][m].r[2].z - pol[j][m].r[1].z;
					float dr = sqrt(dx*dx + dy*dy + dz*dz);
					dx = 2.0*MT_RADIUS*dx/dr;
					dy = 2.0*MT_RADIUS*dy/dr;
					dz = 2.0*MT_RADIUS*dz/dr;
					float dx0 = pol[j][m].r[1].x - pol[j][m].r[0].x;
					float dy0 = pol[j][m].r[1].y - pol[j][m].r[0].y;
					float dz0 = pol[j][m].r[1].z - pol[j][m].r[0].z;
					float dr0 = sqrt(dx0*dx0 + dy0*dy0 + dz0*dz0);
					dx0 = MT_RADIUS*dx0/dr0;
					dy0 = MT_RADIUS*dy0/dr0;
					dz0 = MT_RADIUS*dz0/dr0;
					pol[j][m].r[2].x += dx;
					pol[j][m].r[2].y += dy;
					pol[j][m].r[2].z += dz;
					pol[j][m].r[1].x += dx0;
					pol[j][m].r[1].y += dy0;
					pol[j][m].r[1].z += dz0;
					pol[j][m].length += 2.0*MT_RADIUS;
					length[i] = pol[j][m].length;
					length[i - 1] = pol[j][m].length;
					length[i - 2] = pol[j][m].length;
					float3 r_new = pol[j][m].r[2];
					if (MEM_SHAPE == MEM_SPHERE || MEM_SHAPE == MEM_ELIPSOID){
						float re2 = pow(r_new.x, 2)*pow(ELLIPSE_A, -2) + pow(r_new.y, 2)*pow(ELLIPSE_B, -2) + pow(r_new.z, 2)*pow(ELLIPSE_C, -2); 
						if (re2 > 1.0){
							pol[j][m].state = MT_STATE_DEPOL;
							pol[j][m].r[2].x -= dx;
							pol[j][m].r[2].y -= dy;
							pol[j][m].r[2].z -= dz;
							pol[j][m].r[1].x -= dx0;
							pol[j][m].r[1].y -= dy0;
							pol[j][m].r[1].z -= dz0;
							pol[j][m].length -= 2.0*MT_RADIUS;
							length[i] = pol[j][m].length;
							length[i - 1] = pol[j][m].length;
							length[i - 2] = pol[j][m].length;
						}
					}
					else	if (MEM_SHAPE == MEM_CUBE || MEM_SHAPE == MEM_CUBOID){
						if (r_new.x > ELLIPSE_A/2 || r_new.x < -ELLIPSE_A/2 || r_new.y > ELLIPSE_B/2 || r_new.y < -ELLIPSE_B/2 || r_new.z > ELLIPSE_C/2 || r_new.z < -ELLIPSE_C/2){
							pol[j][m].state = MT_STATE_DEPOL;
							pol[j][m].r[2].x -= dx;
							pol[j][m].r[2].y -= dy;
							pol[j][m].r[2].z -= dz;
							pol[j][m].r[1].x -= dx0;
							pol[j][m].r[1].y -= dy0;
							pol[j][m].r[1].z -= dz0;
							pol[j][m].length -= 2.0*MT_RADIUS;
							length[i] = pol[j][m].length;
							length[i - 1] = pol[j][m].length;
							length[i - 2] = pol[j][m].length;
						}
					}
					for (p = 0; p < KIN_COUNT; p++){
						int n_chrom = 2*floor((CONTOUR_LENGTH/2 - KIN_RADIUS)/(2*KIN_RADIUS) + 0.5);
    					int k = POLE_COUNT*MT_NUM*MAX_MT_LENGTH + (p + 1)*kn - n_chrom - 1;
    					float3 rj = r[k];
						float dxK = rj.x - pol[j][m].r[2].x;
						float dyK = rj.y - pol[j][m].r[2].y;
						float dzK = rj.z - pol[j][m].r[2].z;
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
						if (drK < KIN_RADIUS || drK < limit && ((normal.x*(pol[j][m].r[2].x - right.x) + normal.y*(pol[j][m].r[2].y - right.y) + normal.z*(pol[j][m].r[2].z - right.z) < 0) && (normal.x*(pol[j][m].r[2].x - left.x) + normal.y*(pol[j][m].r[2].y - left.y) + normal.z*(pol[j][m].r[2].z - left.z) > 0))){
							pol[j][m].state = MT_STATE_DEPOL;
							pol[j][m].r[2].x -= dx;
							pol[j][m].r[2].y -= dy;
							pol[j][m].r[2].z -= dz;
							pol[j][m].r[1].x -= dx0;
							pol[j][m].r[1].y -= dy0;
							pol[j][m].r[1].z -= dz0;
							pol[j][m].length -= 2.0*MT_RADIUS;
							length[i] = pol[j][m].length;
							length[i - 1] = pol[j][m].length;
							length[i - 2] = pol[j][m].length;
							minimiz = 1;	
						}
						if (p % 2 == 0){
							int k_next = k + kn;
							float3 rjn = r[k_next];
							float3 r_center;
							r_center.x = (rj.x + rjn.x)/2;
							r_center.y = (rj.y + rjn.y)/2;
							r_center.z = (rj.z + rjn.z)/2;

							if (myDistance(r_center, pol[j][m].r[2]) < 2*KIN_RADIUS){
								pol[j][m].state = MT_STATE_DEPOL;
								pol[j][m].r[2].x -= dx;
								pol[j][m].r[2].y -= dy;
								pol[j][m].r[2].z -= dz;
								pol[j][m].r[1].x -= dx0;
								pol[j][m].r[1].y -= dy0;
								pol[j][m].r[1].z -= dz0;
								pol[j][m].length -= 2.0*MT_RADIUS;
								length[i] = pol[j][m].length;
								length[i - 1] = pol[j][m].length;
								length[i - 2] = pol[j][m].length;
								minimiz = 1;	
							}
						}
						int x_num = CORONA_X;
						int y_num = floor(x_num/RATIO);
						for (int cr = x_num*y_num; cr < kn; cr++){
							int c_id = POLE_COUNT*MT_NUM*MAX_MT_LENGTH + p*kn + cr;
							float3 rc = r[c_id];
							float dxc = rc.x - pol[j][m].r[2].x;
							float dyc = rc.y - pol[j][m].r[2].y;
							float dzc = rc.z - pol[j][m].r[2].z;
							float drc = sqrtf(dxc*dxc + dyc*dyc + dzc*dzc);
							if (drc < 60.0){
								pol[j][m].state = MT_STATE_DEPOL;
								pol[j][m].r[2].x -= dx;
								pol[j][m].r[2].y -= dy;
								pol[j][m].r[2].z -= dz;
								pol[j][m].r[1].x -= dx0;
								pol[j][m].r[1].y -= dy0;
								pol[j][m].r[1].z -= dz0;
								pol[j][m].length -= 2.0*MT_RADIUS;
								length[i] = pol[j][m].length;
								length[i - 1] = pol[j][m].length;
								length[i - 2] = pol[j][m].length;
								minimiz = 1;	
							}
						}
					}
				}
			}
		}
	}
}