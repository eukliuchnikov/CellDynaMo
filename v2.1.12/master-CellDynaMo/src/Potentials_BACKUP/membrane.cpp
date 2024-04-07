#include "membrane.h"		

void membrane(float3* f, float* m_length){
	int m, j; 
	for (j = 0; j < POLE_COUNT; j++){
		for (m = 0; m < MT_NUM; m++){
			if (pol[j][m].state != MT_STATE_DEPOL){
				int MT_nums = pol[j][m].cyl_N;
				float3 r = pol[j][m].r[MT_nums];
				int i = j*MT_NUM*MAX_MT_LENGTH + m*MAX_MT_LENGTH + MT_nums;
				if (MEM_SHAPE == MEM_SPHERE || MEM_SHAPE == MEM_ELIPSOID){
					float re2 = pow(r.x, 2)*pow(ELLIPSE_A, -2) + pow(r.y, 2)*pow(ELLIPSE_B, -2) + pow(r.z, 2)*pow(ELLIPSE_C, -2); 
					if (re2 > 1.0){
						float3 normal;
						normal.x = - 2*r.x*pow(ELLIPSE_A, -2);
						normal.y = - 2*r.y*pow(ELLIPSE_B, -2);
						normal.z = - 2*r.z*pow(ELLIPSE_C, -2);
						float cosa = scalar(normal, r)*pow(length(normal), -1)*pow(length(r), -1);
						//if (cosa < cos(M_PI/6)){
							float dx = pol[j][m].r[MT_nums].x - pol[j][m].r[MT_nums - 1].x;
							float dy = pol[j][m].r[MT_nums].y - pol[j][m].r[MT_nums - 1].y;
							float dz = pol[j][m].r[MT_nums].z - pol[j][m].r[MT_nums - 1].z;
							float dr = sqrt(dx*dx + dy*dy + dz*dz);
							dx = 2.0*MT_RADIUS*dx/dr;
							dy = 2.0*MT_RADIUS*dy/dr;
							dz = 2.0*MT_RADIUS*dz/dr;
							pol[j][m].state = MT_STATE_DEPOL;
							pol[j][m].r[MT_nums].x -= dx;
							pol[j][m].r[MT_nums].y -= dy;
							pol[j][m].r[MT_nums].z -= dz;
							pol[j][m].r[MT_nums - 1].x -= dx;
							pol[j][m].r[MT_nums - 1].y -= dy;
							pol[j][m].r[MT_nums - 1].z -= dz;
							pol[j][m].length -= 2.0*MT_RADIUS;
							m_length[i] = pol[j][m].length;
							m_length[i - 1] = pol[j][m].length - (MT_nums - 1)*1000.0;
	                        f[mds.N - 2 + j].x -= 10*(re2 - 1.0)*r.x/length(r);
	                        f[mds.N - 2 + j].y -= 10*(re2 - 1.0)*r.y/length(r);
	                        f[mds.N - 2 + j].z -= 10*(re2 - 1.0)*r.z/length(r);
						/*}
						else{
							f[i].x -= K_STEFF_M*(re2 - 1.0)*r.x*pow(length(r), -1);
							f[i].y -= K_STEFF_M*(re2 - 1.0)*r.y*pow(length(r), -1);		
							f[i].z -= K_STEFF_M*(re2 - 1.0)*r.z*pow(length(r), -1); 
						}*/		
					}
				}
				else	if (MEM_SHAPE == MEM_CUBE || MEM_SHAPE == MEM_CUBOID){
					if (r.x > ELLIPSE_A/2 || r.x < -ELLIPSE_A/2 || r.y > ELLIPSE_B/2 || r.y < -ELLIPSE_B/2 || r.z > ELLIPSE_C/2 || r.z < -ELLIPSE_C/2){
						float dx = pol[j][m].r[MT_nums].x - pol[j][m].r[MT_nums - 1].x;
						float dy = pol[j][m].r[MT_nums].y - pol[j][m].r[MT_nums - 1].y;
						float dz = pol[j][m].r[MT_nums].z - pol[j][m].r[MT_nums - 1].z;
						float dr = sqrt(dx*dx + dy*dy + dz*dz);
						dx = 2.0*MT_RADIUS*dx/dr;
						dy = 2.0*MT_RADIUS*dy/dr;
						dz = 2.0*MT_RADIUS*dz/dr;
						pol[j][m].state = MT_STATE_DEPOL;
						pol[j][m].r[MT_nums].x -= dx;
						pol[j][m].r[MT_nums].y -= dy;
						pol[j][m].r[MT_nums].z -= dz;
						pol[j][m].r[MT_nums - 1].x -= dx;
						pol[j][m].r[MT_nums - 1].y -= dy;
						pol[j][m].r[MT_nums - 1].z -= dz;
						pol[j][m].length -= 2.0*MT_RADIUS;
						m_length[i] = pol[j][m].length;
						m_length[i - 1] = pol[j][m].length - (MT_nums - 1)*1000.0;
                        /*f[mds.N - 2 + j].x -= 10*(re2 - 1.0)*r.x/length(r);
                        f[mds.N - 2 + j].y -= 10*(re2 - 1.0)*r.y/length(r);
                        f[mds.N - 2 + j].z -= 10*(re2 - 1.0)*r.z/length(r);*/
                    }
				}
			}
		}
	}
}