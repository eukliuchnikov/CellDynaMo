/*
 *pushing module
 *
 *  Created on: Mar 15, 2018
 *	Edited on:	Mar 30, 2019
 *      Author: kliuchnikov
 */

#include "pushing.cuh"		

__global__ void pushing(float3* r, float3* f, int N, int* type, int kn, long long int step, int* n_force, float* f_force, Param* d_parameters){
	int i, j, h;
    int N_pol = d_parameters->pole_n;
    int N_mt = d_parameters->mt_n;
    int N_mt_max = d_parameters->max_mt_n;
    int N_kt = d_parameters->kin_n;
    float mt_rad = d_parameters->mt_r;
    float kt_rad = d_parameters->kt_r;
    float e_lj = d_parameters->e_rep/0.23;
    float zone = d_parameters->zone;
    float big_r = d_parameters->big_r;
    int chrom_cond = d_parameters->chrom;
    int x_num = d_parameters->x_num;
    int y_num = d_parameters->y_num;
    int chrom_bn = chrom_cond*d_parameters->chrom_num;
	i = blockIdx.x*blockDim.x + threadIdx.x;
	if (i < N){
		if (type[i] == PLUS_DET || type[i] == MT_REG || type[i] == LEFT_MINUS || type[i] == RIGHT_MINUS || type[i] == PLUS_ATT || type[i] == PLUS_DET_INVALID){
			float3 ri = r[i];
		    for (h = 0; h < N_kt; h++){
		        j = N_pol*N_mt*N_mt_max + (h + 1)*kn - 1 - chrom_bn;

		        float3 rj = r[j];

		        float dx = rj.x - ri.x;
		        float dy = rj.y - ri.y;
		        float dz = rj.z - ri.z;

		        float dr = sqrtf(dx*dx + dy*dy + dz*dz);

                float high = big_r - kt_rad + zone;
			    float limit = sqrtf(pow(zone, 2) + pow(big_r, 2) - pow(high, 2));
                if (dr < limit + 4*mt_rad){
                    float3 kt1;
				    float3 kt2;
                    if (h % 2 == 0){
					    kt1 = r[N_pol*N_mt*N_mt_max + (h + 1)*kn - 1 - chrom_bn];
					    kt2 = r[N_pol*N_mt*N_mt_max + (h + 2)*kn - 1 - chrom_bn];
				    }
				    else	{
					    kt1 = r[N_pol*N_mt*N_mt_max + h*kn - 1 - chrom_bn];
					    kt2 = r[N_pol*N_mt*N_mt_max + (h + 1)*kn - 1 - chrom_bn];
				    }
                    float3 normal;
			        normal.x = kt2.x - kt1.x;
			        normal.y = kt2.y - kt1.y;
			        normal.z = kt2.z - kt1.z; 
			        float dn = sqrt(normal.x*normal.x + normal.y*normal.y + normal.z*normal.z);
                    float3 delta;
			        delta.x = zone*normal.x/dn;
			        delta.y = zone*normal.y/dn;
			        delta.z = zone*normal.z/dn;

			        float3 right, left;
			        right.x = kt2.x + delta.x;
			        right.y = kt2.y + delta.y;
			        right.z = kt2.z + delta.z;

			        left.x = kt1.x - delta.x;
			        left.y = kt1.y - delta.y;
			        left.z = kt1.z - delta.z;
                    
                    if ((normal.x*(ri.x - right.x) + normal.y*(ri.y - right.y) + normal.z*(ri.z - right.z)) < 0 && (normal.x*(ri.x - left.x) + normal.y*(ri.y - left.y) + normal.z*(ri.z - left.z)) > 0){
                        if (dr > limit - mt_rad){
                            float df = -6*e_lj*powf(limit/dr, 12)*powf(dr, -1)*powf(dr, -1);
                            atomicAdd(&f[i].x, df*dx);
	                        atomicAdd(&f[i].y, df*dy);  		
	                        atomicAdd(&f[i].z, df*dz);

	                        atomicAdd(&f[j].x, -df*dx);
	                        atomicAdd(&f[j].y, -df*dy);		
	                        atomicAdd(&f[j].z, -df*dz);

                             if (type[j] == LEFT_KT  && step == 0){
                                atomicAdd(&n_force[h*5], 1);
                                atomicAdd(&f_force[h*5 + 1], -df*1.66);
                                atomicAdd(&n_force[h*5 + 2], 1); 
                            }
                            else    if (type[j] == RIGHT_KT  && step == 0){
                                atomicAdd(&n_force[h*5], 1);
                                atomicAdd(&f_force[h*5 + 1], -df*1.66);
                                atomicAdd(&n_force[h*5 + 2], 1); 
                            }
                        }       
                    }
                    else    {
                        float minR = 2*mt_rad;
                        float maxR = kt_rad;
                        for (int ks = N_pol*N_mt*N_mt_max + h*kn; ks < N_pol*N_mt*N_mt_max + (h + 1)*kn - 1 - chrom_bn; ks++){
                            float3 rk = r[ks];
                            float dxk = rk.x - ri.x;
		                    float dyk = rk.y - ri.y;
		                    float dzk = rk.z - ri.z;

		                    float drk = sqrtf(dxk*dxk + dyk*dyk + dzk*dzk);
                            if (drk < minR){
                                float dxc = rk.x - rj.x;
		                        float dyc = rk.y - rj.y;
		                        float dzc = rk.z - rj.z;

		                        float drc = sqrtf(dxc*dxc + dyc*dyc + dzc*dzc);
                                minR = drk;
                                maxR = drc;
                            }
                        } 
                        if (minR < 2*mt_rad && dr > maxR - mt_rad){
                    		float df = -6*e_lj*powf(maxR/dr, 12)*powf(dr, -1)*powf(dr, -1);                            
                            atomicAdd(&f[i].x, df*dx);
	                        atomicAdd(&f[i].y, df*dy);  		
	                        atomicAdd(&f[i].z, df*dz);

	                        atomicAdd(&f[j].x, -df*dx);
	                        atomicAdd(&f[j].y, -df*dy);		
	                        atomicAdd(&f[j].z, -df*dz);

                            if (type[j] == LEFT_KT  && step == 0){
                                atomicAdd(&n_force[h*5], 1);
                                atomicAdd(&f_force[h*5 + 1], -df*1.66);
                                atomicAdd(&n_force[h*5 + 2], 1); 
                            }
                            else    if (type[j] == RIGHT_KT  && step == 0){
                                atomicAdd(&n_force[h*5], 1);
                                atomicAdd(&f_force[h*5 + 1], -df*1.66);
                                atomicAdd(&n_force[h*5 + 2], 1); 
                            }
                        }                  
                    }
                }
		    }
            for (int p = 0; p < N_kt; p++){
                for (int ch = 1; ch <= chrom_bn; ch++){
                    int k = N_pol*N_mt*N_mt_max + (p + 1)*kn - ch;

                    float3 rj = r[k];

		            float dx = rj.x - ri.x;
		            float dy = rj.y - ri.y;
		            float dz = rj.z - ri.z;

		            float dr = sqrtf(dx*dx + dy*dy + dz*dz);
                    if (dr < kt_rad + 2*mt_rad && dr > kt_rad - mt_rad){
                        float df = -6*e_lj*powf((kt_rad + 3*mt_rad)/dr, 12)*powf(dr, -1)*powf(dr, -1)/10;
                        //printf("%f\t%f\t%f\t%f\n", e_lj, dr, kt_rad + 3*mt_rad/dr, df);
                        atomicAdd(&f[i].x, df*dx);
                        atomicAdd(&f[i].y, df*dy);  		
                        atomicAdd(&f[i].z, df*dz);

                        atomicAdd(&f[k].x, -df*dx);
                        atomicAdd(&f[k].y, -df*dy);		
                        atomicAdd(&f[k].z, -df*dz);
                        if (type[k - 5 + ch] == LEFT_KT/*  && step == 0*/){
                            atomicAdd(&n_force[p*5], 1);
                            atomicAdd(&f_force[p*5 + 1], -df*1.66);
                            atomicAdd(&n_force[p*5 + 2], 1);   
                                 
                            /*atomicAdd(&n_force[p*5 + 1], 1);
                            atomicAdd(&f_force[p*5 + 1], -df*1.66);
                            atomicAdd(&n_force[p*5 + 2], 1);*/
                        }
                        else    if (type[k - 5 + ch] == RIGHT_KT/*  && step == 0*/){
                            atomicAdd(&n_force[p*5], 1);
                            atomicAdd(&f_force[p*5 + 1], -df*1.66);
                            atomicAdd(&n_force[p*5 + 2], 1);

                            /*atomicAdd(&n_force[p*5], 1);
                            atomicAdd(&f_force[p*5], -df*1.66);
                            atomicAdd(&n_force[p*5 + 3], 1);*/
                        }
                    }
                }
                for (int sh = x_num*y_num; sh < kn; sh++){
                    int k = N_pol*N_mt*N_mt_max + p*kn + sh;
    
                    float3 rj = r[k];

		            float dx = rj.x - ri.x;
		            float dy = rj.y - ri.y;
		            float dz = rj.z - ri.z;

		            float dr = sqrtf(dx*dx + dy*dy + dz*dz);
                    if (dr < 80.0){
                        float3 kt1;
				        float3 kt2;
                        if (p % 2 == 0){
					        kt1 = r[N_pol*N_mt*N_mt_max + (p + 1)*kn - 1 - chrom_bn];
					        kt2 = r[N_pol*N_mt*N_mt_max + (p + 2)*kn - 1 - chrom_bn];
				        
	                        float dx1 = kt1.x - ri.x;
		                    float dy1 = kt1.y - ri.y;
		                    float dz1 = kt1.z - ri.z;

	                        float dx2 = kt2.x - ri.x;
		                    float dy2 = kt2.y - ri.y;
		                    float dz2 = kt2.z - ri.z;
                            
		                    float dr1 = sqrtf(dx1*dx1 + dy1*dy1 + dz1*dz1);
		                    float dr2 = sqrtf(dx2*dx2 + dy2*dy2 + dz2*dz2);

                            if (dr1 < dr2){
                                float df = -6*e_lj*powf(80.0/dr, 12)*powf(dr, -1)*powf(dr, -1)/10;
                                //printf("%f\n, df");
                                if (df > 100.0){
                                    df = 100.0;
                                } 
                                if (df < -100.0){
                                    df = -100.0;
                                }                                                               
                                //printf("%f\t%f\t%f\t%f\n", e_lj, dr, kt_rad + 3*mt_rad/dr, df);
                                atomicAdd(&f[i].x, df*dx);
                                atomicAdd(&f[i].y, df*dy);  		
                                atomicAdd(&f[i].z, df*dz);

                                atomicAdd(&f[N_pol*N_mt*N_mt_max + (p + 1)*kn - 1 - chrom_bn].x, -df*dx);
                                atomicAdd(&f[N_pol*N_mt*N_mt_max + (p + 1)*kn - 1 - chrom_bn].y, -df*dy);		
                                atomicAdd(&f[N_pol*N_mt*N_mt_max + (p + 1)*kn - 1 - chrom_bn].z, -df*dz);
                                
                                atomicAdd(&n_force[p*5], 1);
                                atomicAdd(&f_force[p*5 + 1], -df*1.66);
                                atomicAdd(&n_force[p*5 + 2], 1);   
                            }
                            if (dr2 < dr1){
                                float df = -6*e_lj*powf(80.0/dr, 12)*powf(dr, -1)*powf(dr, -1)/10;
                                if (df > 100.0){
                                    df = 100.0;
                                }
                                if (df < -100.0){
                                    df = -100.0;
                                } 
                                //printf("%f\t%f\t%f\t%f\n", e_lj, dr, kt_rad + 3*mt_rad/dr, df);
                                atomicAdd(&f[i].x, df*dx);
                                atomicAdd(&f[i].y, df*dy);  		
                                atomicAdd(&f[i].z, df*dz);

                                atomicAdd(&f[N_pol*N_mt*N_mt_max + (p + 2)*kn - 1 - chrom_bn].x, -df*dx);
                                atomicAdd(&f[N_pol*N_mt*N_mt_max + (p + 2)*kn - 1 - chrom_bn].y, -df*dy);		
                                atomicAdd(&f[N_pol*N_mt*N_mt_max + (p + 2)*kn - 1 - chrom_bn].z, -df*dz);
                              
                                atomicAdd(&n_force[p*5], 1);
                                atomicAdd(&f_force[p*5 + 1], -df*1.66);
                                atomicAdd(&n_force[p*5 + 2], 1);
                            }
                        }
                    }
                }
            }
        }
    	//__syncthreads();
	}
}		
