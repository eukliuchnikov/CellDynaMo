#include "pulling.cuh"		

__global__ void pulling(float3* r, float3* f, int N, int* type, int kn, long long int step, int* n_force, float* f_force, Param* d_parameters, int* connector, float* length, float* att_l, float* link_l){
	int i, j;
    int N_pol = d_parameters->pole_n;
    int N_mt = d_parameters->mt_n;
    int N_kt = d_parameters->kin_n;
    int N_mt_max = d_parameters->max_mt_n;
    float mt_rad = d_parameters->mt_r;
    float e_mor = d_parameters->e_att;
    int chrom_cond = d_parameters->chrom;
    int chrom_bn = chrom_cond*d_parameters->chrom_num;
    i = blockIdx.x*blockDim.x + threadIdx.x;
	if (i < N){
		if (type[i] == PLUS_ATT){
			float3 ri = r[i];
           
            float3 rp = r[i - 1];
            float dxp = rp.x - ri.x;
            float dyp = rp.y - ri.y;
            float dzp = rp.z - ri.z;

            float drp = sqrtf(dxp*dxp + dyp*dyp + dzp*dzp);
            float3 rpp = r[i - 2];
            float dxpp = rpp.x - rp.x;
            float dypp = rpp.y - rp.y;
            float dzpp = rpp.z - rp.z;

            float drpp = sqrtf(dxpp*dxpp + dypp*dypp + dzpp*dzpp);
            int pol;
            if (i < N_mt*N_mt_max){
                pol = N - 2;
            }
            else    if (i < N_pol*N_mt*N_mt_max){
                pol = N - 1;
            }
            if (drp + drpp > att_l[i] + mt_rad){
                length[i] = att_l[i];
                length[i - 1] = att_l[i];
                length[i - 2] = att_l[i];
                float dx0 = r[i - 2].x - r[pol].x;
			    float dy0 = r[i - 2].y - r[pol].y;
			    float dz0 = r[i - 2].z - r[pol].z;
			    float dr0 = sqrt(dx0*dx0 + dy0*dy0 + dz0*dz0);
                r[i].x = r[i - 1].x - (length[i]/2)*dxp/drp;
                r[i].y = r[i - 1].y - (length[i]/2)*dyp/drp;
                r[i].z = r[i - 1].z - (length[i]/2)*dyp/drp;

                ri = r[i];
            }
            if (connector[i] != 0){
                j = connector[i];
	            float3 rj = r[j];
	            float dx = rj.x - ri.x;
	            float dy = rj.y - ri.y;
	            float dz = rj.z - ri.z;

	            float dr = sqrtf(dx*dx + dy*dy + dz*dz);
              
                float3 r_pole1 = r[N - 2];
                float3 r_pole2 = r[N - 1];

                float dx1 = r_pole1.x - ri.x;
                float dy1 = r_pole1.y - ri.y;
                float dz1 = r_pole1.z - ri.z;

                float dx2 = r_pole2.x - ri.x;
                float dy2 = r_pole2.y - ri.y;
                float dz2 = r_pole2.z - ri.z;
    
                float dr1 = sqrtf(dx1*dx1 + dy1*dy1 + dz1*dz1);
                float dr2 = sqrtf(dx2*dx2 + dy2*dy2 + dz2*dz2);
                if (dr1 > 0.3*MIN_MT_LENGTH && dr2 > 0.3*MIN_MT_LENGTH){
                    float df = e_mor*(dr - link_l[i])*powf(dr, -1);
                    atomicAdd(&f[i].x, df*dx);
                    atomicAdd(&f[i].y, df*dy);		
                    atomicAdd(&f[i].z, df*dz);                   
                        
                    atomicAdd(&f[j].x, -df*dx);
	                atomicAdd(&f[j].y, -df*dy);		
	                atomicAdd(&f[j].z, -df*dz);
                    for (int ks = 0; ks < N_kt; ks++){
                        if (type[j] == LEFT_KT && j == N_pol*N_mt*N_mt_max + (ks + 1)*kn - 1 - chrom_bn  && step == 0){
                            atomicAdd(&n_force[ks*5], 1);
                            atomicAdd(&f_force[ks*5], df*1.66);
                            atomicAdd(&n_force[ks*5 + 1], 1);

                            /*atomicAdd(&n_force[ks*5], 1);
                            atomicAdd(&f_force[ks*5], df*1.66);
                            atomicAdd(&n_force[ks*5 + 4], 1);*/
                        }
                        else    if (type[j] == RIGHT_KT && j == N_pol*N_mt*N_mt_max + (ks + 1)*kn - 1 - chrom_bn  && step == 0){
                            atomicAdd(&n_force[ks*5], 1);
                            atomicAdd(&f_force[ks*5], df*1.66);
                            atomicAdd(&n_force[ks*5 + 1], 1);
                            
                            /*atomicAdd(&n_force[ks*5 + 1], 1);
                            atomicAdd(&f_force[ks*5 + 1], df*1.66);
                            atomicAdd(&n_force[ks*5 + 5], 1);*/
                        }
                    }
                }
            }
	    }
    	//__syncthreads();
	}
}		
