#include "kt_ex_vol.cuh"		

__global__ void excl_vol_kt(float3* r, float3* f, int N, int* type, int kn, Param* d_parameters){
	int i, j;
    int N_pol = d_parameters->pole_n;
    int N_mt = d_parameters->mt_n;
    int N_kt = d_parameters->kin_n;
    int N_mt_max = d_parameters->max_mt_n;
    float mt_rad = d_parameters->mt_r;
    float kt_rad = d_parameters->kt_r;
    float e_lj = d_parameters->e_rep/0.23;
    int chrom_cond = d_parameters->chrom;
    int chrom_bn = chrom_cond*d_parameters->chrom_num;
    int x_num = d_parameters->x_num;
    int y_num = d_parameters->y_num;
    //float e_mor = d_parameters->e_att/0.23;
	i = blockIdx.x*blockDim.x + threadIdx.x;
	if (i < N){
		if (type[i] == LEFT_KT){
			float3 ri = r[i];
            for ((j = N_pol*N_mt*N_mt_max + kn - 1); j < N; j += kn){
                if (i != j && j != (i + kn)){
                    float3 rj = r[j];
                    float dx = rj.x - ri.x;
	                float dy = rj.y - ri.y;
	                float dz = rj.z - ri.z;
                    float dr = sqrtf(dx*dx + dy*dy + dz*dz);
                    if (dr < 2*(kt_rad + mt_rad)){
                        float df = -6*powf((2*(mt_rad + kt_rad)/dr), 12)*powf(dr, -1)*powf(dr, -1)*e_lj;
                        atomicAdd(&f[i].x, df*dx);
	                    atomicAdd(&f[i].y, df*dy);  		
	                    atomicAdd(&f[i].z, df*dz);

	                    atomicAdd(&f[j].x, -df*dx);
	                    atomicAdd(&f[j].y, -df*dy);		
	                    atomicAdd(&f[j].z, -df*dz);
                    }
                }
            }
            int k = i + kn;
            float3 rk = r[k];
            float dx = rk.x - ri.x;
            float dy = rk.y - ri.y;
            float dz = rk.z - ri.z;
            float dr = sqrtf(dx*dx + dy*dy + dz*dz);
            if (dr < 2*(kt_rad)){
                float df = -6*powf((2*(kt_rad)/dr), 12)*powf(dr, -1)*powf(dr, -1)*e_lj;
                atomicAdd(&f[i].x, df*dx);
                atomicAdd(&f[i].y, df*dy);  		
                atomicAdd(&f[i].z, df*dz);

                atomicAdd(&f[k].x, -df*dx);
                atomicAdd(&f[k].y, -df*dy);		
                atomicAdd(&f[k].z, -df*dz);
            }
		}
        if (type[i] == RIGHT_KT){
			float3 ri = r[i];
            for ((j = N_pol*N_mt*N_mt_max + kn - 1); j < N; j += kn){
                if (i != j && j != (i - kn)){
                    float3 rj = r[j];
                    float dx = rj.x - ri.x;
	                float dy = rj.y - ri.y;
	                float dz = rj.z - ri.z;
                    float dr = sqrtf(dx*dx + dy*dy + dz*dz);
                    if (dr < 2*(kt_rad + mt_rad)){
                        float df = -6*powf((2*(mt_rad + kt_rad)/dr), 12)*powf(dr, -1)*powf(dr, -1)*e_lj;
                        atomicAdd(&f[i].x, df*dx);
	                    atomicAdd(&f[i].y, df*dy);  		
	                    atomicAdd(&f[i].z, df*dz);

	                    atomicAdd(&f[j].x, -df*dx);
	                    atomicAdd(&f[j].y, -df*dy);		
	                    atomicAdd(&f[j].z, -df*dz);
                    }
                }
            }
		}
    	if (type[i] == CHROM){
			float3 ri = r[i];
            int ks;
            int c_id = i - N_pol*N_mt*N_mt_max;
            int k_id = int(c_id/kn);
            for (ks = 0; ks < N_kt; ks++){ 
                for (int ch = 1; ch <= chrom_bn + 1; ch++){
                    int k = N_pol*N_mt*N_mt_max + (ks + 1)*kn - ch;
                    if (ks != k_id){
                        float3 rj = r[k];
                        float dx = rj.x - ri.x;
	                    float dy = rj.y - ri.y;
	                    float dz = rj.z - ri.z;
                        float dr = sqrtf(dx*dx + dy*dy + dz*dz);
                        if (dr < 2*kt_rad){
                            float df = -6*powf((2*kt_rad/dr), 12)*powf(dr, -1)*powf(dr, -1)*e_lj;
                            atomicAdd(&f[i].x, df*dx);
	                        atomicAdd(&f[i].y, df*dy);  		
	                        atomicAdd(&f[i].z, df*dz);

	                        atomicAdd(&f[k].x, -df*dx);
	                        atomicAdd(&f[k].y, -df*dy);		
	                        atomicAdd(&f[k].z, -df*dz);
                        }
                    }
                }
            }
            int ic = i - (kn - x_num*y_num) - 0.5*x_num*y_num;
             
            float3 rc = r[ic];
            float dx = rc.x - ri.x;
            float dy = rc.y - ri.y;
            float dz = rc.z - ri.z;
            float dr = sqrtf(dx*dx + dy*dy + dz*dz);
            if (dr < sqrtf(2.0)*kt_rad){
                float df = -6*powf((sqrtf(2.0)*kt_rad/dr), 12)*powf(dr, -1)*powf(dr, -1)*e_lj;
                atomicAdd(&f[i].x, df*dx);
                atomicAdd(&f[i].y, df*dy);  		
                atomicAdd(&f[i].z, df*dz);

                /*atomicAdd(&f[ic].x, -df*dx);
                atomicAdd(&f[ic].y, -df*dy);		
                atomicAdd(&f[ic].z, -df*dz);*/
            }
        }
    	//__syncthreads();
	}
}	

__global__ void chrom_ter(float3* r, float3* f, int N, int* type, int kn, Param* d_parameters){
	int i;
    float pol_dist = d_parameters->pole_dist;
    float pole_rad = d_parameters->pole_r;
    float kt_rad = d_parameters->kt_r;
    //float e_mor = d_parameters->e_att/0.23;
	i = blockIdx.x*blockDim.x + threadIdx.x;
	if (i < N){
		if (type[i] == LEFT_KT || type[i] == RIGHT_KT || type[i] == CHROM){
			float re2 = pow(r[i].x, 2)*pow(pol_dist - 1.0625*pole_rad - 1.5*kt_rad, -2) + pow(r[i].y, 2)*pow(pol_dist - 1.0625*pole_rad - 1.5*kt_rad, -2) + pow(r[i].z, 2)*pow(pol_dist - 1.0625*pole_rad - 1.5*kt_rad, -2); 
			if (re2 > 1.0){
                float dr = sqrtf(r[i].x*r[i].x + r[i].y*r[i].y + r[i].z*r[i].z);
				atomicAdd(&f[i].x, -100000*(re2 - 1.0)*r[i].x*powf(dr, -1));
                atomicAdd(&f[i].y, -100000*(re2 - 1.0)*r[i].y*powf(dr, -1));  		
                atomicAdd(&f[i].z, -100000*(re2 - 1.0)*r[i].z*powf(dr, -1));
			}
		}
	}
}	
__global__ void chrom_memb(float3* r, float3* f, int N, int* type, int kn, Param* d_parameters){
	int i;
    float kt_rad = d_parameters->kt_r;
    float mt_rad = d_parameters->mt_r;
    float a_rad = d_parameters->a;
    float b_rad = d_parameters->b;
    float c_rad = d_parameters->c;
    //float e_mor = d_parameters->e_att/0.23;
	i = blockIdx.x*blockDim.x + threadIdx.x;
	if (i < N){
		if (type[i] == LEFT_KT || type[i] == RIGHT_KT || type[i] == CHROM){
			float re2 = pow(r[i].x, 2)*pow(a_rad - kt_rad, -2) + pow(r[i].y, 2)*pow(b_rad - kt_rad, -2) + pow(r[i].z, 2)*pow(c_rad - kt_rad, -2); 
			if (re2 > 1.0){
                float dr = sqrtf(r[i].x*r[i].x + r[i].y*r[i].y + r[i].z*r[i].z);
				atomicAdd(&f[i].x, -100000*(re2 - 1.0)*r[i].x*powf(dr, -1));
                atomicAdd(&f[i].y, -100000*(re2 - 1.0)*r[i].y*powf(dr, -1));  		
                atomicAdd(&f[i].z, -100000*(re2 - 1.0)*r[i].z*powf(dr, -1));
			}
		}
        if (type[i] == MT_REG || type[i] == PLUS_DET || type[i] == PLUS_ATT || type[i] == PLUS_DET_INVALID){
			float re2 = pow(r[i].x, 2)*pow(a_rad - 2*mt_rad, -2) + pow(r[i].y, 2)*pow(b_rad - 2*mt_rad, -2) + pow(r[i].z, 2)*pow(c_rad - 2*mt_rad, -2); 
			if (re2 > 1.0){
                float dr = sqrtf(r[i].x*r[i].x + r[i].y*r[i].y + r[i].z*r[i].z);
				atomicAdd(&f[i].x, -100000*(re2 - 1.0)*r[i].x*powf(dr, -1));
                atomicAdd(&f[i].y, -100000*(re2 - 1.0)*r[i].y*powf(dr, -1));  		
                atomicAdd(&f[i].z, -100000*(re2 - 1.0)*r[i].z*powf(dr, -1));
			}
		}
	}
}		
