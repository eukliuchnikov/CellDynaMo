#include "harmonic.cuh"		

__global__ void computeHarmonic(float3* r, float3* f, int N, int* type, int kn, float rk, int* harmonicKinCount, int* harmonicKin, float* harmonicKinRadii, float* kt_radius, Param* d_parameters, float* length){
    int N_pol = d_parameters->pole_n;
    int N_mt = d_parameters->mt_n;
    int N_kt = d_parameters->kin_n;
    int N_mt_max = d_parameters->max_mt_n;
    float K_kt = d_parameters->k_kt;
    float K_mt = d_parameters->k_mt;
    int MpK = d_parameters->mHkPm;
    int chrom_cond = d_parameters->chrom;
    int chrom_bn = chrom_cond*d_parameters->chrom_num;
    int x_num = d_parameters->x_num;
    int y_num = d_parameters->y_num;
    float kt_rad = d_parameters->kt_r;
	int i;	
	i = blockIdx.x*blockDim.x + threadIdx.x;	
	if (i < N){
		if (type[i] == MT_REG || type[i] == PLUS_DET || type[i] == PLUS_ATT || type[i] == PLUS_DET_INVALID){
            float r0 = length[i]/2.0;
   			float3 ri = r[i - 1];
			float3 rj = r[i];

			float dx = rj.x - ri.x;
			float dy = rj.y - ri.y;
			float dz = rj.z - ri.z;
	        
			float dr = sqrtf(dx*dx + dy*dy + dz*dz);
            float df = K_mt*(dr - r0)*powf(dr, -1);

		    atomicAdd(&f[i - 1].x, df*dx);
		    atomicAdd(&f[i - 1].y, df*dy);		
		    atomicAdd(&f[i - 1].z, df*dz);
            
			atomicAdd(&f[i].x, -df*dx);
			atomicAdd(&f[i].y, -df*dy);		
			atomicAdd(&f[i].z, -df*dz);
		}
        /*if (type[i] == LEFT_MINUS){
            float r0 = length[i]/2.0;
   			float3 ri = r[N - 2];
			float3 rj = r[i];

			float dx = rj.x - ri.x;
			float dy = rj.y - ri.y;
			float dz = rj.z - ri.z;
	    
			float dr = sqrtf(dx*dx + dy*dy + dz*dz);
            float df = K_mt*(dr - r0)*powf(dr, -1);

		    atomicAdd(&f[N - 2].x, df*dx);
		    atomicAdd(&f[N - 2].y, df*dy);		
		    atomicAdd(&f[N - 2].z, df*dz);
            
			atomicAdd(&f[i].x, -df*dx);
			atomicAdd(&f[i].y, -df*dy);		
			atomicAdd(&f[i].z, -df*dz);
		}
        if (type[i] == RIGHT_MINUS){
            float r0 = length[i]/2.0;
   			float3 ri = r[N - 1];
			float3 rj = r[i];

			float dx = rj.x - ri.x;
			float dy = rj.y - ri.y;
			float dz = rj.z - ri.z;
	    
			float dr = sqrtf(dx*dx + dy*dy + dz*dz);
            float df = K_mt*(dr - r0)*powf(dr, -1);

		    atomicAdd(&f[N - 1].x, df*dx);
		    atomicAdd(&f[N - 1].y, df*dy);		
		    atomicAdd(&f[N - 1].z, df*dz);
            
			atomicAdd(&f[i].x, -df*dx);
			atomicAdd(&f[i].y, -df*dy);		
			atomicAdd(&f[i].z, -df*dz);
		}*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//**************************************************************************************************************************************************************//
//between centres of kinetochores
//**************************************************************************************************************************************************************//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		if (type[i] == LEFT_KT){
			float3 ri = r[i];
			int j = i + kn;
			float3 rj = r[j];

			float dx = rj.x - ri.x;
			float dy = rj.y - ri.y;
			float dz = rj.z - ri.z;

			float dr = sqrtf(dx*dx + dy*dy + dz*dz);
			float df = 0.5*K_kt*(dr - rk)*powf(dr, -1);

			atomicAdd(&f[i].x, df*dx);
			atomicAdd(&f[i].y, df*dy);		
			atomicAdd(&f[i].z, df*dz);

			atomicAdd(&f[j].x, -df*dx);
			atomicAdd(&f[j].y, -df*dy);		
			atomicAdd(&f[j].z, -df*dz);
            
            if (chrom_bn != 0){
                for (int chrom = 1; chrom < chrom_bn + 1; chrom++){
                    float3 rch1 = r[i + chrom];
			        int c1 = i + kn + chrom;
			        float3 rch2 = r[c1];

			        float dx1 = rch2.x - rch1.x;
			        float dy1 = rch2.y - rch1.y;
			        float dz1 = rch2.z - rch1.z;

			        float dr1 = sqrtf(dx1*dx1 + dy1*dy1 + dz1*dz1);
			        float df1 = 0.001*K_kt*(dr1 - rk)*powf(dr1, -1);
                    

                    atomicAdd(&f[i + chrom].x, df1*dx1);
			        atomicAdd(&f[i + chrom].y, df1*dy1);		
			        atomicAdd(&f[i + chrom].z, df1*dz1);

			        atomicAdd(&f[c1].x, -df1*dx1);
			        atomicAdd(&f[c1].y, -df1*dy1);		
			        atomicAdd(&f[c1].z, -df1*dz1);
                }

                for (int chrom = 1; chrom < chrom_bn; chrom += chrom_bn/2){
                    int k = i + chrom;
			        float3 rk = r[k];

			        float dxk = rk.x - ri.x;
			        float dyk = rk.y - ri.y;
			        float dzk = rk.z - ri.z;

			        float drk = sqrtf(dxk*dxk + dyk*dyk + dzk*dzk);
			        float dfk = 200*K_kt*(drk - kt_radius[k - N_pol*N_mt*N_mt_max])*powf(drk, -1);


			        atomicAdd(&f[i].x, dfk*dxk);
			        atomicAdd(&f[i].y, dfk*dyk);		
			        atomicAdd(&f[i].z, dfk*dzk);

			        atomicAdd(&f[k].x, -dfk*dxk);
			        atomicAdd(&f[k].y, -dfk*dyk);		
			        atomicAdd(&f[k].z, -dfk*dzk);
                }
            }
		}
        if (type[i] == RIGHT_KT && chrom_bn != 0){
            float3 ri = r[i];
            for (int chrom = 1; chrom < chrom_bn; chrom += chrom_bn/2){
                int k = i + chrom;
			    float3 rk = r[k];

			    float dxk = rk.x - ri.x;
			    float dyk = rk.y - ri.y;
			    float dzk = rk.z - ri.z;

			    float drk = sqrtf(dxk*dxk + dyk*dyk + dzk*dzk);
			    float dfk = 200*K_kt*(drk - kt_radius[k - N_pol*N_mt*N_mt_max])*powf(drk, -1);

			    atomicAdd(&f[i].x, dfk*dxk);
			    atomicAdd(&f[i].y, dfk*dyk);		
			    atomicAdd(&f[i].z, dfk*dzk);

			    atomicAdd(&f[k].x, -dfk*dxk);
			    atomicAdd(&f[k].y, -dfk*dyk);		
			    atomicAdd(&f[k].z, -dfk*dzk);
            }
		}
        if (type[i] == CHROM){

            float3 ri = r[i];
            int ks;
            for (ks = 0; ks < N_kt; ks ++){            
                if (i < N_pol*N_mt*N_mt_max + (ks + 1)*kn - chrom_bn/2 && i > N_pol*N_mt*N_mt_max + (ks + 1)*kn - chrom_bn || i < N_pol*N_mt*N_mt_max + (ks + 1)*kn && i > N_pol*N_mt*N_mt_max + (ks + 1)*kn - chrom_bn/2){
                    int k = i - 1;
			        float3 rk = r[k];

			        float dxk = rk.x - ri.x;
			        float dyk = rk.y - ri.y;
			        float dzk = rk.z - ri.z;

			        float drk = sqrtf(dxk*dxk + dyk*dyk + dzk*dzk);
			        float dfk = 200*K_kt*(drk - kt_radius[i - N_pol*N_mt*N_mt_max])*powf(drk, -1);
                    
                    //printf("%d\t%d\t%f\t%f\n", i, k, drk, kt_radius[i - N_pol*N_mt*N_mt_max]);
			        atomicAdd(&f[i].x, dfk*dxk);
			        atomicAdd(&f[i].y, dfk*dyk);		
			        atomicAdd(&f[i].z, dfk*dzk);

			        atomicAdd(&f[k].x, -dfk*dxk);
			        atomicAdd(&f[k].y, -dfk*dyk);		
			        atomicAdd(&f[k].z, -dfk*dzk);
                }
            }
        }
        if (type[i] == LEFT_NDC){
			float3 ri = r[i];
            int j, mult, ks;
            for (ks = 0; ks < N_kt; ks += 2){
                if (i < N_pol*N_mt*N_mt_max + (ks + 1)*kn && i >= N_pol*N_mt*N_mt_max + ks*kn){
                    j = N_pol*N_mt*N_mt_max + (ks + 1)*kn - 1 - chrom_bn;
                    mult = ks;
                }
            }
            float3 rj = r[j];

			float dx = rj.x - ri.x;
			float dy = rj.y - ri.y;
			float dz = rj.z - ri.z;

			float dr = sqrtf(dx*dx + dy*dy + dz*dz);
			float df = K_kt*(dr - kt_radius[i - N_pol*N_mt*N_mt_max])*powf(dr, -1);

			atomicAdd(&f[i].x, df*dx);
			atomicAdd(&f[i].y, df*dy);		
			atomicAdd(&f[i].z, df*dz);
            
            int l;          
            for (l = 0; l < harmonicKinCount[i - N_pol*N_mt*N_mt_max]; l++){
				float3 rl = r[N_pol*N_mt*N_mt_max + kn*mult + harmonicKin[(i - N_pol*N_mt*N_mt_max)*MpK + l]];

                float dxl = rl.x - ri.x;
				float dyl = rl.y - ri.y;
				float dzl = rl.z - ri.z;	
				float drl = sqrtf(dxl*dxl + dyl*dyl + dzl*dzl);
               
				float dfl = K_kt*(drl - harmonicKinRadii[(i - N_pol*N_mt*N_mt_max)*MpK + l])*powf(drl, -1)/100;
                             
                atomicAdd(&f[i].x, dfl*dxl);
			    atomicAdd(&f[i].y, dfl*dyl);		
			    atomicAdd(&f[i].z, dfl*dzl);

                atomicAdd(&f[N_pol*N_mt*N_mt_max + kn*mult + harmonicKin[(i - N_pol*N_mt*N_mt_max)*MpK + l]].x, -dfl*dxl);
			    atomicAdd(&f[N_pol*N_mt*N_mt_max + kn*mult + harmonicKin[(i - N_pol*N_mt*N_mt_max)*MpK + l]].y, -dfl*dyl);	
			    atomicAdd(&f[N_pol*N_mt*N_mt_max + kn*mult + harmonicKin[(i - N_pol*N_mt*N_mt_max)*MpK + l]].z, -dfl*dzl);
            }
        }
        if (type[i] == SHELL_LEFT){
			float3 ri = r[i];
            int j, mult, ks;
            for (ks = 0; ks < N_kt; ks += 2){
                if (i < N_pol*N_mt*N_mt_max + (ks + 1)*kn && i >= N_pol*N_mt*N_mt_max + ks*kn){
                    j = N_pol*N_mt*N_mt_max + (ks + 1)*kn - 1 - chrom_bn;
                    mult = ks;
                }
            }

            float3 rj = r[j];

			float dx = rj.x - ri.x;
			float dy = rj.y - ri.y;
			float dz = rj.z - ri.z;

			float dr = sqrtf(dx*dx + dy*dy + dz*dz);
			float df = K_kt*(dr - kt_radius[i - N_pol*N_mt*N_mt_max])*powf(dr, -1);

			atomicAdd(&f[i].x, df*dx);
			atomicAdd(&f[i].y, df*dy);		
			atomicAdd(&f[i].z, df*dz);
            
            int l;          
            for (l = 0; l < harmonicKinCount[i - N_pol*N_mt*N_mt_max]; l++){
				float3 rl = r[N_pol*N_mt*N_mt_max + kn*mult + harmonicKin[(i - N_pol*N_mt*N_mt_max)*MpK + l]];

                float dxl = rl.x - ri.x;
				float dyl = rl.y - ri.y;
				float dzl = rl.z - ri.z;	
				float drl = sqrtf(dxl*dxl + dyl*dyl + dzl*dzl);
               
				float dfl = K_kt*(drl - harmonicKinRadii[(i - N_pol*N_mt*N_mt_max)*MpK + l])*powf(drl, -1);
                             
                atomicAdd(&f[i].x, dfl*dxl);
			    atomicAdd(&f[i].y, dfl*dyl);		
			    atomicAdd(&f[i].z, dfl*dzl);

                atomicAdd(&f[N_pol*N_mt*N_mt_max + kn*mult + harmonicKin[(i - N_pol*N_mt*N_mt_max)*MpK + l]].x, -dfl*dxl);
			    atomicAdd(&f[N_pol*N_mt*N_mt_max + kn*mult + harmonicKin[(i - N_pol*N_mt*N_mt_max)*MpK + l]].y, -dfl*dyl);	
			    atomicAdd(&f[N_pol*N_mt*N_mt_max + kn*mult + harmonicKin[(i - N_pol*N_mt*N_mt_max)*MpK + l]].z, -dfl*dzl);
            }
        }
        if (type[i] == RIGHT_NDC){
			float3 ri = r[i];
            int j, mult, ks;
            for (ks = 1; ks < N_kt; ks += 2){
                if (i < N_pol*N_mt*N_mt_max + (ks + 1)*kn && i >= N_pol*N_mt*N_mt_max + ks*kn){
                    j = N_pol*N_mt*N_mt_max + (ks + 1)*kn - 1 - chrom_bn;
                    mult = ks;
                }
            }
            float3 rj = r[j];

			float dx = rj.x - ri.x;
			float dy = rj.y - ri.y;
			float dz = rj.z - ri.z;

			float dr = sqrtf(dx*dx + dy*dy + dz*dz);
			float df = K_kt*(dr - kt_radius[i - N_pol*N_mt*N_mt_max])*powf(dr, -1);
            //printf("%d\t%f\t%f\t%f\n", i, dr, kt_radius[i - N_pol*N_mt*N_mt_max], df);

			atomicAdd(&f[i].x, df*dx);
			atomicAdd(&f[i].y, df*dy);		
			atomicAdd(&f[i].z, df*dz);

            int l;
            for (l = 0; l < harmonicKinCount[i - N_pol*N_mt*N_mt_max]; l++){
			    float3 rl = r[N_pol*N_mt*N_mt_max + kn*mult + harmonicKin[(i - N_pol*N_mt*N_mt_max)*MpK + l]];
                float dxl = rl.x - ri.x;
			    float dyl = rl.y - ri.y;
			    float dzl = rl.z - ri.z;	
			    float drl = sqrtf(dxl*dxl + dyl*dyl + dzl*dzl);

			    float dfl = K_kt*(drl - harmonicKinRadii[(i - N_pol*N_mt*N_mt_max)*MpK + l])*powf(drl, -1)/100;
                //printf("%d\t%f\t%f\n", i, drl, harmonicKinRadii[(i - N_pol*N_mt*N_mt_max)*MpK + l]);                
                atomicAdd(&f[i].x, dfl*dxl);
		        atomicAdd(&f[i].y, dfl*dyl);		
		        atomicAdd(&f[i].z, dfl*dzl);

                atomicAdd(&f[N_pol*N_mt*N_mt_max + kn*mult + harmonicKin[(i - N_pol*N_mt*N_mt_max)*MpK + l]].x, -dfl*dxl);
		        atomicAdd(&f[N_pol*N_mt*N_mt_max + kn*mult + harmonicKin[(i - N_pol*N_mt*N_mt_max)*MpK + l]].y, -dfl*dyl);		
		        atomicAdd(&f[N_pol*N_mt*N_mt_max + kn*mult + harmonicKin[(i - N_pol*N_mt*N_mt_max)*MpK + l]].z, -dfl*dzl);
            }
        }
        if (type[i] == SHELL_RIGHT){
			float3 ri = r[i];
            int j, mult, ks;
            for (ks = 1; ks < N_kt; ks += 2){
                if (i < N_pol*N_mt*N_mt_max + (ks + 1)*kn && i >= N_pol*N_mt*N_mt_max + ks*kn){
                    j = N_pol*N_mt*N_mt_max + (ks + 1)*kn - 1 - chrom_bn;
                    mult = ks;
                }
            }
            float3 rj = r[j];

			float dx = rj.x - ri.x;
			float dy = rj.y - ri.y;
			float dz = rj.z - ri.z;

			float dr = sqrtf(dx*dx + dy*dy + dz*dz);
			float df = K_kt*(dr - kt_radius[i - N_pol*N_mt*N_mt_max])*powf(dr, -1);
            //printf("%d\t%f\t%f\t%f\n", i, dr, kt_radius[i - N_pol*N_mt*N_mt_max], df);

			atomicAdd(&f[i].x, df*dx);
			atomicAdd(&f[i].y, df*dy);		
			atomicAdd(&f[i].z, df*dz);

            int l;
            for (l = 0; l < harmonicKinCount[i - N_pol*N_mt*N_mt_max]; l++){
			    float3 rl = r[N_pol*N_mt*N_mt_max + kn*mult + harmonicKin[(i - N_pol*N_mt*N_mt_max)*MpK + l]];
                float dxl = rl.x - ri.x;
			    float dyl = rl.y - ri.y;
			    float dzl = rl.z - ri.z;	
			    float drl = sqrtf(dxl*dxl + dyl*dyl + dzl*dzl);

			    float dfl = K_kt*(drl - harmonicKinRadii[(i - N_pol*N_mt*N_mt_max)*MpK + l])*powf(drl, -1);
                //printf("%d\t%f\t%f\n", i, drl, harmonicKinRadii[(i - N_pol*N_mt*N_mt_max)*MpK + l]);                
                atomicAdd(&f[i].x, dfl*dxl);
		        atomicAdd(&f[i].y, dfl*dyl);		
		        atomicAdd(&f[i].z, dfl*dzl);

                atomicAdd(&f[N_pol*N_mt*N_mt_max + kn*mult + harmonicKin[(i - N_pol*N_mt*N_mt_max)*MpK + l]].x, -dfl*dxl);
		        atomicAdd(&f[N_pol*N_mt*N_mt_max + kn*mult + harmonicKin[(i - N_pol*N_mt*N_mt_max)*MpK + l]].y, -dfl*dyl);		
		        atomicAdd(&f[N_pol*N_mt*N_mt_max + kn*mult + harmonicKin[(i - N_pol*N_mt*N_mt_max)*MpK + l]].z, -dfl*dzl);
            }
        }
    //	__syncthreads();
	}
}
