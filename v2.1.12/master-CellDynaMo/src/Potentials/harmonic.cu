#include "harmonic.cuh"		

__global__ void computeHarmonic(float3* r, float3* f, int N, int* type, int kn, float rk, int* harmonicKinCount, int* harmonicKin, float* harmonicKinRadii, float* kt_radius, Param* d_parameters, float* length, int* connect1, int* connect2, float* len1, float* len2, float* f_out){
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
    float mt_rad = d_parameters->mt_r;
	int i;	
	i = blockIdx.x*blockDim.x + threadIdx.x;	
	if (i < N){
        if (type[i] == LIG_MT_ACTIVE){
            //printf("%d\t%d\t%d\n", i, type[i], LIG_MT_ACTIVE);
            float r0 = 50.0;
   			float3 ri = r[i];
			float3 rj = r[i + 1];

			float dx = rj.x - ri.x;
			float dy = rj.y - ri.y;
			float dz = rj.z - ri.z;
	        
			float dr = sqrtf(dx*dx + dy*dy + dz*dz);
            //printf("%d\t%d: %f\t%f\n", i, i - 1, dr, r0);
            if (dr > r0){
                float df = 10*K_mt*(dr - r0)*powf(dr, -1);
                f_out[i] = df;

	            atomicAdd(&f[i].x, df*dx);
	            atomicAdd(&f[i].y, df*dy);		
	            atomicAdd(&f[i].z, df*dz);
                
		        atomicAdd(&f[i + 1].x, -df*dx);
		        atomicAdd(&f[i + 1].y, -df*dy);		
		        atomicAdd(&f[i + 1].z, -df*dz);
            }
            int m1 = connect1[i];
            
            float3 r_m1 = r[m1];
            float dr1 = len1[i];

		    float dx1 = r_m1.x - ri.x;
		    float dy1 = r_m1.y - ri.y;
		    float dz1 = r_m1.z - ri.z;

            float dr_m1 = sqrtf(dx1*dx1 + dy1*dy1 + dz1*dz1);
            //printf("1: %d\t%d\t: %f\t%f\n", i, m1, dr_m1, dr1);
            //printf("%d\t%d: %f\t%f\n", i, i - 1, dr, r0);
            //if (dr_m1 > 2*mt_rad){
            float df1 = 100*K_mt*(dr_m1 - dr1)*powf(dr_m1, -1);
            if (df1 > 100){
                df1 = 100.0;
            }
            if (df1 < -100){
                df1 = -100.0;
            }

            atomicAdd(&f[i].x, df1*dx1);
            atomicAdd(&f[i].y, df1*dy1);		
            atomicAdd(&f[i].z, df1*dz1);
            
            atomicAdd(&f[m1].x, -df1*dx1);
            atomicAdd(&f[m1].y, -df1*dy1);		
            atomicAdd(&f[m1].z, -df1*dz1);

            int m2 = connect2[i];            
            float3 r_m2 = r[m2];
            float dr2 = len2[i];

            float dx2 = r_m2.x - ri.x;
		    float dy2 = r_m2.y - ri.y;
		    float dz2 = r_m2.z - ri.z;

            float dr_m2 = sqrtf(dx2*dx2 + dy2*dy2 + dz2*dz2);
            //printf("2: %d\t%d\t: %f\t%f\n", i, m2, dr_m2, dr2);
            //printf("%d\t%d: %f\t%f\n", i, i - 1, dr, r0);
            //if (dr_m2 > 2*mt_rad){
            float df2 = 100*K_mt*(dr_m2 - dr2)*powf(dr_m2, -1);
            if (df2 > 100){
                df2 = 100.0;
            }
            if (df2 < -100){
                df2 = -100.0;
            }

            atomicAdd(&f[i].x, df2*dx2);
            atomicAdd(&f[i].y, df2*dy2);		
            atomicAdd(&f[i].z, df2*dz2);
            
            atomicAdd(&f[m2].x, -df2*dx2);
            atomicAdd(&f[m2].y, -df2*dy2);		
            atomicAdd(&f[m2].z, -df2*dz2);
            //}

            float3 ri1 = r[connect1[i]];  //a0
		    float3 rj1 = r[connect2[i]];  //a1
            
            float3 intersect;
            float min_dist;
            float dxc, dyc, dzc;

            float3 p_vector1;
            p_vector1.x = rj1.x - ri1.x;  //AB
            p_vector1.y = rj1.y - ri1.y;
            p_vector1.z = rj1.z - ri1.z;

            float dp1 = sqrtf(p_vector1.x*p_vector1.x + p_vector1.y*p_vector1.y + p_vector1.z+ p_vector1.z); //magA

            float3 p_vectorA;
            p_vectorA.x = ri.x - ri1.x;  //AP
            p_vectorA.y = ri.y - ri1.y;
            p_vectorA.z = ri.z - ri1.z;

            float3 p_vectorB;
            p_vectorB.x = ri.x - rj1.x;  //BP
            p_vectorB.y = ri.y - rj1.y;
            p_vectorB.z = ri.z - rj1.z;

            float AB_BP_dot = (p_vector1.x*p_vectorB.x + p_vector1.y*p_vectorB.y + p_vector1.z*p_vectorB.z);
            float AB_AP_dot = (p_vector1.x*p_vectorA.x + p_vector1.y*p_vectorA.y + p_vector1.z*p_vectorA.z);

            if (AB_BP_dot > 0){
                dxc = ri.x - rj1.x;
                dyc = ri.y - rj1.y;
                dzc = ri.z - rj1.z;
                min_dist = sqrtf(dxc*dxc + dyc*dyc + dzc*dzc);
                intersect = rj1;
            }
            else    if (AB_AP_dot < 0){
                dxc = ri.x - ri1.x;
                dyc = ri.y - ri1.y;
                dzc = ri.z - ri1.z;
                min_dist = sqrtf(dxc*dxc + dyc*dyc + dzc*dzc);
                intersect = ri1;
            }
            else    {
                float3 cross;
                cross.x = p_vectorA.y*p_vectorB.z - p_vectorB.y*p_vectorA.z; //cross
                cross.y = p_vectorA.z*p_vectorB.x - p_vectorB.z*p_vectorA.x;
                cross.z = p_vectorA.x*p_vectorB.y - p_vectorB.x*p_vectorA.y;

                float d_cross = sqrtf(cross.x*cross.x + cross.y*cross.y + cross.z*cross.z); //denom
                min_dist = d_cross/dp1;

                float t = AB_AP_dot/powf(dp1, 2);
                intersect.x = ri1.x + p_vector1.x*t;
                intersect.y = ri1.y + p_vector1.y*t;
                intersect.z = ri1.z + p_vector1.z*t;
                dxc = ri.x - intersect.x;
                dyc = ri.y - intersect.y;
                dzc = ri.z - intersect.z;
            }

            float3 d_a;
            d_a.x = intersect.x - ri1.x;
            d_a.y = intersect.y - ri1.y;
            d_a.z = intersect.z - ri1.z;
            float alpha = sqrtf(d_a.x*d_a.x + d_a.y* d_a.y + d_a.z*d_a.z)/dp1; 
            //printf("MT: %d\t%d\t%d\t%f\t%f\n", i, connect1[i], connect2[i], min_dist, mt_rad);
            float dfc = 10*K_mt*(min_dist - mt_rad)*powf(min_dist, -1);
            if (dfc > 100.0){
                dfc = 100.0;
            } 
            if (dfc < -100.0){
                dfc = -100.0;
            } 

            /*atomicAdd(&f[connect1[i]].x, dfc*(1 - alpha)*dxc);
            atomicAdd(&f[connect1[i]].y, dfc*(1 - alpha)*dyc);  		
            atomicAdd(&f[connect1[i]].z, dfc*(1 - alpha)*dzc);

            atomicAdd(&f[connect2[i]].x, dfc*alpha*dxc);
            atomicAdd(&f[connect2[i]].y, dfc*alpha*dyc);  		
            atomicAdd(&f[connect2[i]].z, dfc*alpha*dzc);

            atomicAdd(&f[i].x, -dfc*dxc);
            atomicAdd(&f[i].y, -dfc*dyc);		
            atomicAdd(&f[i].z, -dfc*dzc);*/ 
        }

        if (type[i] == LIG_CH_ACTIVE){
   			float3 ri = r[i];
			float3 rj = r[i + 1];
            int m1 = connect1[i];
           
            float3 r_m1 = r[m1];
            float dr1 = len1[i];

		    float dx1 = r_m1.x - ri.x;
		    float dy1 = r_m1.y - ri.y;
		    float dz1 = r_m1.z - ri.z;

            float dr_m1 = sqrtf(dx1*dx1 + dy1*dy1 + dz1*dz1);

            float df1 = 10*K_mt*(dr_m1 - dr1)*powf(dr_m1, -1);

            atomicAdd(&f[i].x, df1*dx1);
            atomicAdd(&f[i].y, df1*dy1);		
            atomicAdd(&f[i].z, df1*dz1);
            
            atomicAdd(&f[m1].x, -df1*dx1);
            atomicAdd(&f[m1].y, -df1*dy1);		
            atomicAdd(&f[m1].z, -df1*dz1);

            int m2 = connect2[i];            
            float3 r_m2 = r[m2];
            float dr2 = len2[i];
            
            float dx2 = r_m2.x - ri.x;
		    float dy2 = r_m2.y - ri.y;
		    float dz2 = r_m2.z - ri.z;

            float dr_m2 = sqrtf(dx2*dx2 + dy2*dy2 + dz2*dz2);
            //printf("%d\t%d: %f\t%f\n", i, i - 1, dr, r0);

            float df2 = 10*K_mt*(dr_m2 - dr2)*powf(dr_m2, -1);

            atomicAdd(&f[i].x, df2*dx2);
            atomicAdd(&f[i].y, df2*dy2);		
            atomicAdd(&f[i].z, df2*dz2);
            
            atomicAdd(&f[m2].x, -df2*dx2);
            atomicAdd(&f[m2].y, -df2*dy2);		
            atomicAdd(&f[m2].z, -df2*dz2);

            float3 ri1 = r[connect1[i]];  //a0
		    float3 rj1 = r[connect2[i]];  //a1
            
            float3 intersect;
            float min_dist;
            float dxc, dyc, dzc;

            float3 p_vector1;
            p_vector1.x = rj1.x - ri1.x;  //AB
            p_vector1.y = rj1.y - ri1.y;
            p_vector1.z = rj1.z - ri1.z;

            float dp1 = sqrtf(p_vector1.x*p_vector1.x + p_vector1.y*p_vector1.y + p_vector1.z+ p_vector1.z); //magA

            float3 p_vectorA;
            p_vectorA.x = ri.x - ri1.x;  //AP
            p_vectorA.y = ri.y - ri1.y;
            p_vectorA.z = ri.z - ri1.z;

            float3 p_vectorB;
            p_vectorB.x = ri.x - rj1.x;  //BP
            p_vectorB.y = ri.y - rj1.y;
            p_vectorB.z = ri.z - rj1.z;

            float AB_BP_dot = (p_vector1.x*p_vectorB.x + p_vector1.y*p_vectorB.y + p_vector1.z*p_vectorB.z);
            float AB_AP_dot = (p_vector1.x*p_vectorA.x + p_vector1.y*p_vectorA.y + p_vector1.z*p_vectorA.z);

            if (AB_BP_dot > 0){
                dxc = ri.x - rj1.x;
                dyc = ri.y - rj1.y;
                dzc = ri.z - rj1.z;
                min_dist = sqrtf(dxc*dxc + dyc*dyc + dzc*dzc);
                intersect = rj1;
            }
            else    if (AB_AP_dot < 0){
                dxc = ri.x - ri1.x;
                dyc = ri.y - ri1.y;
                dzc = ri.z - ri1.z;
                min_dist = sqrtf(dxc*dxc + dyc*dyc + dzc*dzc);
                intersect = ri1;
            }
            else    {
                float3 cross;
                cross.x = p_vectorA.y*p_vectorB.z - p_vectorB.y*p_vectorA.z; //cross
                cross.y = p_vectorA.z*p_vectorB.x - p_vectorB.z*p_vectorA.x;
                cross.z = p_vectorA.x*p_vectorB.y - p_vectorB.x*p_vectorA.y;

                float d_cross = sqrtf(cross.x*cross.x + cross.y*cross.y + cross.z*cross.z); //denom
                min_dist = d_cross/dp1;

                float t = AB_AP_dot/powf(dp1, 2);
                intersect.x = ri1.x + p_vector1.x*t;
                intersect.y = ri1.y + p_vector1.y*t;
                intersect.z = ri1.z + p_vector1.z*t;
                dxc = ri.x - intersect.x;
                dyc = ri.y - intersect.y;
                dzc = ri.z - intersect.z;
            }

            float3 d_a;
            d_a.x = intersect.x - ri1.x;
            d_a.y = intersect.y - ri1.y;
            d_a.z = intersect.z - ri1.z;
            float alpha = sqrtf(d_a.x*d_a.x + d_a.y* d_a.y + d_a.z*d_a.z)/dp1; 
            //printf("CH: %d\t%d\t%d\t%f\t%f\n", i, connect1[i], connect2[i], min_dist, kt_rad);
            float dfc = 10*K_mt*(min_dist - kt_rad)*powf(min_dist, -1);
            if (dfc > 100.0){
                dfc = 100.0;
            } 
            if (dfc < -100.0){
                dfc = -100.0;
            } 

            atomicAdd(&f[connect1[i]].x, dfc*(1 - alpha)*dxc);
            atomicAdd(&f[connect1[i]].y, dfc*(1 - alpha)*dyc);  		
            atomicAdd(&f[connect1[i]].z, dfc*(1 - alpha)*dzc);

            atomicAdd(&f[connect2[i]].x, dfc*alpha*dxc);
            atomicAdd(&f[connect2[i]].y, dfc*alpha*dyc);  		
            atomicAdd(&f[connect2[i]].z, dfc*alpha*dzc);

            atomicAdd(&f[i].x, -dfc*dxc);
            atomicAdd(&f[i].y, -dfc*dyc);		
            atomicAdd(&f[i].z, -dfc*dzc);
        }
		if (type[i] == MT_REG || type[i] == PLUS_DET || type[i] == PLUS_ATT || type[i] == PLUS_DET_INVALID){
            float r0 = length[i]/(N_mt_max - 1);
            if (r0 > 12*mt_rad){
       			float3 ri = r[i - 1];
			    float3 rj = r[i];

			    float dx = rj.x - ri.x;
			    float dy = rj.y - ri.y;
			    float dz = rj.z - ri.z;
	            
			    float dr = sqrtf(dx*dx + dy*dy + dz*dz);
                //printf("%d\t%d: %f\t%f\n", i, i - 1, dr, r0);
                float df = 100*K_mt*(dr - r0)*powf(dr, -1);

		        atomicAdd(&f[i - 1].x, df*dx);
		        atomicAdd(&f[i - 1].y, df*dy);		
		        atomicAdd(&f[i - 1].z, df*dz);
                
			    atomicAdd(&f[i].x, -df*dx);
			    atomicAdd(&f[i].y, -df*dy);		
			    atomicAdd(&f[i].z, -df*dz);
            }		
        }
        if (type[i] == LEFT_MINUS){
            float r0 = length[i]/(N_mt_max - 1);
            if (r0 > 12*mt_rad){
       			float3 ri = r[N - 2];
			    float3 rj = r[i];

			    float dx = rj.x - ri.x;
			    float dy = rj.y - ri.y;
			    float dz = rj.z - ri.z;
	        
			    float dr = sqrtf(dx*dx + dy*dy + dz*dz);
                float df = 100*K_mt*(dr - r0)*powf(dr, -1);

		        atomicAdd(&f[N - 2].x, df*dx);
		        atomicAdd(&f[N - 2].y, df*dy);		
		        atomicAdd(&f[N - 2].z, df*dz);
                
			    atomicAdd(&f[i].x, -df*dx);
			    atomicAdd(&f[i].y, -df*dy);		
			    atomicAdd(&f[i].z, -df*dz);
            }
		}
        if (type[i] == RIGHT_MINUS){
            float r0 = length[i]/(N_mt_max - 1);
            if (r0 > 12*mt_rad){
       			float3 ri = r[N - 1];
			    float3 rj = r[i];

			    float dx = rj.x - ri.x;
			    float dy = rj.y - ri.y;
			    float dz = rj.z - ri.z;
	        
			    float dr = sqrtf(dx*dx + dy*dy + dz*dz);
                float df = 100*K_mt*(dr - r0)*powf(dr, -1);

		        atomicAdd(&f[N - 1].x, df*dx);
		        atomicAdd(&f[N - 1].y, df*dy);		
		        atomicAdd(&f[N - 1].z, df*dz);
                
			    atomicAdd(&f[i].x, -df*dx);
			    atomicAdd(&f[i].y, -df*dy);		
			    atomicAdd(&f[i].z, -df*dz);
            }
		}
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
			float df = 5.0*K_kt*(dr - rk)*powf(dr, -1);

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
                    //printf("LKT: %d\t%d\t%f\t%f\n", i, k, drk, kt_radius[k - N_pol*N_mt*N_mt_max]);

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
               // printf("RKT: %d\t%d\t%f\t%f\n", i, k, drk, kt_radius[k - N_pol*N_mt*N_mt_max]);

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
                   // printf("CHROM: %d\t%d\t%f\t%f\n", i, k, drk, kt_radius[k - N_pol*N_mt*N_mt_max]);
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
        /*if (type[i] == SHELL_LEFT){
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
        }*/
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
        /*if (type[i] == SHELL_RIGHT){
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
        }*/
    //	__syncthreads();
	}
}
