#include "bending.cuh"		

__global__ void computeAngles(float3* r, float3* f, int N, int* type, int kn, float* kt_cos, float* dyn_cos, Param* d_parameters, float* length){
	int i;
    float temper = d_parameters->temper;
    float pers_l = d_parameters->lp;
    int N_pol = d_parameters->pole_n;
    int N_mt = d_parameters->mt_n;
    int N_kt = d_parameters->kin_n;
    int N_mt_max = d_parameters->max_mt_n;
    float K_kt = d_parameters->k_kt;
    int chrom_cond = d_parameters->chrom;
    int chrom_bn = chrom_cond*d_parameters->chrom_num;
	i = blockIdx.x*blockDim.x + threadIdx.x;
	if (i < N){
		if (type[i] == LEFT_MINUS){
            float k_theta = 10.0*KB*temper*pers_l/1000.0;
            float theta0 = M_PI;
            int j = i + 1;
            int k = N - 2;
			float3 r1 = r[j];
			float3 r2 = r[i];
			float3 r3 = r[k];

			float3 dr12, dr32;
			dr12.x = r1.x - r2.x;
			dr12.y = r1.y - r2.y;
			dr12.z = r1.z - r2.z;

			dr32.x = r3.x - r2.x;
			dr32.y = r3.y - r2.y;
			dr32.z = r3.z - r2.z;

			float r12inv = 1.0f*powf(dr12.x*dr12.x + dr12.y*dr12.y + dr12.z*dr12.z, -0.5);
			float r32inv = 1.0f*powf(dr32.x*dr32.x + dr32.y*dr32.y + dr32.z*dr32.z, -0.5);
			float costheta = (dr12.x*dr32.x + dr12.y*dr32.y + dr12.z*dr32.z)*r12inv*r32inv;
			if (costheta > 1.0f){
				costheta = 1.0f;
			}
			else	if (costheta < -1.0f){
				costheta = -1.0f;
			}
			float sintheta = sqrtf(1.0f - costheta*costheta);
			float theta = acos(costheta);
			float diff = theta - theta0;
   			if (sintheta < 1.e-6){
				if (diff < 0){
					diff *= 2.0f*k_theta;
				}
				else{
					diff *= -2.0f*k_theta;
				}
			}
			else{
				diff *= (-2.0f*k_theta)*powf(sintheta, -1);
			}
			float c1 = diff*r12inv;
			float c2 = diff*r32inv;

			float3 f1, f2, f3;
			f1.x = c1*(dr12.x*(r12inv*costheta) - dr32.x*r32inv);
			f1.y = c1*(dr12.y*(r12inv*costheta) - dr32.y*r32inv);
			f1.z = c1*(dr12.z*(r12inv*costheta) - dr32.z*r32inv);
			f2 = f1;
			f3.x = c2*(dr32.x*(r32inv*costheta) - dr12.x*r12inv);
			f3.y = c2*(dr32.y*(r32inv*costheta) - dr12.y*r12inv);
			f3.z = c2*(dr32.z*(r32inv*costheta) - dr12.z*r12inv);
			f2.x += f3.x;
			f2.y += f3.y;
			f2.z += f3.z;

			f2.x = -f2.x;
			f2.y = -f2.y;
			f2.z = -f2.z;

            atomicAdd(&f[j].x, f1.x);
		    atomicAdd(&f[j].y, f1.y);		
		    atomicAdd(&f[j].z, f1.z);

            /*atomicAdd(&f[i].x, f2.x);
			atomicAdd(&f[i].y, f2.y);		
			atomicAdd(&f[i].z, f2.z);*/
           
            atomicAdd(&f[k].x, f3.x);
		    atomicAdd(&f[k].y, f3.y);		
		    atomicAdd(&f[k].z, f3.z);       
		}
        if (type[i] == RIGHT_MINUS){
            float k_theta = 10.0*KB*temper*pers_l/1000.0;
            float theta0 = M_PI;
            int j = i + 1;
            int k = N - 1;
			float3 r1 = r[j];
			float3 r2 = r[i];
			float3 r3 = r[k];

			float3 dr12, dr32;
			dr12.x = r1.x - r2.x;
			dr12.y = r1.y - r2.y;
			dr12.z = r1.z - r2.z;

			dr32.x = r3.x - r2.x;
			dr32.y = r3.y - r2.y;
			dr32.z = r3.z - r2.z;

			float r12inv = 1.0f*powf(dr12.x*dr12.x + dr12.y*dr12.y + dr12.z*dr12.z, -0.5);
			float r32inv = 1.0f*powf(dr32.x*dr32.x + dr32.y*dr32.y + dr32.z*dr32.z, -0.5);
			float costheta = (dr12.x*dr32.x + dr12.y*dr32.y + dr12.z*dr32.z)*r12inv*r32inv;
			if (costheta > 1.0f){
				costheta = 1.0f;
			}
			else	if (costheta < -1.0f){
				costheta = -1.0f;
			}
			float sintheta = sqrtf(1.0f - costheta*costheta);
			float theta = acos(costheta);
			float diff = theta - theta0;
   			if (sintheta < 1.e-6){
				if (diff < 0){
					diff *= 2.0f*k_theta;
				}
				else{
					diff *= -2.0f*k_theta;
				}
			}
			else{
				diff *= (-2.0f*k_theta)*powf(sintheta, -1);
			}
			float c1 = diff*r12inv;
			float c2 = diff*r32inv;

			float3 f1, f2, f3;
			f1.x = c1*(dr12.x*(r12inv*costheta) - dr32.x*r32inv);
			f1.y = c1*(dr12.y*(r12inv*costheta) - dr32.y*r32inv);
			f1.z = c1*(dr12.z*(r12inv*costheta) - dr32.z*r32inv);
			f2 = f1;
			f3.x = c2*(dr32.x*(r32inv*costheta) - dr12.x*r12inv);
			f3.y = c2*(dr32.y*(r32inv*costheta) - dr12.y*r12inv);
			f3.z = c2*(dr32.z*(r32inv*costheta) - dr12.z*r12inv);
			f2.x += f3.x;
			f2.y += f3.y;
			f2.z += f3.z;

			f2.x = -f2.x;
			f2.y = -f2.y;
			f2.z = -f2.z;

            atomicAdd(&f[j].x, f1.x);
		    atomicAdd(&f[j].y, f1.y);		
		    atomicAdd(&f[j].z, f1.z);

            /*atomicAdd(&f[i].x, f2.x);
			atomicAdd(&f[i].y, f2.y);		
			atomicAdd(&f[i].z, f2.z);*/
           
            atomicAdd(&f[k].x, f3.x);
		    atomicAdd(&f[k].y, f3.y);		
		    atomicAdd(&f[k].z, f3.z);       
		}
        if (type[i] == MT_REG || type[i] == MT_END){
            float k_theta = 10.0*KB*temper*pers_l/1000.0;
            float theta0 = M_PI;
            int j = i + 1;
            int k = i - 1;
			float3 r1 = r[j];
			float3 r2 = r[i];
			float3 r3 = r[k];

			float3 dr12, dr32;
			dr12.x = r1.x - r2.x;
			dr12.y = r1.y - r2.y;
			dr12.z = r1.z - r2.z;

			dr32.x = r3.x - r2.x;
			dr32.y = r3.y - r2.y;
			dr32.z = r3.z - r2.z;

			float r12inv = 1.0f*powf(dr12.x*dr12.x + dr12.y*dr12.y + dr12.z*dr12.z, -0.5);
			float r32inv = 1.0f*powf(dr32.x*dr32.x + dr32.y*dr32.y + dr32.z*dr32.z, -0.5);
			float costheta = (dr12.x*dr32.x + dr12.y*dr32.y + dr12.z*dr32.z)*r12inv*r32inv;
			if (costheta > 1.0f){
				costheta = 1.0f;
			}
			else	if (costheta < -1.0f){
				costheta = -1.0f;
			}
			float sintheta = sqrtf(1.0f - costheta*costheta);
			float theta = acos(costheta);
			float diff = theta - theta0;
   			if (sintheta < 1.e-6){
				if (diff < 0){
					diff *= 2.0f*k_theta;
				}
				else{
					diff *= -2.0f*k_theta;
				}
			}
			else{
				diff *= (-2.0f*k_theta)*powf(sintheta, -1);
			}
			float c1 = diff*r12inv;
			float c2 = diff*r32inv;

			float3 f1, f2, f3;
			f1.x = c1*(dr12.x*(r12inv*costheta) - dr32.x*r32inv);
			f1.y = c1*(dr12.y*(r12inv*costheta) - dr32.y*r32inv);
			f1.z = c1*(dr12.z*(r12inv*costheta) - dr32.z*r32inv);
			f2 = f1;
			f3.x = c2*(dr32.x*(r32inv*costheta) - dr12.x*r12inv);
			f3.y = c2*(dr32.y*(r32inv*costheta) - dr12.y*r12inv);
			f3.z = c2*(dr32.z*(r32inv*costheta) - dr12.z*r12inv);
			f2.x += f3.x;
			f2.y += f3.y;
			f2.z += f3.z;

			f2.x = -f2.x;
			f2.y = -f2.y;
			f2.z = -f2.z;

            atomicAdd(&f[j].x, f1.x);
		    atomicAdd(&f[j].y, f1.y);		
		    atomicAdd(&f[j].z, f1.z);

            atomicAdd(&f[i].x, f2.x);
			atomicAdd(&f[i].y, f2.y);		
			atomicAdd(&f[i].z, f2.z);
            
            if (type[k] != RIGHT_MINUS && type[k] != LEFT_MINUS){
                atomicAdd(&f[k].x, f3.x);
		        atomicAdd(&f[k].y, f3.y);		
		        atomicAdd(&f[k].z, f3.z);       
            }
		}


        if ((type[i] == LEFT_KT || type[i] == RIGHT_KT) && chrom_bn != 0){
            float k_theta = K_kt;
            for (int chrom = 1; chrom < chrom_bn; chrom += chrom_bn/2){
                int j = i + chrom;
                int k = j + 1;

			    float3 r1 = r[k];
			    float3 r2 = r[j];
			    float3 r3 = r[i];

			    float3 dr12, dr32;
			    dr12.x = r1.x - r2.x;
			    dr12.y = r1.y - r2.y;
			    dr12.z = r1.z - r2.z;

			    dr32.x = r3.x - r2.x;
			    dr32.y = r3.y - r2.y;
			    dr32.z = r3.z - r2.z;

			    float r12inv = 1.0f*powf(dr12.x*dr12.x + dr12.y*dr12.y + dr12.z*dr12.z, -0.5);
			    float r32inv = 1.0f*powf(dr32.x*dr32.x + dr32.y*dr32.y + dr32.z*dr32.z, -0.5);
			    float costheta = (dr12.x*dr32.x + dr12.y*dr32.y + dr12.z*dr32.z)*r12inv*r32inv;
			    if (costheta > 1.0f){
				    costheta = 1.0f;
			    }
			    else	if (costheta < -1.0f){
				    costheta = -1.0f;
			    }
			    float sintheta = sqrtf(1.0f - costheta*costheta);
			    float theta = acos(costheta);
			    float diff = theta - kt_cos[j - N_pol*N_mt*N_mt_max];

                dyn_cos[k - N_pol*N_mt*N_mt_max] = theta*180/M_PI;
                //printf ("BETA: %d\t%f\t%f\n", k, theta*180/M_PI, kt_cos[k - N_pol*N_mt*N_mt_max]*180/M_PI);
       			if (sintheta < 1.e-6){
				    if (diff < 0){
					    diff *= 2.0f*357140;
				    }
				    else{
					    diff *= -2.0f*357140;
				    }
			    }
			    else{
				    diff *= (-2.0f*357140)*powf(sintheta, -1);
			    }
			    float c1 = diff*r12inv;
			    float c2 = diff*r32inv;

			    float3 f1, f2, f3;
			    f1.x = c1*(dr12.x*(r12inv*costheta) - dr32.x*r32inv);
			    f1.y = c1*(dr12.y*(r12inv*costheta) - dr32.y*r32inv);
			    f1.z = c1*(dr12.z*(r12inv*costheta) - dr32.z*r32inv);
			    f2 = f1;
			    f3.x = c2*(dr32.x*(r32inv*costheta) - dr12.x*r12inv);
			    f3.y = c2*(dr32.y*(r32inv*costheta) - dr12.y*r12inv);
			    f3.z = c2*(dr32.z*(r32inv*costheta) - dr12.z*r12inv);
			    f2.x += f3.x;
			    f2.y += f3.y;
			    f2.z += f3.z;

			    f2.x = -f2.x;
			    f2.y = -f2.y;
			    f2.z = -f2.z;

                atomicAdd(&f[k].x, f1.x);
		        atomicAdd(&f[k].y, f1.y);		
		        atomicAdd(&f[k].z, f1.z);

                atomicAdd(&f[j].x, f2.x);
			    atomicAdd(&f[j].y, f2.y);		
			    atomicAdd(&f[j].z, f2.z);
               
                atomicAdd(&f[i].x, f3.x);
		        atomicAdd(&f[i].y, f3.y);		
		        atomicAdd(&f[i].z, f3.z);
            }

            int j = i + 1;
            int k = i + chrom_bn/2 + 1;

		    float3 r1 = r[k];
		    float3 r2 = r[i];
		    float3 r3 = r[j];

		    float3 dr12, dr32;
		    dr12.x = r1.x - r2.x;
		    dr12.y = r1.y - r2.y;
		    dr12.z = r1.z - r2.z;

		    dr32.x = r3.x - r2.x;
		    dr32.y = r3.y - r2.y;
		    dr32.z = r3.z - r2.z;

		    float r12inv = 1.0f*powf(dr12.x*dr12.x + dr12.y*dr12.y + dr12.z*dr12.z, -0.5);
		    float r32inv = 1.0f*powf(dr32.x*dr32.x + dr32.y*dr32.y + dr32.z*dr32.z, -0.5);
		    float costheta = (dr12.x*dr32.x + dr12.y*dr32.y + dr12.z*dr32.z)*r12inv*r32inv;
		    if (costheta > 1.0f){
			    costheta = 1.0f;
		    }
		    else	if (costheta < -1.0f){
			    costheta = -1.0f;
		    }
		    float sintheta = sqrtf(1.0f - costheta*costheta);
            float new_sin = (dr12.x*dr32.y - dr12.y*dr32.x)*r12inv*r32inv;
		    float theta = acos(costheta);
		    float diff = theta - kt_cos[j - N_pol*N_mt*N_mt_max];
            if (new_sin >= 0.0){
                dyn_cos[j - N_pol*N_mt*N_mt_max] = theta*180/M_PI;
            }
            if (new_sin < 0.0){
                float val = 360.0 - theta*180/M_PI;
                dyn_cos[j - N_pol*N_mt*N_mt_max] = val;
            }                        
            //printf ("ALPHA: %d\t%f\t%f\n", j, theta*180/M_PI, kt_cos[j - N_pol*N_mt*N_mt_max]*180/M_PI);
   			if (sintheta < 1.e-6){
			    if (diff < 0){
				    diff *= 2.0f*357140;
			    }
			    else{
				    diff *= -2.0f*357140;
			    }
		    }
		    else{
			    diff *= (-2.0f*357140)*powf(sintheta, -1);
		    }
		    float c1 = diff*r12inv;
		    float c2 = diff*r32inv;

		    float3 f1, f2, f3;
		    f1.x = c1*(dr12.x*(r12inv*costheta) - dr32.x*r32inv);
		    f1.y = c1*(dr12.y*(r12inv*costheta) - dr32.y*r32inv);
		    f1.z = c1*(dr12.z*(r12inv*costheta) - dr32.z*r32inv);
		    f2 = f1;
		    f3.x = c2*(dr32.x*(r32inv*costheta) - dr12.x*r12inv);
		    f3.y = c2*(dr32.y*(r32inv*costheta) - dr12.y*r12inv);
		    f3.z = c2*(dr32.z*(r32inv*costheta) - dr12.z*r12inv);
		    f2.x += f3.x;
		    f2.y += f3.y;
		    f2.z += f3.z;

		    f2.x = -f2.x;
		    f2.y = -f2.y;
		    f2.z = -f2.z;

            atomicAdd(&f[k].x, f1.x);
	        atomicAdd(&f[k].y, f1.y);		
	        atomicAdd(&f[k].z, f1.z);

            atomicAdd(&f[i].x, f2.x);
		    atomicAdd(&f[i].y, f2.y);		
		    atomicAdd(&f[i].z, f2.z);
           
            atomicAdd(&f[j].x, f3.x);
	        atomicAdd(&f[j].y, f3.y);		
	        atomicAdd(&f[j].z, f3.z);
            //printf("%d\t%d\t%d\t%f\t%f\n", i, j, k, theta, kt_cos[j - N_pol*N_mt*N_mt_max]); */
        }
        if (type[i] == CHROM){
            int j, k, ks;
            for (ks = 0; ks < N_kt; ks++){
                if (i < N_pol*N_mt*N_mt_max + (ks + 1)*kn - chrom_bn/2 && i > N_pol*N_mt*N_mt_max + (ks + 1)*kn - chrom_bn + 1 || i < N_pol*N_mt*N_mt_max + (ks + 1)*kn && i > N_pol*N_mt*N_mt_max + (ks + 1)*kn - chrom_bn/2 + 1){
                    j = i - 1;
                    k = j - 1;
                    float3 r1 = r[k];
		            float3 r2 = r[j];
		            float3 r3 = r[i];
                  
		            float3 dr12, dr32;
		            dr12.x = r1.x - r2.x;
		            dr12.y = r1.y - r2.y;
		            dr12.z = r1.z - r2.z;

		            dr32.x = r3.x - r2.x;
		            dr32.y = r3.y - r2.y;
		            dr32.z = r3.z - r2.z;

		            float r12inv = 1.0f*powf(dr12.x*dr12.x + dr12.y*dr12.y + dr12.z*dr12.z, -0.5);
		            float r32inv = 1.0f*powf(dr32.x*dr32.x + dr32.y*dr32.y + dr32.z*dr32.z, -0.5);
		            float costheta = (dr12.x*dr32.x + dr12.y*dr32.y + dr12.z*dr32.z)*r12inv*r32inv;
		            if (costheta > 1.0f){
			            costheta = 1.0f;
		            }
		            else	if (costheta < -1.0f){
			            costheta = -1.0f;
		            }
		            float sintheta = sqrtf(1.0f - costheta*costheta);
		            float theta = acos(costheta);
		            float diff = theta - kt_cos[i - N_pol*N_mt*N_mt_max];
                    dyn_cos[i - N_pol*N_mt*N_mt_max] = theta*180/M_PI;
                    //printf ("BETA: %d\t%f\t%f\n", k, theta*180/M_PI, kt_cos[k - N_pol*N_mt*N_mt_max]*180/M_PI);
           			if (sintheta < 1.e-6){
			            if (diff < 0){
				            diff *= 2.0f*357140;
			            }
			            else{
				            diff *= -2.0f*357140;
			            }
		            }
		            else{
			            diff *= (-2.0f*357140)*powf(sintheta, -1);
		            }
		            float c1 = diff*r12inv;
		            float c2 = diff*r32inv;

		            float3 f1, f2, f3;
		            f1.x = c1*(dr12.x*(r12inv*costheta) - dr32.x*r32inv);
		            f1.y = c1*(dr12.y*(r12inv*costheta) - dr32.y*r32inv);
		            f1.z = c1*(dr12.z*(r12inv*costheta) - dr32.z*r32inv);
		            f2 = f1;
		            f3.x = c2*(dr32.x*(r32inv*costheta) - dr12.x*r12inv);
		            f3.y = c2*(dr32.y*(r32inv*costheta) - dr12.y*r12inv);
		            f3.z = c2*(dr32.z*(r32inv*costheta) - dr12.z*r12inv);
		            f2.x += f3.x;
		            f2.y += f3.y;
		            f2.z += f3.z;

		            f2.x = -f2.x;
		            f2.y = -f2.y;
		            f2.z = -f2.z;

                    atomicAdd(&f[k].x, f1.x);
	                atomicAdd(&f[k].y, f1.y);		
	                atomicAdd(&f[k].z, f1.z);

                    atomicAdd(&f[j].x, f2.x);
		            atomicAdd(&f[j].y, f2.y);		
		            atomicAdd(&f[j].z, f2.z);
                   
                    atomicAdd(&f[i].x, f3.x);
	                atomicAdd(&f[i].y, f3.y);		
	                atomicAdd(&f[i].z, f3.z);
                }
            }
        }

        if (type[i] == LEFT_NDC){
            int j, k, ks;
            for (ks = 0; ks < N_kt; ks += 2){
                if (i < N_pol*N_mt*N_mt_max + (ks + 1)*kn && i >= N_pol*N_mt*N_mt_max + ks*kn){
                    j = N_pol*N_mt*N_mt_max + (ks + 1)*kn - 1 - chrom_bn;
                    k = N_pol*N_mt*N_mt_max + (ks + 2)*kn - 1 - chrom_bn;
                }
            }
			float3 r1 = r[i];
			float3 r2 = r[j];
			float3 r3 = r[k];

			float3 dr12, dr32;
			dr12.x = r1.x - r2.x;
			dr12.y = r1.y - r2.y;
			dr12.z = r1.z - r2.z;

			dr32.x = r3.x - r2.x;
			dr32.y = r3.y - r2.y;
			dr32.z = r3.z - r2.z;

			float r12inv = 1.0f*powf(dr12.x*dr12.x + dr12.y*dr12.y + dr12.z*dr12.z, -0.5);
			float r32inv = 1.0f*powf(dr32.x*dr32.x + dr32.y*dr32.y + dr32.z*dr32.z, -0.5);
			float costheta = (dr12.x*dr32.x + dr12.y*dr32.y + dr12.z*dr32.z)*r12inv*r32inv;
			if (costheta > 1.0f){
				costheta = 1.0f;
			}
			else	if (costheta < -1.0f){
				costheta = -1.0f;
			}
			float sintheta = sqrtf(1.0f - costheta*costheta);
			float theta = acos(costheta);
			float diff = theta - kt_cos[i - N_pol*N_mt*N_mt_max]; 
            if (sintheta < 1.e-6){
				if (diff < 0){
					diff *= 2.0f*K_kt*100000;
				}
				else{
					diff *= -2.0f*K_kt*100000;
				}
			}
			else{
				diff *= (-2.0f*K_kt*100000)*powf(sintheta, -1);
			}
			float c1 = diff*r12inv;
			float c2 = diff*r32inv;
    
			float3 f1, f2, f3;
			f1.x = c1*(dr12.x*(r12inv*costheta) - dr32.x*r32inv);
			f1.y = c1*(dr12.y*(r12inv*costheta) - dr32.y*r32inv);
			f1.z = c1*(dr12.z*(r12inv*costheta) - dr32.z*r32inv);
			f2 = f1;
			f3.x = c2*(dr32.x*(r32inv*costheta) - dr12.x*r12inv);
			f3.y = c2*(dr32.y*(r32inv*costheta) - dr12.y*r12inv);
			f3.z = c2*(dr32.z*(r32inv*costheta) - dr12.z*r12inv);
			f2.x += f3.x;
			f2.y += f3.y;
			f2.z += f3.z;

			f2.x = -f2.x;
			f2.y = -f2.y;
			f2.z = -f2.z;

            atomicAdd(&f[i].x, f1.x);
			atomicAdd(&f[i].y, f1.y);		
			atomicAdd(&f[i].z, f1.z);
        }
        if (type[i] == SHELL_LEFT){
            int j, k, ks;
            for (ks = 0; ks < N_kt; ks += 2){
                if (i < N_pol*N_mt*N_mt_max + (ks + 1)*kn && i >= N_pol*N_mt*N_mt_max + ks*kn){
                    j = N_pol*N_mt*N_mt_max + (ks + 1)*kn - 1 - chrom_bn;
                    k = N_pol*N_mt*N_mt_max + (ks + 2)*kn - 1 - chrom_bn;
                }
            }
			float3 r1 = r[i];
			float3 r2 = r[j];
			float3 r3 = r[k];

			float3 dr12, dr32;
			dr12.x = r1.x - r2.x;
			dr12.y = r1.y - r2.y;
			dr12.z = r1.z - r2.z;

			dr32.x = r3.x - r2.x;
			dr32.y = r3.y - r2.y;
			dr32.z = r3.z - r2.z;

			float r12inv = 1.0f*powf(dr12.x*dr12.x + dr12.y*dr12.y + dr12.z*dr12.z, -0.5);
			float r32inv = 1.0f*powf(dr32.x*dr32.x + dr32.y*dr32.y + dr32.z*dr32.z, -0.5);
			float costheta = (dr12.x*dr32.x + dr12.y*dr32.y + dr12.z*dr32.z)*r12inv*r32inv;
			if (costheta > 1.0f){
				costheta = 1.0f;
			}
			else	if (costheta < -1.0f){
				costheta = -1.0f;
			}
			float sintheta = sqrtf(1.0f - costheta*costheta);
			float theta = acos(costheta);
			float diff = theta - kt_cos[i - N_pol*N_mt*N_mt_max]; 
            if (sintheta < 1.e-6){
				if (diff < 0){
					diff *= 2.0f*K_kt*1000000;
				}
				else{
					diff *= -2.0f*K_kt*1000000;
				}
			}
			else{
				diff *= (-2.0f*K_kt*1000000)*powf(sintheta, -1);
			}
			float c1 = diff*r12inv;
			float c2 = diff*r32inv;
    
			float3 f1, f2, f3;
			f1.x = c1*(dr12.x*(r12inv*costheta) - dr32.x*r32inv);
			f1.y = c1*(dr12.y*(r12inv*costheta) - dr32.y*r32inv);
			f1.z = c1*(dr12.z*(r12inv*costheta) - dr32.z*r32inv);
			f2 = f1;
			f3.x = c2*(dr32.x*(r32inv*costheta) - dr12.x*r12inv);
			f3.y = c2*(dr32.y*(r32inv*costheta) - dr12.y*r12inv);
			f3.z = c2*(dr32.z*(r32inv*costheta) - dr12.z*r12inv);
			f2.x += f3.x;
			f2.y += f3.y;
			f2.z += f3.z;

			f2.x = -f2.x;
			f2.y = -f2.y;
			f2.z = -f2.z;

            atomicAdd(&f[i].x, f1.x);
			atomicAdd(&f[i].y, f1.y);		
			atomicAdd(&f[i].z, f1.z);
        }
        else    if (type[i] == RIGHT_NDC){
            int j, k, ks;

            for (ks = 1; ks < N_kt; ks += 2){
                if (i < N_pol*N_mt*N_mt_max + (ks + 1)*kn && i >= N_pol*N_mt*N_mt_max + ks*kn){
                    k = N_pol*N_mt*N_mt_max + ks*kn - 1 - chrom_bn;
                    j = N_pol*N_mt*N_mt_max + (ks + 1)*kn - 1 - chrom_bn;
                }
            }

			float3 r1 = r[i];
			float3 r2 = r[j];
			float3 r3 = r[k];
            
			float3 dr12, dr32;
			dr12.x = r1.x - r2.x;
			dr12.y = r1.y - r2.y;
			dr12.z = r1.z - r2.z;

			dr32.x = r3.x - r2.x;
			dr32.y = r3.y - r2.y;
			dr32.z = r3.z - r2.z;

			float r12inv = 1.0f*powf(dr12.x*dr12.x + dr12.y*dr12.y + dr12.z*dr12.z, -0.5);
			float r32inv = 1.0f*powf(dr32.x*dr32.x + dr32.y*dr32.y + dr32.z*dr32.z, -0.5);
			float costheta = (dr12.x*dr32.x + dr12.y*dr32.y + dr12.z*dr32.z)*r12inv*r32inv;
			if (costheta > 1.0f){
				costheta = 1.0f;
			}
			else	if (costheta < -1.0f){
				costheta = -1.0f;
			}
			float sintheta = sqrtf(1.0f - costheta*costheta);
			float theta = acos(costheta);
			float diff = theta - kt_cos[i - N_pol*N_mt*N_mt_max];
			if (sintheta < 1.e-6){
				if (diff < 0){
					diff *= 2.0f*K_kt*100000;
				}
				else{
					diff *= -2.0f*K_kt*100000;
				}
			}
			else{
				diff *= (-2.0f*K_kt*100000)*powf(sintheta, -1);
			}
			float c1 = diff*r12inv;
			float c2 = diff*r32inv;

			float3 f1, f2, f3;
			f1.x = c1*(dr12.x*(r12inv*costheta) - dr32.x*r32inv);
			f1.y = c1*(dr12.y*(r12inv*costheta) - dr32.y*r32inv);
			f1.z = c1*(dr12.z*(r12inv*costheta) - dr32.z*r32inv);
			f2 = f1;
			f3.x = c2*(dr32.x*(r32inv*costheta) - dr12.x*r12inv);
			f3.y = c2*(dr32.y*(r32inv*costheta) - dr12.y*r12inv);
			f3.z = c2*(dr32.z*(r32inv*costheta) - dr12.z*r12inv);
			f2.x += f3.x;
			f2.y += f3.y;
			f2.z += f3.z;

			f2.x = -f2.x;
			f2.y = -f2.y;
			f2.z = -f2.z;

            atomicAdd(&f[i].x, f1.x);
			atomicAdd(&f[i].y, f1.y);		
			atomicAdd(&f[i].z, f1.z);
        }
        else    if (type[i] == SHELL_RIGHT){
            int j, k;
            for (int ks = 1; ks < N_kt; ks += 2){
                if (i < N_pol*N_mt*N_mt_max + (ks + 1)*kn && i >= N_pol*N_mt*N_mt_max + ks*kn){
                    k = N_pol*N_mt*N_mt_max + ks*kn - 1 - chrom_bn;
                    j = N_pol*N_mt*N_mt_max + (ks + 1)*kn - 1 - chrom_bn;
                }
            }
			float3 r1 = r[i];
			float3 r2 = r[j];
			float3 r3 = r[k];

			float3 dr12, dr32;
			dr12.x = r1.x - r2.x;
			dr12.y = r1.y - r2.y;
			dr12.z = r1.z - r2.z;

			dr32.x = r3.x - r2.x;
			dr32.y = r3.y - r2.y;
			dr32.z = r3.z - r2.z;

			float r12inv = 1.0f*powf(dr12.x*dr12.x + dr12.y*dr12.y + dr12.z*dr12.z, -0.5);
			float r32inv = 1.0f*powf(dr32.x*dr32.x + dr32.y*dr32.y + dr32.z*dr32.z, -0.5);
			float costheta = (dr12.x*dr32.x + dr12.y*dr32.y + dr12.z*dr32.z)*r12inv*r32inv;
			if (costheta > 1.0f){
				costheta = 1.0f;
			}
			else	if (costheta < -1.0f){
				costheta = -1.0f;
			}
			float sintheta = sqrtf(1.0f - costheta*costheta);
			float theta = acos(costheta);
			float diff = theta - kt_cos[i - N_pol*N_mt*N_mt_max];
			if (sintheta < 1.e-6){
				if (diff < 0){
					diff *= 2.0f*K_kt*1000000;
				}
				else{
					diff *= -2.0f*K_kt*1000000;
				}
			}
			else{
				diff *= (-2.0f*K_kt*1000000)*powf(sintheta, -1);
			}
			float c1 = diff*r12inv;
			float c2 = diff*r32inv;

			float3 f1, f2, f3;
			f1.x = c1*(dr12.x*(r12inv*costheta) - dr32.x*r32inv);
			f1.y = c1*(dr12.y*(r12inv*costheta) - dr32.y*r32inv);
			f1.z = c1*(dr12.z*(r12inv*costheta) - dr32.z*r32inv);
			f2 = f1;
			f3.x = c2*(dr32.x*(r32inv*costheta) - dr12.x*r12inv);
			f3.y = c2*(dr32.y*(r32inv*costheta) - dr12.y*r12inv);
			f3.z = c2*(dr32.z*(r32inv*costheta) - dr12.z*r12inv);
			f2.x += f3.x;
			f2.y += f3.y;
			f2.z += f3.z;

			f2.x = -f2.x;
			f2.y = -f2.y;
			f2.z = -f2.z;

            atomicAdd(&f[i].x, f1.x);
			atomicAdd(&f[i].y, f1.y);		
			atomicAdd(&f[i].z, f1.z);
        }
    //	__syncthreads();
	}	
}	
