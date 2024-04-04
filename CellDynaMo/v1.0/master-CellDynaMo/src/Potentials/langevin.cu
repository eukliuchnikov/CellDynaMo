#include "langevin.cuh"
#include "../Chemistry/react.h"
#include "../Math/ht.cu"

__global__ void integrateGPU(float3* r, float3* f, int N, int* type, int kn, float3 l_f, float3 r_f, Param* d_parameters){
	int i;
    float time_step = d_parameters->timestep;
    float viscos = d_parameters->visc;
    float temperatur = d_parameters->temper;
	i = blockIdx.x*blockDim.x + threadIdx.x;
	if (i < N){
        //float4 rf = rforce(i);
        //float pole_gamma = 6*M_PI*viscos*d_parameters->pole_r;
        //float pole_var = sqrtf(KB*temperatur*2.0*time_step*pow(pole_gamma, -1));
        float kt_gamma = 6*M_PI*viscos*d_parameters->ndc_d/2.0;
        //float kt_var = sqrtf(KB*temperatur*2.0*time_step*pow(kt_gamma, -1));
        float ktc_gamma = 6*M_PI*viscos*d_parameters->kt_r;
        float shell_gamma = 6*M_PI*viscos*60.0;
        float ktc_var = sqrtf(KB*temperatur*2.0*time_step*pow(ktc_gamma, -1));
//MT dynamics
		if (type[i] == MT_REG || type[i] == PLUS_DET || type[i] == PLUS_ATT || type[i] == PLUS_DET_INVALID){
            float gamma = 6*M_PI*viscos*d_parameters->mt_r;
            float var = sqrtf(KB*temperatur*2.0*time_step*pow(gamma, -1));
			r[i].x += f[i].x*time_step*pow(gamma, -1)/* + var*rf.x*/;
			r[i].y += f[i].y*time_step*pow(gamma, -1);
			r[i].z += f[i].z*time_step*pow(gamma, -1);        
		}
//Kinetochores
		if (type[i] == LEFT_NDC || type[i] == RIGHT_NDC){
			r[i].x += f[i].x*time_step*pow(kt_gamma, -1);
			r[i].y += f[i].y*time_step*pow(kt_gamma, -1);
			r[i].z += f[i].z*time_step*pow(kt_gamma, -1);            
		}

		else	if (type[i] == LEFT_KT || type[i] == RIGHT_KT || type[i] == CHROM){
			r[i].x += f[i].x*time_step*pow(ktc_gamma, -1);/* + ktc_var*rf.x*///;
			r[i].y += f[i].y*time_step*pow(ktc_gamma, -1);
			r[i].z += f[i].z*time_step*pow(ktc_gamma, -1);
		}
        else	if (type[i] == SHELL_LEFT || type[i] == SHELL_RIGHT){
			r[i].x += f[i].x*time_step*pow(shell_gamma, -1);/* + ktc_var*rf.x*///;
			r[i].y += f[i].y*time_step*pow(shell_gamma, -1);
			r[i].z += f[i].z*time_step*pow(shell_gamma, -1);
		}
		/*else	if (type[i] == LEFT_POLE || type[i] == LEFT_MINUS){
			r[i].x += l_f.x;
			r[i].y += l_f.y;
			r[i].z += l_f.z;
        }
        else	if (type[i] == RIGHT_POLE || type[i] == RIGHT_MINUS){
			r[i].x += r_f.x;
			r[i].y += r_f.y;
			r[i].z += r_f.z;
        }*/
        f[i].x = 0.0f;
    	f[i].y = 0.0f;
	    f[i].z = 0.0f;
	   // __syncthreads();
	}
}

__global__ void CoM(float3* r, float3* f, int N, int* type, int kn){
	int i;
	i = blockIdx.x*blockDim.x + threadIdx.x;
	if (i < N){
//MT dynamics
//Kinetochores
		if (type[i] == LEFT_KT){
			float dx = (r[i].x + r[i + kn].x)/2.0;
			float dy = (r[i].y + r[i + kn].y)/2.0;
			float dz = (r[i].z + r[i + kn].z)/2.0;
            r[i].x -= dx;
            r[i].y -= dy;
            r[i].z -= dz;
            r[i + kn].x -= dx;
            r[i + kn].y -= dy;
            r[i + kn].z -= dz;            
		}
	   // __syncthreads();
	}
}

__global__ void integrateGPU_init(float3* r, float3* f, int N, int* type, int kn, Param* d_parameters){
	int i;
    float time_step = d_parameters->timestep;
    float viscos = d_parameters->visc;
    float temperatur = d_parameters->temper;
	i = blockIdx.x*blockDim.x + threadIdx.x;
    //float mult = 
	if (i < N){
        float kt_gamma = 6*M_PI*viscos*d_parameters->ndc_d/2.0;
        //float kt_var = sqrtf(KB*temperatur*2.0*time_step*pow(kt_gamma, -1));
        float ktc_gamma = 6*M_PI*viscos*d_parameters->kt_r;
        //float ktc_var = sqrtf(KB*temperatur*2.0*time_step*pow(ktc_gamma, -1));
		if (type[i] == LEFT_NDC || type[i] == RIGHT_NDC){
            
			r[i].x += f[i].x*time_step*10*pow(kt_gamma, -1);
			r[i].y += f[i].y*time_step*10*pow(kt_gamma, -1);
			r[i].z += f[i].z*time_step*10*pow(kt_gamma, -1);            
		}

		else	if (type[i] == LEFT_KT || type[i] == RIGHT_KT || type[i] == CHROM){
            if (f[i].x > 1000000000.0){
                f[i].x = 1000000000.0;
            }
            if (f[i].y > 1000000000.0){
                f[i].y = 1000000000.0;
            }
            if (f[i].z > 1000000000.0){
                f[i].z = 1000000000.0;
            }
            if (f[i].x < -1000000000.0){
                f[i].x = -1000000000.0;
            }
            if (f[i].y < -1000000000.0){
                f[i].y = -1000000000.0;
            }
            if (f[i].z < -1000000000.0){
                f[i].z = -1000000000.0;
            }
			r[i].x += f[i].x*time_step*10*pow(ktc_gamma, -1);/* + ktc_var*rf.x*///;
			r[i].y += f[i].y*time_step*10*pow(ktc_gamma, -1);
            r[i].z += f[i].z*time_step*10*pow(ktc_gamma, -1);
		}
        f[i].x = 0.0f;
    	f[i].y = 0.0f;
	    f[i].z = 0.0f;
	   // __syncthreads();
	}
}

float3 kin_force(int index, float3* f, Param* d_parameters){
    float3 result;
    float time_step = d_parameters->timestep;
    float viscos = d_parameters->visc;
    float temperatur = d_parameters->temper;
    float ktc_gamma = 6*M_PI*viscos*d_parameters->kt_r;
    float ktc_var = sqrtf(KB*temperatur*2.0*time_step*pow(ktc_gamma, -1));
    result.x = f[index].x*time_step*pow(ktc_gamma, -1);
    result.y = f[index].y*time_step*pow(ktc_gamma, -1);
    result.z = f[index].z*time_step*pow(ktc_gamma, -1)/* - ktc_var*gasdev(&seed)*/;
    f[index].x = 0.0f;
    f[index].y = 0.0f;
    f[index].z = 0.0f;
    return result;
}

float3 pole_force(int index, float3* f, Param* d_parameters){
    float3 result;
    float time_step = d_parameters->timestep;
    float viscos = d_parameters->visc;
    float temperatur = d_parameters->temper;
    float pole_gamma = 6*M_PI*viscos*d_parameters->pole_r;
    float pole_var = sqrtf(KB*temperatur*2.0*time_step*pow(pole_gamma, -1));
    result.x = f[index].x*time_step*pow(pole_gamma, -1)/* - pole_var*gasdev(&seed)*/;
    result.y = f[index].y*time_step*pow(pole_gamma, -1);
    result.z = f[index].z*time_step*pow(pole_gamma, -1);
    f[index].x = 0.0f;
    f[index].y = 0.0f;
    f[index].z = 0.0f;
    return result;
}

void checkCUDAError(){
	cudaError_t error = cudaGetLastError();
	if(error != cudaSuccess){
		printf("CUDA error: %s \n", cudaGetErrorString(error));
exit(0);
	}
}

void KT_bonds_init(float* d_kin_cos, float* kin_cos, int* d_harmonicKinCount, int* a_harmonicKinCount, float* d_harmonicKinRadii, float* a_harmonicKinRadii, int* d_harmonicKin, int* a_harmonicKin, Param* d_parameters){
    int N_kt = d_parameters->kin_n;
    cudaMalloc((void**)&d_kin_cos, N_kt*kn*sizeof(float));
    cudaMemcpy(d_kin_cos, kin_cos, N_kt*kn*sizeof(float), cudaMemcpyHostToDevice);

    cudaMalloc((void**)&d_harmonicKinCount, N_kt*kn*sizeof(int));
    cudaMemcpy(d_harmonicKinCount, a_harmonicKinCount, N_kt*kn*sizeof(int), cudaMemcpyHostToDevice);
    cudaMalloc((void**)&d_harmonicKinRadii, N_kt*kn*maxHarmonicKinPerMonomer*sizeof(float));
    cudaMemcpy(d_harmonicKinRadii, a_harmonicKinRadii, N_kt*kn*maxHarmonicKinPerMonomer*sizeof(float), cudaMemcpyHostToDevice);
    cudaMalloc((void**)&d_harmonicKin, N_kt*kn*maxHarmonicKinPerMonomer*sizeof(int));
    cudaMemcpy(d_harmonicKin, a_harmonicKin, N_kt*kn*maxHarmonicKinPerMonomer*sizeof(int), cudaMemcpyHostToDevice);
}

void cudaRANDinit(){
   initRand(seed, mds.N);
}

void dynamics(DCD dcd, float time){
	long long int step;
    int blockSize = 512;
	int blockNum = mds.N/blockSize + 1;
    cudaMemcpy(mds.d_type, mds.h_type, mds.N*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(mds.d_connector, mds.h_connector, mds.N*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(mds.d_att_l, mds.h_att_l, mds.N*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(mds.d_link_l, mds.h_link_l, mds.N*sizeof(float), cudaMemcpyHostToDevice);
    for (step = 0; step < 20; step++){
        //printf("%d\n", step);
	    //DCD_step(dcd);
        membrane(mds.h_f, mds.h_length);
		cudaMemcpy(mds.d_r, mds.h_r, mds.N*sizeof(float3), cudaMemcpyHostToDevice);	
        cudaMemcpy(mds.d_f, mds.h_f, mds.N*sizeof(float3), cudaMemcpyHostToDevice);
        cudaMemcpy(mds.d_length, mds.h_length, mds.N*sizeof(float), cudaMemcpyHostToDevice);
        //printf("Harmonic\n");       
        computeHarmonic<<<blockNum, blockSize>>>(mds.d_r, mds.d_f, mds.N, mds.d_type, kn, rk, d_harmonicKinCount, d_harmonicKin, d_harmonicKinRadii, d_kin_rad, d_parameters, mds.d_length);
        //printf("Bending\n");
        computeAngles<<<blockNum, blockSize>>>(mds.d_r, mds.d_f, mds.N, mds.d_type, kn, d_kin_cos, d_dyn_cos, d_parameters, mds.d_length);
        cudaMemcpy(mds.d_att_n, mds.h_att_n, mds.N*sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(mds.d_att_f, mds.h_att_f, mds.N*sizeof(float), cudaMemcpyHostToDevice);   
        //printf("Pushing\n");      
        pushing<<<blockNum, blockSize>>>(mds.d_r, mds.d_f, mds.N, mds.d_type, kn, step, mds.d_att_n, mds.d_att_f, d_parameters);
        //printf("Pulling\n");
        pulling<<<blockNum, blockSize>>>(mds.d_r, mds.d_f, mds.N, mds.d_type, kn, step, mds.d_att_n, mds.d_att_f, d_parameters, mds.d_connector, mds.d_length, mds.d_att_l, mds.d_link_l);
        excl_vol_kt<<<blockNum, blockSize>>>(mds.d_r, mds.d_f, mds.N, mds.d_type, kn, d_parameters);
        cudaMemcpy(mds.h_att_n, mds.d_att_n, mds.N*sizeof(int), cudaMemcpyDeviceToHost);
        cudaMemcpy(mds.h_att_f, mds.d_att_f, mds.N*sizeof(float), cudaMemcpyDeviceToHost);     
        for (int ks = 0; ks < h_parameters->kin_n; ks ++){
            if (step == 0 && (mds.h_att_n[ks*5] != 0)){
                printf("KT#%d FORCES:\tN_total\tN_pull\tN_push\tF_pull\tF_push\n%f\t%d\t%d\t%d\t%f\t%f\n", ks, time, mds.h_att_n[ks*5], mds.h_att_n[ks*5 + 1], mds.h_att_n[ks*5 + 2], mds.h_att_f[ks*5], mds.h_att_f[ks*5 + 1]);
            }
            mds.h_att_n[ks*5] = 0;
            mds.h_att_n[ks*5 + 1] = 0;
            mds.h_att_f[ks*5] = 0.0;
            mds.h_att_f[ks*5 + 1] = 0.0;
            mds.h_att_n[ks*5 + 2] = 0;
            /*mds.h_att_n[ks*5 + 3] = 0;
            mds.h_att_n[ks*5 + 4] = 0;
            mds.h_att_n[ks*5 + 5] = 0;*/
        } 
        int ln = mds.N - 2;
        int rn = mds.N - 1;  

        float3 left_pole = pole_force(ln, mds.h_f, h_parameters);
        float3 right_pole = pole_force(rn, mds.h_f, h_parameters);
        //printf("Integr\n");
		integrateGPU<<<blockNum, blockSize>>>(mds.d_r, mds.d_f, mds.N, mds.d_type, kn, left_pole, right_pole, d_parameters);
        //CoM<<<blockNum, blockSize>>>(mds.d_r, mds.d_f, mds.N, mds.d_type, kn);
        cudaMemcpy(mds.h_f, mds.d_f, mds.N*sizeof(float3), cudaMemcpyDeviceToHost);
        cudaMemcpy(mds.h_r, mds.d_r, mds.N*sizeof(float3), cudaMemcpyDeviceToHost);
        cudaMemcpy(mds.h_length, mds.d_length, mds.N*sizeof(float), cudaMemcpyDeviceToHost);
        cudaMemcpy(dyn_cos, d_dyn_cos, KIN_COUNT*kn*sizeof(int), cudaMemcpyDeviceToHost);
        checkCUDAError();
	}
}

void minimiz(DCD dcd){
	long long int step;
    int blockSize = 512;
	int blockNum = mds.N/blockSize + 1;
    cudaMemcpy(mds.d_type, mds.h_type, mds.N*sizeof(int), cudaMemcpyHostToDevice);
    for (step = 0; step < 30000; step++){
        if (step % 10000 == 0){
	        DCD_step(dcd);
        }
        membrane(mds.h_f, mds.h_length);
		cudaMemcpy(mds.d_r, mds.h_r, mds.N*sizeof(float3), cudaMemcpyHostToDevice);	
        cudaMemcpy(mds.d_f, mds.h_f, mds.N*sizeof(float3), cudaMemcpyHostToDevice);
        computeHarmonic<<<blockNum, blockSize>>>(mds.d_r, mds.d_f, mds.N, mds.d_type, kn, rk, d_harmonicKinCount, d_harmonicKin, d_harmonicKinRadii, d_kin_rad, d_parameters, mds.d_length);
        computeAngles<<<blockNum, blockSize>>>(mds.d_r, mds.d_f, mds.N, mds.d_type, kn, d_kin_cos, d_dyn_cos, d_parameters, mds.d_length);
        cudaMemcpy(mds.d_att_n, mds.h_att_n, mds.N*sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(mds.d_att_f, mds.h_att_f, mds.N*sizeof(float), cudaMemcpyHostToDevice);   
        pushing<<<blockNum, blockSize>>>(mds.d_r, mds.d_f, mds.N, mds.d_type, kn, step, mds.d_att_n, mds.d_att_f, d_parameters);
        excl_vol_kt<<<blockNum, blockSize>>>(mds.d_r, mds.d_f, mds.N, mds.d_type, kn, d_parameters);
        float radi = CH_TER;
        chrom_ter<<<blockNum, blockSize>>>(mds.d_r, mds.d_f, mds.N, mds.d_type, kn, d_parameters, radi);
		integrateGPU_init<<<blockNum, blockSize>>>(mds.d_r, mds.d_f, mds.N, mds.d_type, kn, d_parameters);
        cudaMemcpy(mds.h_f, mds.d_f, mds.N*sizeof(float3), cudaMemcpyDeviceToHost);
        cudaMemcpy(mds.h_r, mds.d_r, mds.N*sizeof(float3), cudaMemcpyDeviceToHost);
        cudaMemcpy(mds.h_length, mds.d_length, mds.N*sizeof(float), cudaMemcpyDeviceToHost);
        cudaMemcpy(dyn_cos, d_dyn_cos, KIN_COUNT*kn*sizeof(int), cudaMemcpyDeviceToHost);
        checkCUDAError();
	}
    for (step = 0; step < 30000; step++){
        if (step % 10000 == 0){
	        DCD_step(dcd);
        }
        membrane(mds.h_f, mds.h_length);
		cudaMemcpy(mds.d_r, mds.h_r, mds.N*sizeof(float3), cudaMemcpyHostToDevice);	
        cudaMemcpy(mds.d_f, mds.h_f, mds.N*sizeof(float3), cudaMemcpyHostToDevice);
        computeHarmonic<<<blockNum, blockSize>>>(mds.d_r, mds.d_f, mds.N, mds.d_type, kn, rk, d_harmonicKinCount, d_harmonicKin, d_harmonicKinRadii, d_kin_rad, d_parameters, mds.d_length);
        computeAngles<<<blockNum, blockSize>>>(mds.d_r, mds.d_f, mds.N, mds.d_type, kn, d_kin_cos, d_dyn_cos, d_parameters, mds.d_length);
        cudaMemcpy(mds.d_att_n, mds.h_att_n, mds.N*sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(mds.d_att_f, mds.h_att_f, mds.N*sizeof(float), cudaMemcpyHostToDevice);   
        pushing<<<blockNum, blockSize>>>(mds.d_r, mds.d_f, mds.N, mds.d_type, kn, step, mds.d_att_n, mds.d_att_f, d_parameters);
        excl_vol_kt<<<blockNum, blockSize>>>(mds.d_r, mds.d_f, mds.N, mds.d_type, kn, d_parameters);
        float radi = 0.625*CH_TER;
        chrom_ter<<<blockNum, blockSize>>>(mds.d_r, mds.d_f, mds.N, mds.d_type, kn, d_parameters, radi);
		integrateGPU_init<<<blockNum, blockSize>>>(mds.d_r, mds.d_f, mds.N, mds.d_type, kn, d_parameters);
        cudaMemcpy(mds.h_f, mds.d_f, mds.N*sizeof(float3), cudaMemcpyDeviceToHost);
        cudaMemcpy(mds.h_r, mds.d_r, mds.N*sizeof(float3), cudaMemcpyDeviceToHost);
        cudaMemcpy(mds.h_length, mds.d_length, mds.N*sizeof(float), cudaMemcpyDeviceToHost);
        cudaMemcpy(dyn_cos, d_dyn_cos, KIN_COUNT*kn*sizeof(int), cudaMemcpyDeviceToHost);
        checkCUDAError();
	}
    for (step = 0; step < 30000; step++){
        if (step % 10000 == 0){
	        DCD_step(dcd);
        }
        membrane(mds.h_f, mds.h_length);
		cudaMemcpy(mds.d_r, mds.h_r, mds.N*sizeof(float3), cudaMemcpyHostToDevice);	
        cudaMemcpy(mds.d_f, mds.h_f, mds.N*sizeof(float3), cudaMemcpyHostToDevice);
        computeHarmonic<<<blockNum, blockSize>>>(mds.d_r, mds.d_f, mds.N, mds.d_type, kn, rk, d_harmonicKinCount, d_harmonicKin, d_harmonicKinRadii, d_kin_rad, d_parameters, mds.d_length);
        computeAngles<<<blockNum, blockSize>>>(mds.d_r, mds.d_f, mds.N, mds.d_type, kn, d_kin_cos, d_dyn_cos, d_parameters, mds.d_length);
        cudaMemcpy(mds.d_att_n, mds.h_att_n, mds.N*sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(mds.d_att_f, mds.h_att_f, mds.N*sizeof(float), cudaMemcpyHostToDevice);   
        pushing<<<blockNum, blockSize>>>(mds.d_r, mds.d_f, mds.N, mds.d_type, kn, step, mds.d_att_n, mds.d_att_f, d_parameters);
        excl_vol_kt<<<blockNum, blockSize>>>(mds.d_r, mds.d_f, mds.N, mds.d_type, kn, d_parameters);
        float radi = CH_TER;
        chrom_ter<<<blockNum, blockSize>>>(mds.d_r, mds.d_f, mds.N, mds.d_type, kn, d_parameters, radi);
		integrateGPU_init<<<blockNum, blockSize>>>(mds.d_r, mds.d_f, mds.N, mds.d_type, kn, d_parameters);
        cudaMemcpy(mds.h_f, mds.d_f, mds.N*sizeof(float3), cudaMemcpyDeviceToHost);
        cudaMemcpy(mds.h_r, mds.d_r, mds.N*sizeof(float3), cudaMemcpyDeviceToHost);
        cudaMemcpy(mds.h_length, mds.d_length, mds.N*sizeof(float), cudaMemcpyDeviceToHost);
        cudaMemcpy(dyn_cos, d_dyn_cos, KIN_COUNT*kn*sizeof(int), cudaMemcpyDeviceToHost);
        checkCUDAError();
	}
}
