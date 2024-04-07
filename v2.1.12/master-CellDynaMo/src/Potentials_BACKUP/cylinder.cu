#include "cylinder.cuh"		

__global__ void cyl_cyl(float3* r, float3* f, int N, int* type, Param* d_parameters){
    int N_pol = d_parameters->pole_n;
    int N_mt = d_parameters->mt_n;
    int N_mt_max = d_parameters->max_mt_n;
    float e_lj = d_parameters->e_rep/0.23;
    float mt_rad = d_parameters->mt_r;
	int i = blockIdx.x*blockDim.x + threadIdx.x;	
	if (i < N_mt*N_mt_max){
        if (type[i] == MT_REG || type[i] == MT_END){
		    float3 ri1 = r[i];  //a0
		    float3 rj1 = r[i + 1];  //a1
            
            float3 intersect1, intersect2;
            float min_dist;
            float dx, dy, dz;

            float3 p_vector1, p1_norm;
            p_vector1.x = rj1.x - ri1.x;  //A
            p_vector1.y = rj1.y - ri1.y;
            p_vector1.z = rj1.z - ri1.z;

            float dp1 = sqrtf(p_vector1.x*p_vector1.x + p_vector1.y*p_vector1.y + p_vector1.z+ p_vector1.z); //magA
            p1_norm.x = p_vector1.x/dp1; //_A
            p1_norm.y = p_vector1.y/dp1;
            p1_norm.z = p_vector1.z/dp1;

            for (int j = N_mt*N_mt_max; j < N_pol*N_mt*N_mt_max; j++){
                if ((type[j] == MT_REG || type[j] == MT_END)){
		            float3 ri2 = r[j]; //b0
		            float3 rj2 = r[j + 1];  //b1

                    float3 p_vector2, p2_norm;
                    p_vector2.x = rj2.x - ri2.x; //B
                    p_vector2.y = rj2.y - ri2.y;
                    p_vector2.z = rj2.z - ri2.z;

                    float dp2 = sqrtf(p_vector2.x*p_vector2.x + p_vector2.y*p_vector2.y + p_vector2.z+ p_vector2.z); //magB
                    p2_norm.x = p_vector2.x/dp2; //_B
                    p2_norm.y = p_vector2.y/dp2;
                    p2_norm.z = p_vector2.z/dp2;

                    float3 cross;
                    cross.x = p1_norm.y*p2_norm.z - p2_norm.y*p1_norm.z; //cross
                    cross.y = p1_norm.z*p2_norm.x - p2_norm.z*p1_norm.x;
                    cross.z = p1_norm.x*p2_norm.y - p2_norm.x*p1_norm.y;

                    float d_cross = sqrtf(cross.x*cross.x + cross.y*cross.y + cross.z*cross.z); //denom
                    d_cross = d_cross*d_cross;
                    if (d_cross == 0){
                        float d0 = p1_norm.x*(ri2.x - ri1.x) + p1_norm.y*(ri2.y - ri1.y) + p1_norm.z*(ri2.z - ri1.z);
                        float d1 = p1_norm.x*(rj2.x - ri1.x) + p1_norm.y*(rj2.y - ri1.y) + p1_norm.z*(rj2.z - ri1.z);
                        if (d0 <= 0 && 0 >= d1){
                            if (sqrtf(d0*d0) < sqrtf(d1*d1)){
                                intersect1 = ri1;
                                intersect2 = ri2;
                                dx = intersect2.x - intersect1.x;
                                dy = intersect2.y - intersect1.y;
                                dz = intersect2.z - intersect1.z;
                                min_dist = sqrtf(dx*dx + dy*dy + dz*dz);
                            }
                            else    {
                                intersect1 = ri1;
                                intersect2 = rj2;
                                
                                dx = intersect2.x - intersect1.x;
                                dy = intersect2.y - intersect1.y;
                                dz = intersect2.z - intersect1.z;
                                min_dist = sqrtf(dx*dx + dy*dy + dz*dz);
                            }
                        }
                        else if (d0 >= dp1 && dp1 <= d1){
                            if (sqrtf(d0*d0) < sqrtf(d1*d1)){
                                intersect1 = rj1;
                                intersect2 = ri2;
                                dx = intersect2.x - intersect1.x;
                                dy = intersect2.y - intersect1.y;
                                dz = intersect2.z - intersect1.z;
                                min_dist = sqrtf(dx*dx + dy*dy + dz*dz);
                            }
                            else    {
                                intersect1 = rj1;
                                intersect2 = rj2;
                                dx = intersect2.x - intersect1.x;
                                dy = intersect2.y - intersect1.y;
                                dz = intersect2.z - intersect1.z;
                                min_dist = sqrtf(dx*dx + dy*dy + dz*dz);
                            }
                        }
                    }
                    else {
                        float3 t;
                        t.x = ri2.x - ri1.x;
                        t.y = ri2.y - ri1.y;
                        t.z = ri2.z - ri1.z;

                        float det0 = t.x*(p2_norm.y*cross.z - cross.y*p2_norm.z) + t.y*(p2_norm.z*cross.x - cross.z*p2_norm.x) + t.z*(p2_norm.x*cross.y - cross.x*p2_norm.y);
                        float det1 = t.x*(p1_norm.y*cross.z - cross.y*p1_norm.z) + t.y*(p1_norm.z*cross.x - cross.z*p1_norm.x) + t.z*(p1_norm.x*cross.y - cross.x*p1_norm.y);
                        float t0 = det0/d_cross;
                        float t1 = det1/d_cross;

                        float3 p1, p2;
                        p1.x = ri1.x + t0*p1_norm.x;
                        p1.y = ri1.y + t0*p1_norm.y;
                        p1.z = ri1.z + t0*p1_norm.z;

                        p2.x = ri2.x + t1*p2_norm.x;
                        p2.y = ri2.y + t1*p2_norm.y;
                        p2.z = ri2.z + t1*p2_norm.z;

                        if (t0 < 0){
                            p1.x = ri1.x;
                            p1.y = ri1.y;
                            p1.z = ri1.z;
                        }
                        else if (t0 > dp1){
                            p1.x = ri2.x;
                            p1.y = ri2.y;
                            p1.z = ri2.z;
                        }
                        if (t1 < 0){
                            p2.x = rj1.x;
                            p2.y = rj1.y;
                            p2.z = rj1.z;
                        }
                        else if (t1 > dp2){
                            p2.x = rj2.x;
                            p2.y = rj2.y;
                            p2.z = rj2.z;
                        }
                        float dot;
                        if (t0 < 0 || t0 > dp1){
                            dot = p2_norm.x*(p1.x - rj1.x) + p2_norm.y*(p1.y - rj1.y) + p2_norm.z*(p1.z - rj1.z);
                            if (dot < 0){
                                dot = 0.0;
                            }
                            else if (dot > dp2){
                                dot = dp2;
                            }
                            p2.x = rj1.x + dot*p2_norm.x;
                            p2.y = rj1.y + dot*p2_norm.y;
                            p2.z = rj1.z + dot*p2_norm.z;
                        }

                        if (t1 < 0 || t1 > dp2){
                            dot = p1_norm.x*(p2.x - ri1.x) + p1_norm.y*(p2.y - ri1.y) + p1_norm.z*(p2.z - ri1.z);
                            if (dot < 0){
                                dot = 0.0;
                            }
                            else if (dot > dp1){
                                dot = dp1;
                            }
                            p1.x = ri1.x + dot*p1_norm.x;
                            p1.y = ri1.y + dot*p1_norm.y;
                            p1.z = ri1.z + dot*p1_norm.z;
                        }

                        intersect1 = p1;
                        intersect2 = p2;
                        dx = p2.x - p1.x;
                        dy = p2.y - p1.y;
                        dz = p2.z - p1.z;
                        min_dist = sqrtf(dx*dx + dy*dy + dz*dz);
                    }

                    if (min_dist < 2.0*mt_rad){
                        float3 d_a;
                        d_a.x = intersect1.x - ri1.x;
                        d_a.y = intersect1.y - ri1.y;
                        d_a.z = intersect1.z - ri1.z;
                        float alpha = sqrtf(d_a.x*d_a.x + d_a.y* d_a.y + d_a.z*d_a.z)/dp1; 

                        float3 d_b;
                        d_b.x = intersect2.x - ri2.x;
                        d_b.y = intersect2.y - ri2.y;
                        d_b.z = intersect2.z - ri2.z; 
                        float betha = sqrtf(d_b.x*d_b.x + d_b.y* d_b.y + d_b.z*d_b.z)/dp2;

                        float df = -6*e_lj*powf((2.0*mt_rad/min_dist), 12)*powf(min_dist, -1)*powf(min_dist, -1);
                        atomicAdd(&f[i].x, df*(1 - alpha)*dx);
                        atomicAdd(&f[i].y, df*(1 - alpha)*dy);  		
                        atomicAdd(&f[i].z, df*(1 - alpha)*dz);

                        atomicAdd(&f[i + 1].x, df*alpha*dx);
                        atomicAdd(&f[i + 1].y, df*alpha*dy);  		
                        atomicAdd(&f[i + 1].z, df*alpha*dz);

                        atomicAdd(&f[j].x, -df*(1 - betha)*dx);
                        atomicAdd(&f[j].y, -df*(1 - betha)*dy);		
                        atomicAdd(&f[j].z, -df*(1 - betha)*dz); 

                        atomicAdd(&f[j + 1].x, -df*betha*dx);
                        atomicAdd(&f[j + 1].y, -df*betha*dy);  		
                        atomicAdd(&f[j + 1].z, -df*betha*dz);
                    } 
                }
            }
        }
    }
}

__global__ void cyl_sphere(float3* r, float3* f, int N, int* type, int kn, Param* d_parameters){
    int N_pol = d_parameters->pole_n;
    int N_mt = d_parameters->mt_n;
    int N_kt = d_parameters->kin_n;
    int N_mt_max = d_parameters->max_mt_n;    
    float e_lj = d_parameters->e_rep/0.23;
    float mt_rad = d_parameters->mt_r;
    float kt_rad = d_parameters->kt_r;
    int chrom_cond = d_parameters->chrom;
    int x_num = d_parameters->x_num;
    int y_num = d_parameters->y_num;
    int chrom_bn = chrom_cond*d_parameters->chrom_num;
	int i = blockIdx.x*blockDim.x + threadIdx.x;	
	if (i < 2*N_pol*N_mt*N_mt_max){
        if (type[i] == MT_REG || type[i] == MT_END || type[i] == LEFT_MINUS || type[i] == RIGHT_MINUS){
            printf("CYL: %d\t%d\n", i, N_pol);
		    float3 ri1 = r[i];  //a0
		    float3 rj1 = r[i + 1];  //a1
            
            float3 intersect;
            float min_dist;
            float dx, dy, dz;

            float3 p_vector1, p1_norm;
            p_vector1.x = rj1.x - ri1.x;  //A
            p_vector1.y = rj1.y - ri1.y;
            p_vector1.z = rj1.z - ri1.z;

            float dp1 = sqrtf(p_vector1.x*p_vector1.x + p_vector1.y*p_vector1.y + p_vector1.z+ p_vector1.z); //magA
            p1_norm.x = p_vector1.x/dp1; //_A
            p1_norm.y = p_vector1.y/dp1;
            p1_norm.z = p_vector1.z/dp1;
            
            for (int h = 0; h < N_kt; h++){
		        int s = N_pol*N_mt*N_mt_max + (h + 1)*kn - 1 - chrom_bn;
                float3 rs = r[s];

                float3 p_vectorS, pS_Scale;
                p_vectorS.x = rs.x - ri1.x;  //A
                p_vectorS.y = rs.y - ri1.y;
                p_vectorS.z = rs.z - ri1.z;

                pS_Scale.x = p_vectorS.x/dp1;
                pS_Scale.y = p_vectorS.y/dp1;
                pS_Scale.z = p_vectorS.z/dp1;
                float dpS = sqrtf(p_vectorS.x*p_vectorS.x + p_vectorS.y*p_vectorS.y + p_vectorS.z+ p_vectorS.z); //magA

                float t = p1_norm.x*pS_Scale.x + p1_norm.y*pS_Scale.y + p1_norm.z*pS_Scale.z;
                if (t < 0.0){
                    t = 0.0;
                }
                if (t > 1.0){
                    t = 1.0;
                }
                intersect.x = p_vector1.x*t;
                intersect.y = p_vector1.y*t;
                intersect.z = p_vector1.z*t;
                dx = intersect.x - p_vectorS.x;
                dy = intersect.y - p_vectorS.y;
                dz = intersect.z - p_vectorS.z;
                min_dist = sqrtf(dx*dx + dy*dy + dz*dz);
               
                if (min_dist < mt_rad + kt_rad){
                    float3 d_a;
                    d_a.x = intersect.x - ri1.x;
                    d_a.y = intersect.y - ri1.y;
                    d_a.z = intersect.z - ri1.z;
                    float alpha = sqrtf(d_a.x*d_a.x + d_a.y* d_a.y + d_a.z*d_a.z)/dp1; 

                   // float df = -6*e_lj*powf(((mt_rad + kt_rad)/min_dist), 12)*powf(min_dist, -1)*powf(min_dist, -1);
                    float df = -6*10.0*powf(((mt_rad + kt_rad)/min_dist), 12)*powf(min_dist, -1)*powf(min_dist, -1)/10000;
                    atomicAdd(&f[i].x, -df*(1 - alpha)*dx);
                    atomicAdd(&f[i].y, -df*(1 - alpha)*dy);  		
                    atomicAdd(&f[i].z, -df*(1 - alpha)*dz);

                    atomicAdd(&f[i + 1].x, -df*alpha*dx);
                    atomicAdd(&f[i + 1].y, -df*alpha*dy);  		
                    atomicAdd(&f[i + 1].z, -df*alpha*dz);

                    atomicAdd(&f[s].x, df*dx);
                    atomicAdd(&f[s].y, df*dy);		
                    atomicAdd(&f[s].z, df*dz); 
                }
            }
            printf("check: %d\t%d", chrom_bn, N_kt);
            for (int p = 0; p < N_kt; p++){
                for (int ch = 1; ch <= chrom_bn; ch++){
                    int s = N_pol*N_mt*N_mt_max + (p + 1)*kn - ch;
                    printf("SAS:%d\n", s);
                    float3 rs = r[s];

                    float3 p_vectorS, pS_Scale;
                    p_vectorS.x = rs.x - ri1.x;  //A
                    p_vectorS.y = rs.y - ri1.y;
                    p_vectorS.z = rs.z - ri1.z;

                    pS_Scale.x = p_vectorS.x/dp1;
                    pS_Scale.y = p_vectorS.y/dp1;
                    pS_Scale.z = p_vectorS.z/dp1;
                    float dpS = sqrtf(p_vectorS.x*p_vectorS.x + p_vectorS.y*p_vectorS.y + p_vectorS.z+ p_vectorS.z); //magA

                    float t = p1_norm.x*pS_Scale.x + p1_norm.y*pS_Scale.y + p1_norm.z*pS_Scale.z;
                    if (t < 0.0){
                        t = 0.0;
                    }
                    if (t > 1.0){
                        t = 1.0;
                    }
                    intersect.x = p_vector1.x*t;
                    intersect.y = p_vector1.y*t;
                    intersect.z = p_vector1.z*t;
                    dx = intersect.x - p_vectorS.x;
                    dy = intersect.y - p_vectorS.y;
                    dz = intersect.z - p_vectorS.z;
                    min_dist = sqrtf(dx*dx + dy*dy + dz*dz);
                   
                    if (min_dist < mt_rad + kt_rad){
                        printf("CHROM\n");
                        float3 d_a;
                        d_a.x = intersect.x - ri1.x;
                        d_a.y = intersect.y - ri1.y;
                        d_a.z = intersect.z - ri1.z;
                        float alpha = sqrtf(d_a.x*d_a.x + d_a.y* d_a.y + d_a.z*d_a.z)/dp1; 

                        //float df = -6*e_lj*powf(((mt_rad + kt_rad)/min_dist), 12)*powf(min_dist, -1)*powf(min_dist, -1);
                        float df = -6*10.0*powf(((mt_rad + kt_rad)/min_dist), 12)*powf(min_dist, -1)*powf(min_dist, -1)/10000;
                        atomicAdd(&f[i].x, -df*(1 - alpha)*dx);
                        atomicAdd(&f[i].y, -df*(1 - alpha)*dy);  		
                        atomicAdd(&f[i].z, -df*(1 - alpha)*dz);

                        atomicAdd(&f[i + 1].x, -df*alpha*dx);
                        atomicAdd(&f[i + 1].y, -df*alpha*dy);  		
                        atomicAdd(&f[i + 1].z, -df*alpha*dz);

                        atomicAdd(&f[s].x, df*dx);
                        atomicAdd(&f[s].y, df*dy);		
                        atomicAdd(&f[s].z, df*dz); 
                    }
                }
                for (int sh = x_num*y_num; sh < kn; sh++){
                    int s = N_pol*N_mt*N_mt_max + p*kn + sh;
                    //printf("SHELL: %d\n", s);
                    
                    float3 rs = r[s];

                    float3 p_vectorS, pS_Scale;
                    p_vectorS.x = rs.x - ri1.x;  //A
                    p_vectorS.y = rs.y - ri1.y;
                    p_vectorS.z = rs.z - ri1.z;

                    pS_Scale.x = p_vectorS.x/dp1;
                    pS_Scale.y = p_vectorS.y/dp1;
                    pS_Scale.z = p_vectorS.z/dp1;
                    float dpS = sqrtf(p_vectorS.x*p_vectorS.x + p_vectorS.y*p_vectorS.y + p_vectorS.z+ p_vectorS.z); //magA

                    float t = p1_norm.x*pS_Scale.x + p1_norm.y*pS_Scale.y + p1_norm.z*pS_Scale.z;
                    if (t < 0.0){
                        t = 0.0;
                    }
                    if (t > 1.0){
                        t = 1.0;
                    }
                    intersect.x = p_vector1.x*t;
                    intersect.y = p_vector1.y*t;
                    intersect.z = p_vector1.z*t;
                    dx = intersect.x - p_vectorS.x;
                    dy = intersect.y - p_vectorS.y;
                    dz = intersect.z - p_vectorS.z;
                    min_dist = sqrtf(dx*dx + dy*dy + dz*dz);
                   
                    if (min_dist < mt_rad + kt_rad){
                        float3 d_a;
                        d_a.x = intersect.x - ri1.x;
                        d_a.y = intersect.y - ri1.y;
                        d_a.z = intersect.z - ri1.z;
                        float alpha = sqrtf(d_a.x*d_a.x + d_a.y* d_a.y + d_a.z*d_a.z)/dp1; 

                        //float df = -6*e_lj*powf(((mt_rad + kt_rad)/min_dist), 12)*powf(min_dist, -1)*powf(min_dist, -1);
                        float df = -6*10.0*powf(((mt_rad + kt_rad)/min_dist), 12)*powf(min_dist, -1)*powf(min_dist, -1)/10000;
                        atomicAdd(&f[i].x, -df*(1 - alpha)*dx);
                        atomicAdd(&f[i].y, -df*(1 - alpha)*dy);  		
                        atomicAdd(&f[i].z, -df*(1 - alpha)*dz);

                        atomicAdd(&f[i + 1].x, -df*alpha*dx);
                        atomicAdd(&f[i + 1].y, -df*alpha*dy);  		
                        atomicAdd(&f[i + 1].z, -df*alpha*dz);

                        /*atomicAdd(&f[s].x, df*dx);
                        atomicAdd(&f[s].y, df*dy);		
                        atomicAdd(&f[s].z, df*dz); */
                    }
                }
            }      
        }
    }                    
}

/*__global__ void end_force(float3* r, float3* f, int N){

	int i = blockIdx.x*blockDim.x + threadIdx.x;	
	if (i == 0){        
		atomicAdd(&f[i].y, 50.0);		
	}
	if (i == 1){        
		atomicAdd(&f[i].y, 50.0);		
	}
}*/
