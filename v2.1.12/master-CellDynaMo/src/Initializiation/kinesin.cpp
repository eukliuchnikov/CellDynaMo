#include "kinesin.h"

void kinesin_init(float3* r, int* type){
	int i, h, k;
	for (h = 0; h < POLE_COUNT; h++){	
		for (k = 0; k < MT_NUM; k++){
			i = POLE_COUNT*MT_NUM*MAX_MT_LENGTH + KIN_COUNT*kn + 2*h*MT_NUM + k*2;
			r[i].x = POLE_DISTANCE;
			r[i].y = 0.0;
			r[i].z = 0.0;
			r[i + 1].x = POLE_DISTANCE;
			r[i + 1].y = 0.0;
			r[i + 1].z = 0.0;

			type[i] = LIG_MT_INACTIVE;
			type[i + 1] = LIG_CH_INACTIVE;
		}
	}
}