#include "phos-dephos.h"

void PhosDephos(int* reactant, int* product, int reactType, int prodType, int l){
	reactant[l]--;
	product[l]++;
	int rseed_kin = rseed + 123;
	int h, m;
	int countNDC = 0;
	int x_num = CORONA_X;
	int y_num = floor(x_num/RATIO);
	int n_chrom = 2*floor((CONTOUR_LENGTH/2 - KIN_RADIUS)/(2*KIN_RADIUS) + 0.5);
	for (h = 0; h < KIN_COUNT; h++){	
		for (m = 0; m < x_num*y_num; m++){
			if (kt[h].cell[m] == l && kt[h].type[m] == reactType){
				countNDC++;
				kt[h].rand[m] = countNDC;
			}
		}
	}
	float c = countNDC*ran2(&rseed_kin);
	for (h = 0; h < KIN_COUNT; h++){	
		for (m = 0; m < x_num*y_num; m++){
			if (kt[h].cell[m] == l && kt[h].type[m] == reactType){
				if (kt[h].rand[m] - 1 < c && c <= kt[h].rand[m]){
					kt[h].type[m] = prodType;
				}
			}
		}
	}	
}