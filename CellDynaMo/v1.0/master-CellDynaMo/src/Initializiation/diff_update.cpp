/*
 *implicit diffusion initiation module
 *
 *  Created on: Mar 15, 2018
 *	Edited on:	Mar 30, 2019
 *      Author: kliuchnikov
 */

#include "diff_update.h"

void diffusionBar(int* cell_type, int* neighbors, int* neighbornum, int* neighbornum_X, int* neighbornum_Y, int* neighbornum_Z, float3* SVcenter){
	int i, j, k;
	for (i = 1; i < num_x_cell + 1; i++){
		for (j = 1; j < num_y_cell + 1; j++){
			for (k = 1; k < num_z_cell + 1; k++){
				int counter = 0;
				int l = COORintoARRAY(i, j, k);
				if (cell_type[l] == DIFFUSION_SPACE  || cell_type[l] == DIFFUSION_MIX_SPACE){
					//X-direction
					if (i != 1 && i != num_x_cell){
						int right = COORintoARRAY(i + 1, j, k);
						int left = COORintoARRAY(i - 1, j, k);
						if ((cell_type[right] == DIFFUSION_SPACE  || cell_type[right] == DIFFUSION_MIX_SPACE) && (cell_type[left] == DIFFUSION_SPACE  || cell_type[left] == DIFFUSION_MIX_SPACE)){
							neighbors[l*6 + counter] = right;
							counter++;
							neighbors[l*6 + counter] = left;
							counter++; 
							neighbornum_X[l] = 2;
						}
						else	if ((cell_type[right] == DIFFUSION_SPACE  || cell_type[right] == DIFFUSION_MIX_SPACE) && (cell_type[left] != DIFFUSION_SPACE  && cell_type[left] != DIFFUSION_MIX_SPACE)){
							neighbors[l*6 + counter] = right;
							counter++;
							neighbornum_X[l] = 1;					
						}
						else	if ((cell_type[right] != DIFFUSION_SPACE  && cell_type[right] != DIFFUSION_MIX_SPACE) && (cell_type[left] == DIFFUSION_SPACE  || cell_type[left] == DIFFUSION_MIX_SPACE)){
							neighbors[l*6 + counter] = left;
							counter++;
							neighbornum_X[l] = 1;					
						}
					}
					else	if (i == 1){
						int right = COORintoARRAY(i + 1, j, k);
						if (cell_type[right] == DIFFUSION_SPACE || cell_type[right] == DIFFUSION_MIX_SPACE){
							neighbors[l*6 + counter] = right;
							counter++;
							neighbornum_X[l] = 1;
						}
					}
					else	if (i == num_x_cell){
						int left = COORintoARRAY(i - 1, j, k);
						if (cell_type[left] == DIFFUSION_SPACE || cell_type[left] == DIFFUSION_MIX_SPACE){
							neighbors[l*6 + counter] = left;
							counter++;
							neighbornum_X[l] = 1;					
						}
					}
					//Y-direction
					if (j != 1 && j != num_y_cell){
						int forw = COORintoARRAY(i, j + 1, k);
						int back = COORintoARRAY(i, j - 1, k);
						if ((cell_type[forw] == DIFFUSION_SPACE || cell_type[forw] == DIFFUSION_MIX_SPACE) && (cell_type[back] == DIFFUSION_SPACE || cell_type[back] == DIFFUSION_MIX_SPACE)){
							neighbors[l*6 + counter] = forw;
							counter++;
							neighbors[l*6 + counter] = back;
							counter++; 
							neighbornum_Y[l] = 2;
						}
						else	if ((cell_type[forw] == DIFFUSION_SPACE || cell_type[forw] == DIFFUSION_MIX_SPACE) && (cell_type[back] != DIFFUSION_SPACE && cell_type[back] != DIFFUSION_MIX_SPACE)){
							neighbors[l*6 + counter] = forw;
							counter++;
							neighbornum_Y[l] = 1;					
						}
						else	if ((cell_type[forw] != DIFFUSION_SPACE && cell_type[forw] != DIFFUSION_MIX_SPACE) && (cell_type[back] == DIFFUSION_SPACE || cell_type[back] == DIFFUSION_MIX_SPACE)){
							neighbors[l*6 + counter] = back;
							counter++;
							neighbornum_Y[l] = 1;					
						}
					}
					else	if (j == 1){
						int forw = COORintoARRAY(i, j + 1, k);
						if (cell_type[forw] == DIFFUSION_SPACE || cell_type[forw] == DIFFUSION_MIX_SPACE){
							neighbors[l*6 + counter] = forw;
							counter++;
							neighbornum_Y[l] = 1;
						}
					}
					else	if (j == num_y_cell){
						int back = COORintoARRAY(i, j - 1, k);
						if (cell_type[back] == DIFFUSION_SPACE || cell_type[back] == DIFFUSION_MIX_SPACE){
							neighbors[l*6 + counter] = back;
							counter++;
							neighbornum_Y[l] = 1;					
						}
					}
					//Z-direction
					if (k != 1 && k != num_z_cell){
						int up = COORintoARRAY(i, j, k  + 1);
						int down = COORintoARRAY(i, j, k - 1);
						if ((cell_type[up] == DIFFUSION_SPACE || cell_type[up] == DIFFUSION_MIX_SPACE) && (cell_type[down] == DIFFUSION_SPACE || cell_type[down] == DIFFUSION_MIX_SPACE)){
							neighbors[l*6 + counter] = up;
							counter++;
							neighbors[l*6 + counter] = down;
							counter++; 
							neighbornum_Z[l] = 2;
						}
						else	if ((cell_type[up] == DIFFUSION_SPACE || cell_type[up] == DIFFUSION_MIX_SPACE) && (cell_type[down] != DIFFUSION_SPACE && cell_type[down] != DIFFUSION_MIX_SPACE)){
							neighbors[l*6 + counter] = up;
							counter++;
							neighbornum_Z[l] = 1;					
						}
						else	if ((cell_type[up] != DIFFUSION_SPACE && cell_type[up] != DIFFUSION_MIX_SPACE) && (cell_type[down] == DIFFUSION_SPACE || cell_type[down] == DIFFUSION_MIX_SPACE)){
							neighbors[l*6 + counter] = down;
							counter++;
							neighbornum_Z[l] = 1;					
						}
					}
					else	if (k == 1){
						int up = COORintoARRAY(i, j, k  + 1);
						if (cell_type[up] == DIFFUSION_SPACE || cell_type[up] == DIFFUSION_MIX_SPACE){
							neighbors[l*6 + counter] = up;
							counter++;
							neighbornum_Z[l] = 1;
						}
					}
					else	if (k == num_z_cell){
						int down = COORintoARRAY(i, j, k - 1);
						if (cell_type[down] == DIFFUSION_SPACE || cell_type[down] == DIFFUSION_MIX_SPACE){
							neighbors[l*6 + counter] = down;
							counter++;
							neighbornum_Z[l] = 1;					
						}
					}
				}	
				neighbornum[l] = neighbornum_X[l] + neighbornum_Y[l] + neighbornum_Z[l];

				SVcenter[l].x = (sv_size/2.0 + i*sv_size);
				SVcenter[l].y = (sv_size/2.0 + j*sv_size);
				SVcenter[l].z = (sv_size/2.0 + k*sv_size);
				SVcenter[l] = COORtransfer(SVcenter[l]);
			}
		}
	}
}

void AdiffusionBar(int* cell_type, int* neighbors, int* neighbornum, int* neighbornum_X, int* neighbornum_Y, int* neighbornum_Z, float3* SVcenter){
	int i, j, k;
	for (i = 1; i < num_x_cell + 1; i++){
		for (j = 1; j < num_y_cell + 1; j++){
			for (k = 1; k < num_z_cell + 1; k++){
				int counter = 0;
				int l = COORintoARRAY(i, j, k);
				if (cell_type[l] == AUA_SPACE || cell_type[l] == DIFFUSION_MIX_SPACE){
					//X-direction
					if (i != 1 && i != num_x_cell){
						int right = COORintoARRAY(i + 1, j, k);
						int left = COORintoARRAY(i - 1, j, k);
						if ((cell_type[right] == AUA_SPACE || cell_type[right] == DIFFUSION_MIX_SPACE) && (cell_type[left] == AUA_SPACE  || cell_type[left] == DIFFUSION_MIX_SPACE)){
							neighbors[l*6 + counter] = right;
							counter++;
							neighbors[l*6 + counter] = left;
							counter++; 
							neighbornum_X[l] = 2;
						}
						else	if ((cell_type[right] == AUA_SPACE || cell_type[right] == DIFFUSION_MIX_SPACE) && (cell_type[left] != AUA_SPACE  && cell_type[left] != DIFFUSION_MIX_SPACE)){
							neighbors[l*6 + counter] = right;
							counter++;
							neighbornum_X[l] = 1;					
						}
						else	if ((cell_type[right] != AUA_SPACE  && cell_type[right] != DIFFUSION_MIX_SPACE) && (cell_type[left] == AUA_SPACE  || cell_type[left] == DIFFUSION_MIX_SPACE)){
							neighbors[l*6 + counter] = left;
							counter++;
							neighbornum_X[l] = 1;					
						}
					}
					else	if (i == 1){
						int right = COORintoARRAY(i + 1, j, k);
						if (cell_type[right] == AUA_SPACE || cell_type[right] == DIFFUSION_MIX_SPACE){
							neighbors[l*6 + counter] = right;
							counter++;
							neighbornum_X[l] = 1;
						}
					}
					else	if (i == num_x_cell){
						int left = COORintoARRAY(i - 1, j, k);
						if (cell_type[left] == AUA_SPACE || cell_type[left] == DIFFUSION_MIX_SPACE){
							neighbors[l*6 + counter] = left;
							counter++;
							neighbornum_X[l] = 1;					
						}
					}
					//Y-direction
					if (j != 1 && j != num_y_cell){
						int forw = COORintoARRAY(i, j + 1, k);
						int back = COORintoARRAY(i, j - 1, k);
						if ((cell_type[forw] == AUA_SPACE || cell_type[forw] == DIFFUSION_MIX_SPACE) && (cell_type[back] == AUA_SPACE || cell_type[back] == DIFFUSION_MIX_SPACE)){
							neighbors[l*6 + counter] = forw;
							counter++;
							neighbors[l*6 + counter] = back;
							counter++; 
							neighbornum_Y[l] = 2;
						}
						else	if ((cell_type[forw] == AUA_SPACE || cell_type[forw] == DIFFUSION_MIX_SPACE) && (cell_type[back] != AUA_SPACE && cell_type[back] != DIFFUSION_MIX_SPACE)){
							neighbors[l*6 + counter] = forw;
							counter++;
							neighbornum_Y[l] = 1;					
						}
						else	if ((cell_type[forw] != AUA_SPACE && cell_type[forw] != DIFFUSION_MIX_SPACE) && (cell_type[back] == AUA_SPACE || cell_type[back] == DIFFUSION_MIX_SPACE)){
							neighbors[l*6 + counter] = back;
							counter++;
							neighbornum_Y[l] = 1;					
						}
					}
					else	if (j == 1){
						int forw = COORintoARRAY(i, j + 1, k);
						if (cell_type[forw] == AUA_SPACE || cell_type[forw] == DIFFUSION_MIX_SPACE){
							neighbors[l*6 + counter] = forw;
							counter++;
							neighbornum_Y[l] = 1;
						}
					}
					else	if (j == num_y_cell){
						int back = COORintoARRAY(i, j - 1, k);
						if (cell_type[back] == AUA_SPACE || cell_type[back] == DIFFUSION_MIX_SPACE){
							neighbors[l*6 + counter] = back;
							counter++;
							neighbornum_Y[l] = 1;					
						}
					}
					//Z-direction
					if (k != 1 && k != num_z_cell){
						int up = COORintoARRAY(i, j, k  + 1);
						int down = COORintoARRAY(i, j, k - 1);
						if ((cell_type[up] == AUA_SPACE || cell_type[up] == DIFFUSION_MIX_SPACE) && (cell_type[down] == AUA_SPACE || cell_type[down] == DIFFUSION_MIX_SPACE)){
							neighbors[l*6 + counter] = up;
							counter++;
							neighbors[l*6 + counter] = down;
							counter++; 
							neighbornum_Z[l] = 2;
						}
						else	if ((cell_type[up] == AUA_SPACE || cell_type[up] == DIFFUSION_MIX_SPACE) && (cell_type[down] != AUA_SPACE && cell_type[down] != DIFFUSION_MIX_SPACE)){
							neighbors[l*6 + counter] = up;
							counter++;
							neighbornum_Z[l] = 1;					
						}
						else	if ((cell_type[up] != AUA_SPACE && cell_type[up] != DIFFUSION_MIX_SPACE) && (cell_type[down] == AUA_SPACE || cell_type[down] == DIFFUSION_MIX_SPACE)){
							neighbors[l*6 + counter] = down;
							counter++;
							neighbornum_Z[l] = 1;					
						}
					}
					else	if (k == 1){
						int up = COORintoARRAY(i, j, k  + 1);
						if (cell_type[up] == AUA_SPACE || cell_type[up] == DIFFUSION_MIX_SPACE){
							neighbors[l*6 + counter] = up;
							counter++;
							neighbornum_Z[l] = 1;
						}
					}
					else	if (k == num_z_cell){
						int down = COORintoARRAY(i, j, k - 1);
						if (cell_type[down] == AUA_SPACE || cell_type[down] == DIFFUSION_MIX_SPACE){
							neighbors[l*6 + counter] = down;
							counter++;
							neighbornum_Z[l] = 1;					
						}
					}
				}	
				neighbornum[l] = neighbornum_X[l] + neighbornum_Y[l] + neighbornum_Z[l];

				SVcenter[l].x = (sv_size/2.0 + i*sv_size);
				SVcenter[l].y = (sv_size/2.0 + j*sv_size);
				SVcenter[l].z = (sv_size/2.0 + k*sv_size);
				SVcenter[l] = COORtransfer(SVcenter[l]);
			}
		}
	}
}