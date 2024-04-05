/*
 * visualization module
 *
 *  Created on: Mar 15, 2018
 *	Edited on:	Mar 30, 2019
 *      Author: kliuchnikov
 */

#include "visualization.h"

#define cimg_display 0

#include "../TCO_Writers/CImg.h"

using namespace cimg_library;

void membrane_visualization(){ 
	if (MEM_SHAPE == MEM_SPHERE || MEM_SHAPE == MEM_ELIPSOID){ 
		float3 cell_center;
		float dtheta = 10.0;
		float dphi = 10.0;
		XYZ membraneXYZ;
		membraneXYZ.atomCount = (int)((2.0*180.0/dtheta + 1)*(180.0/dphi + 1));
		membraneXYZ.atoms = (XYZAtom*)calloc(membraneXYZ.atomCount, sizeof(XYZAtom));
		cell_center.x = 0; 
		cell_center.y = 0; 
		cell_center.z = 0;
		float theta;
		float phi;
		int n = 0;
		for (theta = 0; theta < 2*M_PI; theta += 10.0*M_PI/180.0){
			for (phi = 0; phi < M_PI; phi += 10.0*M_PI/180.0){
				membraneXYZ.atoms[n].name = 'M';
				membraneXYZ.atoms[n].x = cell_center.x + ELLIPSE_A*sin(theta)*cos(phi);
				membraneXYZ.atoms[n].y = cell_center.y + ELLIPSE_B*sin(theta)*sin(phi);
				membraneXYZ.atoms[n].z = cell_center.z + ELLIPSE_C*cos(theta);
				n++;
			}
		}
		membraneXYZ.atomCount = n;
		const char* xyz_name = config.output_membrane_xyz.c_str();
	  	char s2[config.output_membrane_xyz.size()+1];
	  	strcpy(s2, xyz_name);
		writeXYZ(xyz_name, &membraneXYZ);
		PSF psfMEM;
		psfMEM.natom = membraneXYZ.atomCount;
		psfMEM.atoms = (PSFAtom*)calloc(psfMEM.natom, sizeof(PSFAtom));
		psfMEM.nbond = psfMEM.natom + 666;
		psfMEM.bonds = (PSFBond*)calloc(psfMEM.nbond, sizeof(PSFBond));
		psfMEM.ntheta = 0;
		psfMEM.nphi = 0;
		psfMEM.nimphi = 0;
		psfMEM.ncmap = 0;
		psfMEM.nnb = 0;
		int o, w, e;
		for (o = 0; o < psfMEM.natom; o++){
			psfMEM.atoms[o].id = o + 1;
			if (o < psfMEM.natom){
				sprintf(psfMEM.atoms[o].name, "MEM");
				sprintf(psfMEM.atoms[o].type, "MEM");
				sprintf(psfMEM.atoms[o].segment, "MEM");
			}
			sprintf(psfMEM.atoms[o].resName, "Mem");
			psfMEM.atoms[o].resid = o + 1;
			psfMEM.atoms[o].m = 1.0;
			psfMEM.atoms[o].q = 0.0;
		}
		w = 0;
		for (e = 0; e < psfMEM.natom; e++){
			if ((e + 1)%19 != 0){
				psfMEM.bonds[w].i = e + 1;
				psfMEM.bonds[w].j = e + 2;
				w++;
				psfMEM.bonds[w].i = e + 1;
				psfMEM.bonds[w].j = e + 20;
				w++;
			}
			else{
				psfMEM.bonds[w].i = e + 1;
				psfMEM.bonds[w].j = e + 20;
				w++;	
			}
		}
		char filenameMEM[1024];
		const char* psf_name = config.output_membrane_psf.c_str();
	  	char s3[config.output_membrane_psf.size()+1];
	  	strcpy(s3, psf_name);
		sprintf(filenameMEM, psf_name);
		writePSF(filenameMEM, &psfMEM);
	}
	else	if (MEM_SHAPE == MEM_CUBE || MEM_SHAPE == MEM_CUBOID){
		XYZ membraneXYZ;
		membraneXYZ.atomCount = 8;
		membraneXYZ.atoms = (XYZAtom*)calloc(membraneXYZ.atomCount, sizeof(XYZAtom));

		membraneXYZ.atoms[0].name = 'M';
		membraneXYZ.atoms[0].x = -ELLIPSE_A/2;
		membraneXYZ.atoms[0].y = -ELLIPSE_B/2;
		membraneXYZ.atoms[0].z = -ELLIPSE_C/2;

		membraneXYZ.atoms[1].name = 'M';
		membraneXYZ.atoms[1].x = -ELLIPSE_A/2;
		membraneXYZ.atoms[1].y = -ELLIPSE_B/2;
		membraneXYZ.atoms[1].z = ELLIPSE_C/2;

		membraneXYZ.atoms[2].name = 'M';
		membraneXYZ.atoms[2].x = -ELLIPSE_A/2;
		membraneXYZ.atoms[2].y = ELLIPSE_B/2;
		membraneXYZ.atoms[2].z = ELLIPSE_C/2;

		membraneXYZ.atoms[3].name = 'M';
		membraneXYZ.atoms[3].x = -ELLIPSE_A/2;
		membraneXYZ.atoms[3].y = ELLIPSE_B/2;
		membraneXYZ.atoms[3].z = -ELLIPSE_C/2;

		membraneXYZ.atoms[4].name = 'M';
		membraneXYZ.atoms[4].x = ELLIPSE_A/2;
		membraneXYZ.atoms[4].y = ELLIPSE_B/2;
		membraneXYZ.atoms[4].z = -ELLIPSE_C/2;

		membraneXYZ.atoms[5].name = 'M';
		membraneXYZ.atoms[5].x = ELLIPSE_A/2;
		membraneXYZ.atoms[5].y = -ELLIPSE_B/2;
		membraneXYZ.atoms[5].z = -ELLIPSE_C/2;

		membraneXYZ.atoms[6].name = 'M';
		membraneXYZ.atoms[6].x = ELLIPSE_A/2;
		membraneXYZ.atoms[6].y = -ELLIPSE_B/2;
		membraneXYZ.atoms[6].z = ELLIPSE_C/2;

		membraneXYZ.atoms[7].name = 'M';
		membraneXYZ.atoms[7].x = ELLIPSE_A/2;
		membraneXYZ.atoms[7].y = ELLIPSE_B/2;
		membraneXYZ.atoms[7].z = ELLIPSE_C/2;

		const char* xyz_name = config.output_membrane_xyz.c_str();
	  	char s2[config.output_membrane_xyz.size()+1];
	  	strcpy(s2, xyz_name);
		writeXYZ(xyz_name, &membraneXYZ);
		PSF psfMEM;
		psfMEM.natom = membraneXYZ.atomCount;
		psfMEM.atoms = (PSFAtom*)calloc(psfMEM.natom, sizeof(PSFAtom));
		psfMEM.nbond = 28;
		psfMEM.bonds = (PSFBond*)calloc(psfMEM.nbond, sizeof(PSFBond));
		psfMEM.ntheta = 0;
		psfMEM.nphi = 0;
		psfMEM.nimphi = 0;
		psfMEM.ncmap = 0;
		psfMEM.nnb = 0;
		int o, w, e;
		for (o = 0; o < psfMEM.natom; o++){
			psfMEM.atoms[o].id = o + 1;
			if (o < psfMEM.natom){
				sprintf(psfMEM.atoms[o].name, "MEM");
				sprintf(psfMEM.atoms[o].type, "MEM");
				sprintf(psfMEM.atoms[o].segment, "MEM");
			}
			sprintf(psfMEM.atoms[o].resName, "Mem");
			psfMEM.atoms[o].resid = o + 1;
			psfMEM.atoms[o].m = 1.0;
			psfMEM.atoms[o].q = 0.0;
		}
		w = 0;
		for (e = 0; e < 4; e++){
			for (int e2 = 0; e2 < 4; e2++){
				if (e < e2){
					psfMEM.bonds[w].i = e + 1;
					psfMEM.bonds[w].j = e2 + 1;
					w++;
				}
			}
		}
		for (e = 4; e < 8; e++){
			for (int e2 = 4; e2 < 8; e2++){
				if (e < e2){
					psfMEM.bonds[w].i = e + 1;
					psfMEM.bonds[w].j = e2 + 1;
					w++;
				}
			}
		}
		for (e = 0; e < 4; e += 3){
			for (int e2 = 4; e2 < 6; e2++){
				if (e < e2){
					psfMEM.bonds[w].i = e + 1;
					psfMEM.bonds[w].j = e2 + 1;
					w++;
				}
			}
		}
		for (e = 1; e < 3; e++){
			for (int e2 = 6; e2 < 8; e2++){
				if (e < e2){
					psfMEM.bonds[w].i = e + 1;
					psfMEM.bonds[w].j = e2 + 1;
					w++;
				}
			}
		}
		for (e = 0; e < 2; e++){
			for (int e2 = 5; e2 < 7; e2++){
				if (e < e2){
					psfMEM.bonds[w].i = e + 1;
					psfMEM.bonds[w].j = e2 + 1;
					w++;
				}
			}
		}
		for (e = 2; e < 4; e++){
			for (int e2 = 4; e2 < 8; e2 += 3){
				if (e < e2){
					psfMEM.bonds[w].i = e + 1;
					psfMEM.bonds[w].j = e2 + 1;
					w++;
				}
			}
		}
		char filenameMEM[1024];
		const char* psf_name = config.output_membrane_psf.c_str();
	  	char s3[config.output_membrane_psf.size()+1];
	  	strcpy(s3, psf_name);
		sprintf(filenameMEM, psf_name);
		writePSF(filenameMEM, &psfMEM);
	}
}

int initPSF(KT* kt){
	int h, k;
	int kbonds = 0;
	PSF psf;
	psf.natom = mds.N;
	psf.atoms = (PSFAtom*)calloc(psf.natom, sizeof(PSFAtom));
	psf.nbond = POLE_COUNT*MT_NUM*(MAX_MT_LENGTH - 1) + KIN_COUNT*(kt[0].n - 2) + KIN_COUNT/2;
	//psf.nbond = POLE_COUNT*MT_NUM*(MAX_MT_LENGTH - 1) + KIN_COUNT*4 + KIN_COUNT/2;
	psf.bonds = (PSFBond*)calloc(psf.nbond, sizeof(PSFBond));
	psf.ntheta = 0;
	psf.nphi = 0;
	psf.nimphi = 0;
	psf.ncmap = 0;
	psf.nnb = 0;
	int n_chrom = 2*floor((CONTOUR_LENGTH/2 - KIN_RADIUS)/(2*KIN_RADIUS) + 0.5);
	int i, j, m, mi, b, l;
	int x_num = CORONA_X;
	int y_num = floor(x_num/RATIO);
	for (i = 0; i < psf.natom; i++){
		psf.atoms[i].id = i + 1;
		if (i < POLE_COUNT*MT_NUM*MAX_MT_LENGTH){
			sprintf(psf.atoms[i].name, "MT");
			sprintf(psf.atoms[i].type, "MT");
			sprintf(psf.atoms[i].segment, "MT");
		}
		else	{
			for (int ks = 0; ks < KIN_COUNT; ks +=2){
				if (i < POLE_COUNT*MT_NUM*MAX_MT_LENGTH + ks*kt[0].n + x_num*y_num && i > POLE_COUNT*MT_NUM*MAX_MT_LENGTH + ks*kt[0].n - 1){
					sprintf(psf.atoms[i].name, "KTL");
					sprintf(psf.atoms[i].type, "KTL");
					sprintf(psf.atoms[i].segment, "KTL");
				}
				if (i < POLE_COUNT*MT_NUM*MAX_MT_LENGTH + (ks + 1)*kt[0].n - n_chrom - 1 && i >= POLE_COUNT*MT_NUM*MAX_MT_LENGTH + ks*kt[0].n + x_num*y_num){
					sprintf(psf.atoms[i].name, "KTLS");
					sprintf(psf.atoms[i].type, "KTLS");
					sprintf(psf.atoms[i].segment, "KTLS");
				}
				if (i == POLE_COUNT*MT_NUM*MAX_MT_LENGTH + (ks + 1)*kt[0].n - n_chrom - 1){
					sprintf(psf.atoms[i].name, "KTLC");
					sprintf(psf.atoms[i].type, "KTLC");
					sprintf(psf.atoms[i].segment, "KTLC");
				}
				if (i > POLE_COUNT*MT_NUM*MAX_MT_LENGTH + (ks + 1)*kt[0].n - n_chrom - 1 && i < POLE_COUNT*MT_NUM*MAX_MT_LENGTH + (ks + 1)*kt[0].n){
					sprintf(psf.atoms[i].name, "CH");
					sprintf(psf.atoms[i].type, "CH");
					sprintf(psf.atoms[i].segment, "CH");
				}
				if (i < POLE_COUNT*MT_NUM*MAX_MT_LENGTH + (ks + 1)*kt[0].n + x_num*y_num && i > POLE_COUNT*MT_NUM*MAX_MT_LENGTH + (ks + 1)*kt[0].n - 1){
					sprintf(psf.atoms[i].name, "KTR");
					sprintf(psf.atoms[i].type, "KTR");
					sprintf(psf.atoms[i].segment, "KTR");
				}
				if (i < POLE_COUNT*MT_NUM*MAX_MT_LENGTH + (ks + 2)*kt[0].n - n_chrom - 1 && i >= POLE_COUNT*MT_NUM*MAX_MT_LENGTH + (ks + 1)*kt[0].n + x_num*y_num){
					sprintf(psf.atoms[i].name, "KTRS");
					sprintf(psf.atoms[i].type, "KTRS");
					sprintf(psf.atoms[i].segment, "KTRS");
				}
				if (i == POLE_COUNT*MT_NUM*MAX_MT_LENGTH + (ks + 2)*kt[0].n - n_chrom - 1){
					sprintf(psf.atoms[i].name, "KTRC");
					sprintf(psf.atoms[i].type, "KTRC");
					sprintf(psf.atoms[i].segment, "KTRC");
				}
				if (i > POLE_COUNT*MT_NUM*MAX_MT_LENGTH + (ks + 2)*kt[0].n - n_chrom - 1 && i < POLE_COUNT*MT_NUM*MAX_MT_LENGTH + (ks + 2)*kt[0].n){
					sprintf(psf.atoms[i].name, "CH");
					sprintf(psf.atoms[i].type, "CH");
					sprintf(psf.atoms[i].segment, "CH");
				}
			}
		}
		if (i == mds.N - 1 || i == mds.N - 2){
			sprintf(psf.atoms[i].name, "Pole");
			sprintf(psf.atoms[i].type, "Pole");
			sprintf(psf.atoms[i].segment, "Pole");
		}
		sprintf(psf.atoms[i].resName, "Cell");
		psf.atoms[i].resid = i + 1;
		psf.atoms[i].m = 1.0;
		psf.atoms[i].q = 0.0;
	}
	float r0 = 2.0*MT_RADIUS;
//MT's bonds
	b = 0;
	for (l = 0; l < POLE_COUNT; l++){
		for (m = 0; m < MT_NUM; m++){
			for (mi = 0; mi < MAX_MT_LENGTH - 1; mi++){
				psf.bonds[b].i = l*MT_NUM*MAX_MT_LENGTH + m*MAX_MT_LENGTH + mi + 1;
				psf.bonds[b].j = l*MT_NUM*MAX_MT_LENGTH + m*MAX_MT_LENGTH + mi + 2;
				b++;
			}
		}
	}
	for (int ks = 0; ks < KIN_COUNT; ks ++){
		for (k = 0; k < kt[0].n - n_chrom - 2; k++){
			psf.bonds[b].i = POLE_COUNT*MT_NUM*MAX_MT_LENGTH + ks*kt[0].n + k + 1;
			psf.bonds[b].j = POLE_COUNT*MT_NUM*MAX_MT_LENGTH + ks*kt[0].n + k + 2;
			b++;
		}
	}
	//CENTERS OF MASS
	for (int ks = 0; ks < KIN_COUNT; ks += 2){
		psf.bonds[b].i = POLE_COUNT*MT_NUM*MAX_MT_LENGTH + (ks + 1)*kt[0].n - n_chrom;
		psf.bonds[b].j = POLE_COUNT*MT_NUM*MAX_MT_LENGTH + (ks + 2)*kt[0].n - n_chrom;
		b++;
	}
	if (CHROM_ARMS == 1){
		for (int ks = 0; ks < KIN_COUNT; ks ++){
			psf.bonds[b].i = POLE_COUNT*MT_NUM*MAX_MT_LENGTH + (ks + 1)*kt[0].n - n_chrom;
			psf.bonds[b].j = POLE_COUNT*MT_NUM*MAX_MT_LENGTH + (ks + 1)*kt[0].n - n_chrom + 1;
			b++;
			psf.bonds[b].i = POLE_COUNT*MT_NUM*MAX_MT_LENGTH + (ks + 1)*kt[0].n - n_chrom;
			psf.bonds[b].j = POLE_COUNT*MT_NUM*MAX_MT_LENGTH + (ks + 1)*kt[0].n - n_chrom/2 + 1;
			b++;
			for (int cn = 0; cn < n_chrom/2 - 1; cn ++){
				psf.bonds[b].i = POLE_COUNT*MT_NUM*MAX_MT_LENGTH + (ks + 1)*kt[0].n - cn;
				psf.bonds[b].j = POLE_COUNT*MT_NUM*MAX_MT_LENGTH + (ks + 1)*kt[0].n - cn - 1;
				b++;
			}
			for (int cn = n_chrom/2; cn < n_chrom - 1; cn ++){
				psf.bonds[b].i = POLE_COUNT*MT_NUM*MAX_MT_LENGTH + (ks + 1)*kt[0].n - cn;
				psf.bonds[b].j = POLE_COUNT*MT_NUM*MAX_MT_LENGTH + (ks + 1)*kt[0].n - cn - 1;
				b++;
			}
		}
	}
	b--;
	char filename[1024];
	const char* psf_name = config.outputpsf.c_str();
  	char s3[config.outputpsf.size()+1];
  	strcpy(s3, psf_name);
	sprintf(filename, psf_name);
	writePSF(filename, &psf);
}

void DCD_step(DCD dcd){
	int i;
	for (i = 0; i < mds.N; i++){
		dcd.frame.X[i] = mds.h_r[i].x;
		dcd.frame.Y[i] = mds.h_r[i].y;
		dcd.frame.Z[i] = mds.h_r[i].z;
	}
	dcdWriteFrame(dcd);
}

void map(int* cell_type, int* aurora_B_conc, int* aurora_A_conc, int* NDC80_0_conc, char* gillespie, int number){
	int img_scale = 10;
	CImg<unsigned char>mapImg(num_x_cell*img_scale, num_y_cell*img_scale, 1, 3);
	sprintf(gillespie, "output/map/map%d.png", number);
	int i, j, k, h;
	for (i = 1; i < num_x_cell + 1; i++){
		for (j = 1; j < num_y_cell + 1; j++){
			for (k = 1; k < num_z_cell + 1; k++){
				int3 cell;
				cell.x = i;
				cell.y = j;
				cell.z = k;
				float3 r;
				r = get_coord(cell);
				h = COORintoARRAY(i, j, k);
				if (k == (num_z_cell + 1)/2){
					for(int s = 0; s < img_scale; s++){
						for(int d = 0; d < img_scale; d++){
							int y = (j - 1)*img_scale + s;
							int x = (i - 1)*img_scale + d;
							
							if (cell_type[h] == OUTER_SPACE){
								mapImg(x, y, 0, 0) = 0;
								mapImg(x, y, 0, 1) = 0;
								mapImg(x, y, 0, 2) = 0;
							}
							if (cell_type[h] == INNER_SPACE){
								mapImg(x, y, 0, 0) = 255;
								mapImg(x, y, 0, 1) = 255;
								mapImg(x, y, 0, 2) = 255;
							}
							if (cell_type[h] == DIFFUSION_SPACE){
								mapImg(x, y, 0, 0) = 150;
								mapImg(x, y, 0, 1) = 230;
								mapImg(x, y, 0, 2) = 112;
							}
							if (aurora_B_conc[h] > 0){
								mapImg(x, y, 0, 0) = 255;
								mapImg(x, y, 0, 1) = 0;
								mapImg(x, y, 0, 2) = 0;
							}
							if (aurora_A_conc[h] > 0){
								mapImg(x, y, 0, 0) = 255;
								mapImg(x, y, 0, 1) = 255;
								mapImg(x, y, 0, 2) = 0;
							}
							if (aurora_A_conc[h] > 0 && aurora_B_conc[h] > 0){
								mapImg(x, y, 0, 0) = 255;
								mapImg(x, y, 0, 1) = 128;
								mapImg(x, y, 0, 2) = 0;
							}
							if (NDC80_0_conc[h] > 0){
								mapImg(x, y, 0, 0) = 0;
								mapImg(x, y, 0, 1) = 0;
								mapImg(x, y, 0, 2) = 255;
							}
							/*if (pow((r.x - 0.5*float(XBOX)), 2)/pow(ELLIPSE_A, 2) + pow((r.y - 0.5*float(YBOX)), 2)/pow(ELLIPSE_B, 2) + pow((r.z - 0.5*float(ZBOX)), 2)/pow(ELLIPSE_C, 2) > 1){
								mapImg(x, y, 0, 0) = 0;
								mapImg(x, y, 0, 1) = 0;
								mapImg(x, y, 0, 2) = 0;
							}*/
						}
					}
				}
			}
		}
	}
	mapImg.rotate(90);
	mapImg.transpose();
	mapImg.save_png(gillespie);
}

void Au_out(float time, int* cell_type, FILE* outp, int* aurora_B_conc, int* aurora_A_conc){
	outp = fopen("output/dat/aurora.dat", "a");
	fprintf(outp, "*******************START*******************\n");
	fprintf(outp, "%f\n", time);
	fprintf(outp, "x\ty\tz\tcell_type\tAurora_B\tAurora_A\n");
	int i, j, k, h;
	for (i = 1; i < num_x_cell + 1; i++){
		for (j = 1; j < num_y_cell + 1; j++){
			for (k = 1; k < num_z_cell + 1; k++){
				h = COORintoARRAY(i, j, k);
				if (aurora_B_conc[h] != 0 || aurora_A_conc[h] != 0){
					fprintf(outp, "%d\t%d\t%d\t%d\t%d\t%d\n", i, j, k, cell_type[h], aurora_B_conc[h], aurora_A_conc[h]);
				}
			}
		}
	}
	fprintf(outp, "********************END********************\n");
	fclose(outp);
}

void ndc_out(float time, FILE* outp, int* NDC80_0_conc, int* NDC80_1_conc, int* NDC80_2_conc, int* NDC80_3_conc, int* NDC80_4_conc, int* NDC80_5_conc, int* NDC80_6_conc, int* NDC80_7_conc){
	outp = fopen("output/dat/ndc.dat", "a");
	fprintf(outp, "*******************START*******************\n");
	fprintf(outp, "%f\n", time);
	fprintf(outp, "x\ty\tz\tndc0\tndc1\tndc2\tndc3\tndc4\tndc5\tndc6\tndc7\n");
	int i, j, k, h;
	for (i = 1; i < num_x_cell + 1; i++){
		for (j = 1; j < num_y_cell + 1; j++){
			for (k = 1; k < num_z_cell + 1; k++){
				h = COORintoARRAY(i, j, k);
				if (NDC80_0_conc[h] != 0 || NDC80_1_conc[h] != 0 || NDC80_2_conc[h] != 0 || NDC80_3_conc[h] != 0 || NDC80_4_conc[h] != 0 || NDC80_5_conc[h] != 0 || NDC80_6_conc[h] != 0 || NDC80_7_conc[h] != 0){
					fprintf(outp, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", i, j, k, NDC80_0_conc[h], NDC80_1_conc[h], NDC80_2_conc[h], NDC80_3_conc[h], NDC80_4_conc[h], NDC80_5_conc[h], NDC80_6_conc[h], NDC80_7_conc[h]);
				}
			}
		}
	}
	fprintf(outp, "********************END********************\n");
	fclose(outp);
}

void att_out(float time, FILE* outp, int triger, int mt_id, int kt_id){
	outp = fopen("output/dat/att.dat", "a");
	if (triger == ATTACH){
		fprintf(outp, "%f\t%d\t%d\tATTACHMENT\n", time, mt_id, kt_id);
	}
	if (triger == DETACH){
		fprintf(outp, "%f\t%d\t%d\tDETACHMENT\n", time, mt_id, kt_id);
	}
	fclose(outp);
}