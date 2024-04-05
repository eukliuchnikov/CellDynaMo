/*
 * transfer some parameters to GPU
 *
 *  Created on: Mar 15, 2018
 *	Edited on:	Mar 30, 2019
 *      Author: kliuchnikov
 */

#include "parameters.h"

Param* h_parameters;
Param* d_parameters;

void parameters(Config& config){
	h_parameters = (Param*)calloc(1, sizeof(Param));
	h_parameters->lp = config.lp;
	h_parameters->timestep = config.dt;
	h_parameters->k_mt = config.k_spring_mt;
	h_parameters->k_kt = config.k_spring_kt;
	h_parameters->e_rep = config.eps_rep*41816.667;
	h_parameters->e_att = config.eps_att*18.518;
	h_parameters->temper = config.temp;
	h_parameters->pole_r = config.pole_radius;
	h_parameters->mt_n = config.mt_num;
	h_parameters->max_mt_n = MAX_MT_LENGTH;
	h_parameters->mt_r = config.mt_radius;
	h_parameters->kt_r = config.kin_radius;
	h_parameters->mHkPm = config.maxHarmonicKinPerMonomer1;
	h_parameters->pole_n = config.pole_count;
	h_parameters->kin_n = config.kin_count;
	h_parameters->chrom = config.chrom_arms;
	h_parameters->chrom_num = 2*floor((config.contour_length/2 - config.kin_radius)/(2*config.kin_radius) + 0.5);
	h_parameters->ndc_l = config.ndc_length;
	h_parameters->visc = config.viscosity;
	h_parameters->ndc_d = config.distance_ndc;
	h_parameters->zone = config.kin_radius - 200.0*(1.0 - config.beta);
	h_parameters->x_num = config.corona_x_num;
	h_parameters->y_num = config.corona_x_num/config.ratio;



	if (config.beta >= 0.925){
		h_parameters->big_r = 5000.0;
	}
	else	if (config.beta >= 0 && config.beta < 0.925){
		h_parameters->big_r = config.kin_radius/(1.0 - config.beta);
	}	
	cudaMalloc((void**)&d_parameters, 1*sizeof(Param));
    cudaMemcpy(d_parameters, h_parameters, 1*sizeof(Param), cudaMemcpyHostToDevice);
}
