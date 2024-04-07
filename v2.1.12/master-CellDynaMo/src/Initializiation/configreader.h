#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

using namespace std;

struct Config{
    float acc_dt;
    int device;
    int rseed;
    string traj_name;
	float total_t;
	float sv_size1;
	float aurora_diff;
	float k_phos;
	float k_dephos;
	float k_attach;
	float k_detach_0;
	float k_detach_1;
	float k_detach_2;
	float k_detach_3;
	float k_detach_4;
	float k_detach_5;
	float k_detach_6;
	float k_detach_7;
	float k_growth;
	float k_shorten;
	float k_catastrophe;
	float k_rescue;	
	float aurora_conc;
	float ph_conc;
	float tub_conc;
	float dt;
	float lp;
	float k_spring_mt;
	float k_spring_kt;
	float k_stiff_mem;
	float eps_rep;
	float eps_att;
	float temp;
	float viscosity;
	float ellipse_a;
	float ellipse_b;
	float ellipse_c;
	float pole_distance;
	float pole_radius;
	int pole_count;
	int mt_num;
	int max_mt_length;
	int start_length;
	float mt_radius;
	float kin_radius;
	float beta;
	int kin_count;
    string arm_cond;
    int chrom_arms;
    float contour_length;
	float distance_ndc;
	int maxHarmonicKinPerMonomer1;

	float* kinetochore_move_x;
	float* kinetochore_move_y;
	float* kinetochore_move_z;
	float* x_angle;
	float* y_angle;
	float* z_angle;
    string kin_shape_cond;
    int kin_shape;
    float s_corona;
    int corona_x_num;
    float ratio;

	float ndc_length;
    int ndc_per_mt;
	string outputpsf;
	string outputdcd;
	string output_membrane_psf;
	string output_membrane_xyz;
	float dcd_freq;
    string map_cond;
    int map;
    string shape_cond;
    int shape;
};

extern Config config;

void loadConfig(Config& config);
