#pragma once
#include "vector_types.h"

typedef struct {
	float x, y, z;
	int w;
    float c, rad;
} kin;

typedef struct {
	int N;
	float3* h_r;
	float3* h_f;
	int* h_type;
    float3* d_r;
    float3* d_f;
    int* d_type;
    float* h_length;
    float* d_length;
    int* h_connector;
    int* d_connector;
    float* h_att_l;
    float* d_att_l;
    int* h_att_n;
    int* d_att_n;
    float* h_att_f;
    float* d_att_f;
    float* h_link_l;
    float* d_link_l;
    int* h_connect_1;
    int* d_connect_1;
    int* h_connect_2;
    int* d_connect_2;
    float* h_len_1;
    float* d_len_1;
    float* h_len_2;
    float* d_len_2;
    float3* h_dom_center;
    float3* d_dom_center;
    float* h_f_out;
    float* d_f_out;
} MDSystem;

typedef struct{
	float3* r;	//xyz-coordinate
	float length;		//quantity of MT-beads
	int state;	//grow or shortenin
    int attach;
	int rand;
    int NDCn;
    int cyl_N;
    int kinesin;
} MT;

typedef struct{
	float3* r;
	int n;
	int* state;
	int centerID;
	int* type;
	int* cell;
	int* react;
	int* rand;
    float* cos;
    float* radius;
    int* attach;
    int* attachN;
} KT;
