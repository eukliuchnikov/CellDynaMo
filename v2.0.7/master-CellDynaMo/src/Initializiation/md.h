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
} MDSystem;

typedef struct{
	float3* r;	//xyz-coordinate
	float length;		//quantity of MT-beads
	int state;	//grow or shortenin
    int attach;
	int rand;
    int NDCn;
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
} KT;
