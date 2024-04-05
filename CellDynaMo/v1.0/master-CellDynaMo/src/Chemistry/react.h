#pragma once
#include "../Initializiation/init.h"
#include "../Math/mat.h"
#include "diffusion.h"
#include "phos-dephos.h"
#include "attach-detach.h"
#include "mt_growth.h"
#include "mt_short.h"
#include "mt_cat.h"
#include "mt_resc.h"

void reaction(float time, FILE* outp_att, int reactnum[35], int* Nab, int* Naa, int* Nn0, int* Nn1, int* Nn2, int* Nn3, int* Nn4, int* Nn5, int* Nn6, int* Nn7, int* Nn0_A, int* Nn1_A, int* Nn2_A, int* Nn3_A, int* Nn4_A, int* Nn5_A, int* Nn6_A, int* Nn7_A, int* Nmtpol, int* Nmtdepol, int* Nmtdet, int & minimiz);
