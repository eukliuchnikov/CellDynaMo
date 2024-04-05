#include "configreader.h"
#include <cuda_runtime_api.h>
#include <cuda.h>

#define DEVICE                      config.device
#define rseed                       config.rseed

#define TOTAL_T                     config.total_t
#define sv_size                     config.sv_size1
//#define tau        	                config.sv_size1*config.sv_size1/(2*config.acc_dt*config.k_shorten)  //timestep, s
#define tau        	                1.0/(config.acc_dt*config.k_shorten)  //timestep, s
#define trans                       pow(10, 7)/6.02  //transfer 1/(s*muM) to nm3/s

#define DaB                         config.aurora_diff
#define K_PHOS                      config.k_phos
#define K_DEPHOS                    config.k_dephos
#define K_ATTACH                    config.k_attach
#define K_DETACH_0                  config.k_detach_0
#define K_DETACH_1                  config.k_detach_1
#define K_DETACH_2                  config.k_detach_2
#define K_DETACH_3                  config.k_detach_3
#define K_DETACH_4                  config.k_detach_4
#define K_DETACH_5                  config.k_detach_5
#define K_DETACH_6                  config.k_detach_6
#define K_DETACH_7                  config.k_detach_7
#define K_GROWTH                    config.k_growth
#define K_SHORTEN                   config.k_shorten
#define K_CAT                       config.k_catastrophe
#define K_RESCUE                    config.k_rescue
#define AU_CONC                     config.aurora_conc
#define PH_CONC                     config.ph_conc
#define TUB_CONC                    config.tub_conc
#define TIMESTEP                    config.dt
#define K_STEFF_M                   config.k_stiff_mem
#define VISCOSITY                   config.viscosity
#define ELLIPSE_A                   config.ellipse_a
#define ELLIPSE_B                   config.ellipse_b
#define ELLIPSE_C                   config.ellipse_c

#define KT_SHAPE                    config.kin_shape
#define KT_SPHERE                   1
#define KT_RECT                     0

#define POLE_DISTANCE               config.pole_distance
#define POLE_RADIUS                 config.pole_radius
#define POLE_COUNT                  2
#define MT_NUM                      config.mt_num
#define MAX_MT_LENGTH               config.max_mt_length
#define MIN_MT_LENGTH               5
#define START_LENGTH                config.start_length
#define MT_RADIUS                   config.mt_radius
#define KIN_RADIUS                  config.kin_radius
#define BETA                        config.beta
#define CORONA                      config.s_corona
#define CORONA_X                    config.corona_x_num
#define RATIO                       config.ratio
//#define NO_ATT_ZONE                 config.kin_radius - config.kin_radius/(1 - config.beta) + sqrtf(config.kin_radius/(1 - config.beta)*config.kin_radius/(1 - config.beta) - config.kin_radius*config.kin_radius)
#define NO_ATT_ZONE                 config.kin_radius - 200*(1 - config.beta)
#define KIN_COUNT                   2
#define CHROM_ARMS                  config.chrom_arms
#define CONTOUR_LENGTH              config.contour_length
#define DISTANCE_NDC                config.distance_ndc
#define maxHarmonicKinPerMonomer    config.maxHarmonicKinPerMonomer1
#define NDC_LENGTH                  config.ndc_length
#define NDC_PER_MT                  config.ndc_per_mt
#define DCD_FREQ                    config.dcd_freq
#define MAP                         config.map
#define MEM_SHAPE                   config.shape
#define MEM_SPHERE                  1
#define MEM_CUBE                    2
#define MEM_ELIPSOID                3
#define MEM_CUBOID                  4
#define ARMOR                       config.armor
#define AURORA_A                    config.aurora
#define RAND_PLACE                  config.rand_place
#define CH_TER                      config.chrom_ter
#define RING                        config.ring
#define RING_RAD                    config.ring_rad

#define KB 						    8.314462e-3		//kJ/K/mol
#define WIDTH                       0.07            //Morse potential's width

#define MT_STATE_POL			    1				//type of plus-end with ability to assemble
#define MT_STATE_DEPOL			    2				//type of plus-end with ability to disassemble	
#define DETACH					    0
#define	ATTACH					    1

#define MT_NO		                0	//MTs in the (0,0,0) point
#define MT_REG		                11	//MTs which are not plus-end as well as minus-end
#define LEFT_POLE	                -1	
#define RIGHT_POLE	                1
#define PLUS_DET_INVALID            100
#define PLUS_DET	                10	//Plus-end that is not attached to NDC80
#define PLUS_ATT	                -10	//Attached plus-end	
#define LEFT_MINUS	                -12
#define RIGHT_MINUS	                12

#define LEFT_NDC	                -21
#define RIGHT_NDC	                21
#define LEFT_KT		                -20
#define RIGHT_KT	                20
#define CHROM                       22
#define SHELL_LEFT                  -23
#define SHELL_RIGHT                 23                       

#define NO_NDC                      10
#define NDC_0	                    0
#define NDC_1	                    1
#define NDC_2	                    2
#define NDC_3	                    3
#define NDC_4	                    4
#define NDC_5	                    5
#define NDC_6	                    6
#define NDC_7	                    7

#define XBOX                        2*ELLIPSE_A + sv_size
#define YBOX                        2*ELLIPSE_B + sv_size
#define ZBOX                        2*ELLIPSE_C + sv_size

#define num_x_cell                  int((XBOX)*pow(config.sv_size1, -1)) 
#define num_y_cell                  int((YBOX)*pow(config.sv_size1, -1))
#define num_z_cell                  int((ZBOX)*pow(config.sv_size1, -1))

#define SVn                         num_x_cell*num_y_cell*num_z_cell

#define OUTER_SPACE                 0
#define INNER_SPACE                 1
#define DIFFUSION_SPACE             2
#define AUA_SPACE                   3
#define DIFFUSION_MIX_SPACE         4

struct Param{
    float timestep;
    float lp;
    float k_mt;
    float k_kt;
    float e_rep;
    float e_att;
    float temper;
    float visc;
    float pole_r;
    int mt_n;
    int max_mt_n;
    float mt_r;
    float kt_r;
    int mHkPm;
    int pole_n;
    int kin_n;
    int chrom;
    float chrom_num;
    float ndc_l;
    float ndc_d;
    float zone;
    float big_r;
    int x_num;
    int y_num;
};

extern Param* h_parameters;
extern Param* d_parameters;

void parameters(Config& config);
