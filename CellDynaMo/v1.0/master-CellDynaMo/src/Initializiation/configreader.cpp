#include "configreader.h"

Config config;

void loadConfig(Config& config) {
    ifstream fin("config/conf.conf");
    string line;
    int chrom_n;
    while (getline(fin, line)) {
        istringstream sin(line.substr(line.find(" ") + 1));
        if (line.find("device") != -1)
            sin >> config.device;
        else if (line.find("rseed") != -1)
            sin >> config.rseed;
        else if (line.find("traj_name ") != -1)
            sin >> config.traj_name;
        else if (line.find("total_t") != -1)
            sin >> config.total_t;
        else if (line.find("dt_accuracy") != -1)
            sin >> config.acc_dt;
        else if (line.find("kin_count") != -1){
            sin >> config.kin_count;
            chrom_n = config.kin_count/2;
            config.kinetochore_move_x = (float*)calloc(chrom_n, sizeof(float));
            config.kinetochore_move_y = (float*)calloc(chrom_n, sizeof(float));
            config.kinetochore_move_z = (float*)calloc(chrom_n, sizeof(float));
            config.x_angle = (float*)calloc(chrom_n, sizeof(float));
            config.y_angle = (float*)calloc(chrom_n, sizeof(float));
            config.z_angle = (float*)calloc(chrom_n, sizeof(float));
        }
        if (line.find("outputpsf") != -1){
            sin >> config.outputpsf;
        	config.outputpsf.replace(config.outputpsf.find("<traj_name>"),11,config.traj_name);
        }
        else if (line.find("outputdcd") != -1){
            sin >> config.outputdcd;
            config.outputdcd.replace(config.outputdcd.find("<traj_name>"),11,config.traj_name);
        }
        else if (line.find("output_membrane_psf") != -1){
            sin >> config.output_membrane_psf;
            config.output_membrane_psf.replace(config.output_membrane_psf.find("<traj_name>"),11,config.traj_name);
        }
        else if (line.find("output_membrane_xyz") != -1){
            sin >> config.output_membrane_xyz;
            config.output_membrane_xyz.replace(config.output_membrane_xyz.find("<traj_name>"),11,config.traj_name);
        }
        else if (line.find("dcd_freq") != -1)
            sin >> config.dcd_freq;
        else if (line.find("contour_length") != -1)
            sin >> config.contour_length;
        else if (line.find("implicit_map") != -1){
            sin >> config.map_cond;
        	if (config.map_cond == "yes"){
				config.map = 1;
			}
			else	if (config.map_cond == "no"){
				config.map = 0;
			}
		}
        else if (line.find("chrom_arms") != -1){
            sin >> config.arm_cond;
            if (config.arm_cond == "yes"){
                config.chrom_arms = 1;
            }
            else    if (config.arm_cond == "no"){
                config.chrom_arms = 0;
            }
        }
        else if (line.find("armor") != -1){
            sin >> config.armor_cond;
            if (config.armor_cond == "yes"){
                config.armor = 1;
            }
            else    if (config.armor_cond == "no"){
                config.armor = 0;
            }
        }
        else if (line.find("membrane_shape") != -1){
            sin >> config.shape_cond;
            if (config.shape_cond == "sphere"){
                config.shape = 1;
            }
            else    if (config.shape_cond == "cube"){
                config.shape = 2;
            }
            else    if (config.shape_cond == "elipsoid"){
                config.shape = 3;
            }
            else    if (config.shape_cond == "cuboid"){
                config.shape = 4;
            }
        }
        else if (line.find("auroraA") != -1){
            sin >> config.aur_cond;
            if (config.aur_cond == "yes"){
                config.aurora = 1;
            }
            else    if (config.aur_cond == "no"){
                config.aurora = 0;
            }
        }
        else if (line.find("random_place") != -1){
            sin >> config.rand_cond;
            if (config.rand_cond == "yes"){
                config.rand_place = 1;
            }
            else    if (config.rand_cond == "no"){
                config.rand_place = 0;
            }
        }
        else if (line.find("chrom_territory") != -1)
            sin >> config.chrom_ter;
        else if (line.find("ring_radius") != -1){
            sin >> config.ring_rad;
        }
        else if (line.find("ring") != -1){
            sin >> config.ring_cond;
            if (config.ring_cond == "yes"){
                config.ring = 1;
            }
            else    if (config.ring_cond == "no"){
                config.ring = 0;
            }
        }
    }
    ifstream fin_new("config/conf.conf");
    string line_new;
    while (getline(fin_new, line_new)) {
        istringstream sin(line_new.substr(line_new.find(" ") + 1));
        if (line_new.find("kinetochore1_move_x") != -1)
            sin >> config.kinetochore_move_x[0];
        else if (line_new.find("kinetochore1_move_y") != -1)
            sin >> config.kinetochore_move_y[0];
        else if (line_new.find("kinetochore1_move_z") != -1)
            sin >> config.kinetochore_move_z[0];
        else if (line_new.find("x1_angle") != -1)
            sin >> config.x_angle[0];
        else if (line_new.find("y1_angle") != -1)
            sin >> config.y_angle[0];
        else if (line_new.find("z1_angle") != -1)
            sin >> config.z_angle[0];
        else if (line_new.find("kinetochore2_move_x") != -1)
            sin >> config.kinetochore_move_x[1];
        else if (line_new.find("kinetochore2_move_y") != -1)
            sin >> config.kinetochore_move_y[1];
        else if (line_new.find("kinetochore2_move_z") != -1)
            sin >> config.kinetochore_move_z[1];
        else if (line_new.find("x2_angle") != -1)
            sin >> config.x_angle[1];
        else if (line_new.find("y2_angle") != -1)
            sin >> config.y_angle[1];
        else if (line_new.find("z2_angle") != -1)
            sin >> config.z_angle[1];
        else if (line_new.find("kinetochore3_move_x") != -1)
            sin >> config.kinetochore_move_x[2];
        else if (line_new.find("kinetochore3_move_y") != -1)
            sin >> config.kinetochore_move_y[2];
        else if (line_new.find("kinetochore3_move_z") != -1)
            sin >> config.kinetochore_move_z[2];
        else if (line_new.find("x3_angle") != -1)
            sin >> config.x_angle[2];
        else if (line_new.find("y3_angle") != -1)
            sin >> config.y_angle[2];
        else if (line_new.find("z3_angle") != -1)
            sin >> config.z_angle[2];
    }
    ifstream fin_f("config/force.conf");
    string line_f;
    while (getline(fin_f, line_f)) {
        istringstream sin(line_f.substr(line_f.find(" ") + 1));
        if (line_f.find("lang_ts") != -1)
            sin >> config.dt;
        else if (line_f.find("mt_persist_length") != -1)
            sin >> config.lp;
        else if (line_f.find("k_spring_mt") != -1)
            sin >> config.k_spring_mt;
        else if (line_f.find("k_spring_kt") != -1)
            sin >> config.k_spring_kt;
        else if (line_f.find("k_stiff_mem") != -1)
            sin >> config.k_stiff_mem;
        else if (line_f.find("f_push") != -1)
            sin >> config.eps_rep;
        else if (line_f.find("f_pull") != -1)
            sin >> config.eps_att;
        else if (line_f.find("temp") != -1)
            sin >> config.temp;
        else if (line_f.find("viscosity") != -1)
            sin >> config.viscosity;
    }
    ifstream fin_c("config/chemistry.conf");
    string line_c;
    while (getline(fin_c, line_c)) {
        istringstream sin(line_c.substr(line_c.find(" ") + 1));
        if (line_c.find("aurora_diff") != -1)
            sin >> config.aurora_diff;
        else if (line_c.find("k_phos") != -1)
            sin >> config.k_phos;
        else if (line_c.find("k_dephos") != -1)
            sin >> config.k_dephos;
        else if (line_c.find("k_attach") != -1)
            sin >> config.k_attach;
        else if (line_c.find("k_detach_0") != -1)
            sin >> config.k_detach_0;
        else if (line_c.find("k_detach_1") != -1)
            sin >> config.k_detach_1;
        else if (line_c.find("k_detach_2") != -1)
            sin >> config.k_detach_2;
        else if (line_c.find("k_detach_3") != -1)
            sin >> config.k_detach_3;
        else if (line_c.find("k_detach_4") != -1)
            sin >> config.k_detach_4;
        else if (line_c.find("k_detach_5") != -1)
            sin >> config.k_detach_5;
        else if (line_c.find("k_detach_6") != -1)
            sin >> config.k_detach_6;
        else if (line_c.find("k_detach_7") != -1)
            sin >> config.k_detach_7;
        else if (line_c.find("k_growth") != -1)
            sin >> config.k_growth;
        else if (line_c.find("k_shorten") != -1)
            sin >> config.k_shorten;
        else if (line_c.find("k_catastrophe") != -1)
            sin >> config.k_catastrophe;
        else if (line_c.find("k_rescue") != -1)
            sin >> config.k_rescue;
        else if (line_c.find("aurora_concentration") != -1)
            sin >> config.aurora_conc;
        else if (line_c.find("phosphatase_concentration") != -1)
            sin >> config.ph_conc;
        else if (line_c.find("tubulin_concentration") != -1)
            sin >> config.tub_conc;
    }
    ifstream fin_p("config/param.conf");
    string line_p;
    while (getline(fin_p, line_p)) {
        istringstream sin(line_p.substr(line_p.find(" ") + 1));
        if (line_p.find("sv_size") != -1)
            sin >> config.sv_size1;
        else if (line_p.find("ellipse_a") != -1)
            sin >> config.ellipse_a;
        else if (line_p.find("ellipse_b") != -1)
            sin >> config.ellipse_b;
        else if (line_p.find("ellipse_c") != -1)
            sin >> config.ellipse_c;
        else if (line_p.find("h_length") != -1)
            sin >> config.ellipse_a;
        else if (line_p.find("h_width") != -1)
            sin >> config.ellipse_b;
        else if (line_p.find("h_high") != -1)
            sin >> config.ellipse_c;
        else if (line_p.find("sphere_radius") != -1){
            sin >> config.ellipse_a;
        	config.ellipse_b = config.ellipse_a;
        	config.ellipse_c = config.ellipse_a;
        }
        else if (line_p.find("h_size") != -1){
            sin >> config.ellipse_a;
        	config.ellipse_b = config.ellipse_a;
        	config.ellipse_c = config.ellipse_a;
        }
        else if (line_p.find("pole_distance") != -1)
            sin >> config.pole_distance;
        else if (line_p.find("pole_radius") != -1)
            sin >> config.pole_radius;
        else if (line_p.find("pole_count") != -1)
            sin >> config.pole_count;
        else if (line_p.find("mt_num") != -1)
            sin >> config.mt_num;
        else if (line_p.find("max_mt_length") != -1)
            sin >> config.max_mt_length;
        else if (line_p.find("start_length") != -1)
            sin >> config.start_length;
        else if (line_p.find("mt_radius") != -1)
            sin >> config.mt_radius;
        else if (line_p.find("kin_radius") != -1)
            sin >> config.kin_radius;
        else if (line_p.find("beta") != -1)
            sin >> config.beta;
        else if (line_p.find("distance_ndc	  ") != -1)
            sin >> config.distance_ndc;
        else if (line_p.find("maxHarmonicKinPerMonomer") != -1)
            sin >> config.maxHarmonicKinPerMonomer1;
        else if (line_p.find("ndc_length") != -1)
            sin >> config.ndc_length;
        else if (line_p.find("ndc_per_mt") != -1)
            sin >> config.ndc_per_mt;
        else if (line_p.find("corona_S") != -1)
            sin >> config.s_corona;
        else if (line_p.find("corona_XY_ratio") != -1)
            sin >> config.ratio;
        else if (line_p.find("corona_X_N") != -1)
            sin >> config.corona_x_num;
        else if (line_p.find("kin_shape") != -1){
            sin >> config.kin_shape_cond;
        	if (config.kin_shape_cond == "sphere"){
				config.kin_shape = 1;
			}
			else	if (config.kin_shape_cond == "rect"){
				config.kin_shape = 0;
			}
		}
    }
}
