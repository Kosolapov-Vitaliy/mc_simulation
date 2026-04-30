#include <mc_method.h>
#include <iostream>


double CalcFRCoef(double dz, double n_i, double n_t) {
    double alpha_i = std::acos(std::fabs(dz));
    if (alpha_i == 0) {
        return ((n_t - n_i) / (n_t + n_i)) * ((n_t - n_i) / (n_t + n_i));
    }
    else if (alpha_i > 0 && (n_i*std::sin(alpha_i))/n_t < 1) {

        double alpha_t = std::asin((n_i / n_t) * std::sin(alpha_i));
        double sum_sin =std::sin(alpha_i + alpha_t);
        double asg_sin = std::sin(alpha_i - alpha_t);
        double sum_tg = std::tan(alpha_i + alpha_t);
        double asg_tg = std::tan(alpha_i - alpha_t);
        return (1 / 2) * (((asg_sin * asg_sin) / (sum_sin * sum_sin)) + ((asg_tg * asg_tg) / (sum_tg * sum_tg)));
    }
    else {                 //if (alpha_i > std::asin(n_i / n_t))
        return 1;
    }
}

void Refract(double& dx, double& dy, double& dz, double n_i, double n_t) {
    double cos_alpha_i = std::fabs(dz);
    double sin_alpha_i = std::sqrt(1.0 - cos_alpha_i * cos_alpha_i);
    double sin_alpha_t = (n_i / n_t) * sin_alpha_i;
    double cos_alpha_t = std::sqrt(1.0 - sin_alpha_t * sin_alpha_t);
    double sign_dz = dz/std::fabs(dz);
    double scale = sin_alpha_t / sin_alpha_i;
    dx *= scale;
    dy *= scale;
    dz = sign_dz * cos_alpha_t;
    double norm = std::sqrt(dx * dx + dy * dy + dz * dz);
    if (norm > 1e-12) {
        dx /= norm;
        dy /= norm;
        dz /= norm;
    }
}

void RunOneIterMCM(const Biotissue& biotissue, Photon& photon, RNGenerate& generator,
    std::vector<Coordinate>& pathway, double n_external, double n_depth) {
    int layer_idx = 0;
    double z_min = photon.z;
    double photon_start = photon.z;
    pathway.push_back(Coordinate(photon.x, photon.y, photon.z));
    double start_weight = photon.weight;
    int i = 1;
    bool in_tissue = true;
    while (photon.Alive(start_weight)&&in_tissue) {
        Layer cur_layer = biotissue[layer_idx];
        double step = generator.LengthGenerate(cur_layer.l);
        double fi = generator.FiGenerate();
        double cos_tetta = generator.CosTettaGenerate(cur_layer.g);
        double sin_tetta = std::sqrt(1 - cos_tetta * cos_tetta);
        double cos_fi = std::cos(fi);
        double sin_fi = std::sin(fi);
        double z_max = cur_layer.thickness + z_min;
        double new_dx, new_dy, new_dz;
        if (std::fabs(photon.dz) > 0.99999)
        {
            double sign_dz = photon.dz / (std::fabs(photon.dz));
            new_dx = cos_fi*sin_tetta;
            new_dy = sin_fi*sin_tetta;
            new_dz = sign_dz*cos_tetta;
        }
        else {
            new_dx = CalcDX(sin_tetta, cos_tetta, sin_fi, cos_fi, photon.dx, photon.dy, photon.dz);
            new_dy = CalcDY(sin_tetta, cos_tetta, sin_fi, cos_fi, photon.dx, photon.dy, photon.dz);
            new_dz = CalcDZ(sin_tetta, cos_tetta, cos_fi, photon.dz);
        } 
        double norm = std::sqrt((new_dx * new_dx) + new_dy * new_dy + new_dz * new_dz);
        new_dx /= norm; new_dy /= norm; new_dz /= norm;
        photon.dx = new_dx; photon.dy = new_dy; photon.dz = new_dz;
        photon.x = photon.x + step * photon.dx;
        photon.y = photon.y + step * photon.dy;
        photon.z = photon.z + step * photon.dz;
        photon.weight -= (photon.weight * cur_layer.mu_a * cur_layer.l);
        if (photon.z <= z_min || photon.z >= z_max) {
            bool border = true;
            while (border) {
                double nec_z = 0.0;
                double frcoef =0.0;
                bool flag = (photon.z <= z_min);
                bool out = false;
                cur_layer = biotissue[layer_idx];                
                if (flag) {
                    nec_z = z_min;
                    if ((layer_idx > 0)){
                        frcoef = CalcFRCoef(photon.dz, cur_layer.n, biotissue[layer_idx - 1].n);
                    }
                    else {
                        frcoef = CalcFRCoef(photon.dz, cur_layer.n, n_external); //Граница со слоем откуда идут фотоны, будем считать, что воздух
                        out = true;
                    }
                }
                else  {
                    nec_z = z_max;
                    if ((layer_idx + 1 < biotissue.GetLayerCount())) {
                        frcoef = CalcFRCoef(photon.dz, cur_layer.n, biotissue[layer_idx + 1].n);
                    }
                    else {
                        frcoef = CalcFRCoef(photon.dz, cur_layer.n, n_depth); //Граница со слоем который идёт дальше вглубь, будем считать, что тоже воздух
                        out = true;
                    }
                }
                double back_step = (photon.z - nec_z) / photon.dz;
                photon.z = nec_z;
                photon.x = photon.x - (photon.dx * back_step);
                photon.y = photon.y - (photon.dy * back_step);
                double norm = std::sqrt((photon.dx * photon.dx) + (photon.dy * photon.dy) + (photon.dz * photon.dz));
                photon.dz = photon.dz / norm;
                photon.dx = photon.dx / norm;
                photon.dy = photon.dy / norm;
                double ksi = generator.KsiGenerate();
                if (ksi <= frcoef) {
                    photon.dz = -photon.dz;
                    photon.dx = photon.dx;
                    photon.dy = photon.dy;
                    photon.x = photon.x + back_step * photon.dx;
                    photon.y = photon.y + back_step * photon.dy;
                    photon.z = photon.z + back_step * photon.dz;
                }
                else if (flag&&(!(out))){
                    double prev_l = cur_layer.l;
                    double prev_n = cur_layer.n;
                    layer_idx--;
                    z_min = z_min - biotissue[layer_idx].thickness;
                    double new_step = (back_step * biotissue[layer_idx].l) / prev_l;
                    Refract(photon.dx, photon.dy, photon.dz, prev_n, biotissue[layer_idx].n);
                    photon.x = photon.x + new_step * photon.dx;
                    photon.y = photon.y + new_step * photon.dy;
                    photon.z = photon.z + new_step * photon.dz;
                }
                else if ((!flag) && (!(out))) {
                    double prev_l = cur_layer.l;
                    double prev_n = cur_layer.n;
                    z_min = z_max;
                    layer_idx++;
                    double new_step = (back_step * biotissue[layer_idx].l) / prev_l;
                    Refract(photon.dx, photon.dy, photon.dz, prev_n, biotissue[layer_idx].n);
                    photon.x = photon.x + back_step * photon.dx;
                    photon.y = photon.y + back_step * photon.dy;
                    photon.z = photon.z + back_step * photon.dz;
                }
                else {
                    in_tissue = false;
                    border = false;
                }
                z_max = cur_layer.thickness + z_min;
                if (!(photon.z <= z_min || photon.z >= z_max)) {
                    border = false;
                }

            }            
        }
        pathway.push_back(Coordinate(photon.x, photon.y, photon.z));
        i++;
    }
}

void RunSimulation(const Biotissue& biotissue, const Photon& photon, int num_photons,
    std::vector<std::vector<Coordinate>>& trajectorys, double n_external, double n_depth) {
    RNGenerate generator = RNGenerate();
    for (int i = 0; i < num_photons; i++) {
        std::vector<Coordinate> cur_pathway;
        Photon cur_photon = photon;
        RunOneIterMCM(biotissue, cur_photon, generator, cur_pathway, n_external, n_depth);
        trajectorys.push_back(cur_pathway);
    }
}

