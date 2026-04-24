#include <mc_method.h>


double CalcFRCoef(double dz, double n_i, double n_t) {
    double alpha_i = std::acos(std::fabs(dz));
    double alpha_t = std::asin((n_i / n_t) * std::sin(alpha_i));
    if (alpha_i == 0) {
        return ((n_t - n_i) / (n_t + n_i)) * ((n_t - n_i) / (n_t + n_i));
    }
    else if (alpha_i > 0 && alpha_i <= std::acos(n_i / n_t)) {
        double sum_sin = alpha_i + alpha_t;
        double asg_sin = alpha_i - alpha_t;
        double sum_tg = alpha_i + alpha_t;
        double asg_tg = alpha_i - alpha_t;
        return (1 / 2) * (((asg_sin * asg_sin) / (sum_sin * sum_sin)) + ((asg_tg * asg_tg) / (sum_tg * sum_tg)));
    }
    else if (alpha_i > std::acos(n_i / n_t)) {
        return 1;
    }
}

void RunOneIterMCM(const Biotissue& biotissue, Photon& photon, RNGenerate& generator,
    std::vector<Coordinate>& pathway) {
    int layer_idx = 0;
    pathway.push_back(Coordinate(photon.x, photon.y, photon.z));
    while (photon.Alive()) {
        const Layer cur_layer = biotissue[layer_idx];
        double step = generator.LengthGenerate(cur_layer.l);
        double fi = generator.FiGenerate();
        double cos_tetta = generator.CosTettaGenerate(cur_layer.g);
        double sin_tetta = std::sqrt(1 - cos_tetta * cos_tetta);
        double cos_fi = std::cos(fi);
        double sin_fi = std::sin(fi);
        if (std::fabs(photon.dz) > 0.99999)
        {
            double sign_dz = photon.dz / (std::fabs(photon.dz));
            photon.dx = cos_fi*sin_tetta;
            photon.dy = sin_fi*sin_tetta;
            photon.dz = sign_dz*cos_tetta;
        }
        else {
            photon.dx = CalcDX(sin_tetta, cos_tetta, sin_fi, cos_fi, photon.dx, photon.dy, photon.dz);
            photon.dy = CalcDY(sin_tetta, cos_tetta, sin_fi, cos_fi, photon.dx, photon.dy, photon.dz);
            photon.dz = CalcDZ(sin_tetta, cos_tetta, cos_fi, photon.dz);
        }
        photon.x = photon.x + step * photon.dx;
        photon.y = photon.y + step * photon.dy;
        photon.z = photon.z + step * photon.dz;
        photon.weight = photon.weight * cur_layer.mu_a * cur_layer.l;
        if (photon.z <= cur_layer.z_bot || photon.z >= cur_layer.z_top) {
            double nec_z =0.0;
            bool flag = false;
            bool up = false;
            if (biotissue.PhotonInTissue(photon)) {
                flag = true;
            }
            if (photon.z <= cur_layer.z_bot) {
                nec_z = cur_layer.z_bot;
            } else if (photon.z >= cur_layer.z_top) {
                nec_z = cur_layer.z_top;
                up = true;
            }
            double back_step = (photon.z - nec_z) / photon.dz;
            photon.z = nec_z;
            photon.x = photon.x - (photon.dx * back_step);
            photon.y = photon.y - (photon.dy * back_step);
            double frcoef = 0.0;
            double ksi = generator.KsiGenerate();
            if (flag) {
                if (up) {
                    frcoef = CalcFRCoef(photon.dz, cur_layer.n, biotissue[layer_idx + 1].n);
                    if (ksi <= frcoef) {
                        photon.dz = -photon.dz;
                    }
                    else {
                        layer_idx++;
                    }
                }
                else {
                    frcoef = CalcFRCoef(photon.dz, cur_layer.n, biotissue[layer_idx - 1].n);
                    if (ksi <= frcoef) {
                        photon.dz = -photon.dz;
                    }
                    else {
                        layer_idx--;
                    }
                }
            }
            else {
                frcoef = CalcFRCoef(photon.dz, cur_layer.n, 1.0);
                if (ksi <= frcoef) {
                    photon.dz = -photon.dz;
                }
                else {
                    photon.weight = 0;
                }
            }
        }
        pathway.push_back(Coordinate(photon.x, photon.y, photon.z));
    }
}

void RunSimulation(const Biotissue& biotissue, const Photon& photon, int num_photons,
    std::vector<std::vector<Coordinate>>& trajectorys) {
    RNGenerate generator = RNGenerate();
    for (int i = 0; i < num_photons; i++) {
        std::vector<Coordinate> cur_pathway;
        Photon cur_photon = photon;
        RunOneIterMCM(biotissue, cur_photon, generator, cur_pathway);
        trajectorys.push_back(cur_pathway);
    }
}