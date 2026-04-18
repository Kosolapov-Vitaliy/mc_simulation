#include <mc_method.h>
double SQRTcoef(double dz) {
    return std::sqrt(1 - dz * dz);
}

double CalcDX(double sin_tt, double cos_tt, double sin_fi, double cos_fi,
    double dx, double dy, double dz) {
    return ((sin_tt / SQRTcoef(dz)) * (dx * dz * cos_fi - dy * sin_fi)) + dx * cos_tt;
}

double CalcDY(double sin_tt, double cos_tt, double sin_fi, double cos_fi,
    double dx, double dy, double dz) {
    return ((sin_tt / SQRTcoef(dz)) * (dy * dz * cos_fi + dx * sin_fi)) + dy * cos_tt;
}
double CalcDZ(double sin_tt, double cos_tt, double cos_fi, double dz) {
    return -sin_tt * cos_fi * SQRTcoef(dz) + dz * cos_tt;
}

void RunMCM(const Biotissue& biotissue, Photon& photon, RNGenerate& generator) {
    int layer_idx = 0;
    while (photon.Alive() && photon.weight > photon.weight*0.0001) {
        const Layer cur_layer = biotissue[layer_idx];
        double step = generator.LengthGenerate(cur_layer.l);
        double fi = generator.FiGenerate();
        double cos_tetta = generator.CosTettaGenerate(cur_layer.g);
        double sin_tetta = std::sqrt(1 - cos_tetta * cos_tetta);
        double cos_fi = std::cos(fi);
        double sin_fi = std::sin(fi);
        if (std::abs(photon.dz) > 0.99999)
        {
            double sign_dz = photon.dz / (std::abs(photon.dz));
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
        

    }
}