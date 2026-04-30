#ifndef MC_METHOD_H
#define MC_METHOD_H

#include <photon.h>
#include <coordinate.h>
#include <biotissue.h>
#include <rund_num_generate.h>


inline double SQRTcoef(double dz) {
    return std::sqrt(1 - dz * dz);
}

inline double CalcDX(double sin_tt, double cos_tt, double sin_fi, double cos_fi,
    double dx, double dy, double dz) {
    return ((sin_tt / SQRTcoef(dz)) * (dx * dz * cos_fi - dy * sin_fi)) + dx * cos_tt;
}

inline double CalcDY(double sin_tt, double cos_tt, double sin_fi, double cos_fi,
    double dx, double dy, double dz) {
    return ((sin_tt / SQRTcoef(dz)) * (dy * dz * cos_fi + dx * sin_fi)) + dy * cos_tt;
}
inline double CalcDZ(double sin_tt, double cos_tt, double cos_fi, double dz) {
    return -sin_tt * cos_fi * SQRTcoef(dz) + dz * cos_tt;
}
double CalcFRCoef(double dz, double n_i, double n_t);
void RunOneIterMCM(const Biotissue& biotissue, Photon& photon, RNGenerate& generator,
    std::vector<Coordinate>& pathway, double n_external, double n_depth);
void RunSimulation(const Biotissue& biotissue,const Photon& photon, int num_photons,
    std::vector<std::vector<Coordinate>>& trajectorys, double n_external, double n_depth);
#endif // !MC_METHOD_H
