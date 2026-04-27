#ifndef PHOTON_H
#define PHOTON_H

struct Photon {
    double x;
    double y;
    double z;
    double dx;
    double dy;
    double dz;
    double weight;
    Photon() :x(0.0), y(0.0), z(0.0), 
        dx(0.0), dy(0.0), dz(0.0), weight(1.0) {}
    Photon(double n_weight):x(0.0), y(0.0), z(0.0),
        dx(0.0), dy(0.0), dz(0.0), weight(n_weight) {}
    Photon(double n_x, double n_y, double n_z,
        double n_dx, double n_dy, double n_dz,
        double n_weight):
        x(n_x), y(n_y), z(n_z),
        dx(n_dx), dy(n_dy), dz(n_dz),
        weight(n_weight) {}
    Photon(const Photon& p) :x(p.x), y(p.y), z(p.z),
        dx(p.dx), dy(p.dy), dz(p.dz), weight(p.weight) {}
    bool Alive(double start_weight) { return weight >= start_weight * 0.0001; }
};

#endif // !PHOTON_H
