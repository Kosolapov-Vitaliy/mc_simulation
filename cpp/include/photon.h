#ifndef PHOTON_H
#define PHOTON_H

class Photon {
private:
    double x;
    double y;
    double z;
    double weight;
public:
    Photon() :x(0.0), y(0.0), z(0.0), weight(1.0) {}
    Photon(double n_weight):x(0.0), y(0.0), z(0.0), weight(n_weight) {}
    Photon(double n_x, double n_y, double n_z, double n_weight):
        x(n_x), y(n_y), z(n_z), weight(n_weight) {}
    Photon(const Photon& p) :x(p.x), y(p.y), z(p.z), weight(p.weight) {}
    double GetX() { return x; }
    double GetY() { return y; }
    double GetZ() { return z; }
    double GetWeight() { return weight; }
    bool Alive() { return weight != 0; }
    void SetX(double n_x) { x = n_x; }
    void SetX(double n_y) { x = n_y; }
    void SetX(double n_z) { x = n_z; }
    void SetX(double n_weight) { x = n_weight; }
};

#endif // !PHOTON_H
