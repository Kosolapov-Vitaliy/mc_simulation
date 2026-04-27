#ifndef LAYER_H
#define LAYER_H

struct Layer {
    double mu_s;
    double mu_a;
    double l;
    double g;
    double n;
    double thickness;
    Layer(double n_mu_s, double n_mu_a, double n_g,
        double n_n, double n_thickness) :mu_s(n_mu_s),
        mu_a(n_mu_a) , g(n_g) , n(n_n) , thickness(n_thickness), l(1.0/(n_mu_a + n_mu_s)){}
};

#endif // !LAYER_H
