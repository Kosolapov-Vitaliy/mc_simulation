#ifndef LAYER_H
#define LAYER_H

struct Layer {
    double mu_s;
    double mu_a;
    double l;
    double g;
    double n;
    double z_top;
    double z_bot;
    Layer(double n_mu_s, double n_mu_a, double n_g,
        double n_n, double n_z_top, double n_z_bot) :mu_s(n_mu_s),
        mu_a(n_mu_a) , g(n_g) , n(n_n) , z_top(n_z_top), z_bot(n_z_bot), l(1/(n_mu_a*n_mu_s)){}
};

#endif // !LAYER_H
