#include <rund_num_generate.h>
#include <numbers>

RNGenerate::RNGenerate(){
    std::random_device rd;
    gen = std::mt19937(rd());
}

double RNGenerate::LengthGenerate(double a_len) {
    std::uniform_real_distribution<> dist(0, 1);
    double ksi = dist(gen);
    double res_len = -std::log(1-ksi)*a_len;
    return res_len;
}

double RNGenerate::FiGenerate() {
    std::uniform_real_distribution<> dist(0, 2*std::numbers::pi);
    double fi = dist(gen);
    return fi;
}

double RNGenerate::CosTettaGenerate(double g) {
    std::uniform_real_distribution<> dist(0, 1);
    double ksi = dist(gen);
    double cos_tetta = 0;
    if (g > 0)
    {
        double g_p2 = g * g;
        double tFrag = (1 - g_p2) / (1 - g + 2 * g * ksi);
        cos_tetta = (1 / (2 * g)) * (1 + g_p2 - (tFrag * tFrag));
    }
    else if (g == 0)
    {
        cos_tetta = 2 * ksi - 1;
    }
    return cos_tetta;
}

double RNGenerate::KsiGenerate() {
    std::uniform_real_distribution<> dist(0, 1);
    double ksi = dist(gen);
    return ksi;
}