#ifndef RUND_NUM_GENERATE_H
#define RUND_NUM_GENERATE_H

#include <iostream>
#include <cmath>
#include <random>


class RNGenerate {
private:
    std::random_device rd;
public:
    double LengthGenerate(double a_len);
    double FiGenerate();
    double CosTettaGenerate(double g);
    double KsiGenerate();
};

#endif // !RUND_NUM_GENERATE_H
