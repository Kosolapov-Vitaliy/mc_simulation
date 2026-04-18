#ifndef MC_METHOD_H
#define MC_METHOD_H

#include <photon.h>
#include <biotissue.h>
#include <rund_num_generate.h>

class MCMethod {
private:

public:

};

void RunMCM(const Biotissue& biotissue, Photon& photon, RNGenerate& generator);
int CrossBoundary(const Biotissue& biotissue, Photon& photon, RNGenerate& generator,
    double remaining_step, int layer_idx, double& new_step, bool& is_out);
#endif // !MC_METHOD_H
