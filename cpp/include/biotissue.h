#ifndef BIOTISSUE_H
#define BIOTISSUE_H

#include <layer.h>
#include <photon.h>
#include <vector>

class Biotissue {
private:
    std::vector<Layer> tissue;
public:
    Biotissue() = default;
    Biotissue(const Biotissue& bio) { tissue = bio.tissue; }
    const Biotissue& operator=(const Biotissue& bio) { tissue = bio.tissue; }
    void AddLayer(const Layer& layer) { tissue.push_back(layer); }
    bool PhotonInTissue(const Photon& photon, const double photon_start) const;
    Layer operator[](int i) const { return tissue[i]; }
};

#endif // !BIOTISSUE_H
