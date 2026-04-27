#include <biotissue.h>

bool Biotissue::PhotonInTissue(const Photon& photon, const double photon_start) const {
    double z_min = photon_start;
    for (int i = 0; i < tissue.size(); i++) {
        if (photon.z>=z_min&&photon.z<=(tissue[i].thickness+z_min)) {
            return true;
        }
        z_min += tissue[i].thickness;
    }
    return false;
}