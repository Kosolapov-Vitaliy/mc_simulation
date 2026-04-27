#include <biotissue.h>

bool Biotissue::PhotonInTissue(const Photon& photon) const {
    for (int i = 0; i < tissue.size(); i++) {
        if (photon.z>=tissue[i].z_bot&&photon.z<=tissue[i].z_top) {
            return true;
        }
    }
    return false;
}