#include <biotissue.h>

bool Biotissue::PhotonInTissue(const Photon& photon) {
    for (int i = 0; i < tissue.size(); i++) {
        if (photon.GetZ()>=tissue[i].z_bot&&photon.GetZ()<=tissue[i].z_top) {
            return true;
        }
    }
    return false;
}