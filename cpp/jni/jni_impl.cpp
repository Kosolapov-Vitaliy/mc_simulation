#include <jni.h>
#include "com_mcgui_javagui_MCSimulationJNI.h"
#include "biotissue.h"
#include "mc_method.h"
#include <vector>

JNIEXPORT jobjectArray JNICALL Java_com_mcgui_javagui_MCSimulationJNI_runSimulation
(JNIEnv* env, jobject obj, jdoubleArray jlayersParams, jlong numPhotons,
    jdouble sourceX, jdouble sourceY, jdouble sourceZ,
    jdouble sourceDx, jdouble sourceDy, jdouble sourceDz,
    jdouble sourceWeight) {

    // 1. Получаем параметры слоёв из Java-массива
    jsize len = env->GetArrayLength(jlayersParams);
    jdouble* params = env->GetDoubleArrayElements(jlayersParams, nullptr);
    // Параметры: 6 чисел на слой [mua, mus, g, n, z_top, z_bot]
    int numLayers = len / 6;
    Biotissue tissue;
    for (int i = 0; i < numLayers; ++i) {
        double mua = params[i * 6];
        double mus = params[i * 6 + 1];
        double g = params[i * 6 + 2];
        double n = params[i * 6 + 3];
        double z_top = params[i * 6 + 4];
        double z_bot = params[i * 6 + 5];
        // В вашем конструкторе Layer: Layer(mus, mua, g, n, z_top, z_bot)
        Layer layer(mus, mua, g, n, z_top, z_bot);
        tissue.AddLayer(layer);
    }
    env->ReleaseDoubleArrayElements(jlayersParams, params, JNI_ABORT);

    // 2. Создаём начальный фотон
    Photon photon(sourceX, sourceY, sourceZ, sourceDx, sourceDy, sourceDz, sourceWeight);

    // 3. Запускаем симуляцию
    std::vector<std::vector<Coordinate>> trajectories;
    RunSimulation(tissue, photon, (int)numPhotons, trajectories);

    // 4. Преобразуем trajectories в jobjectArray (массив массивов double)
    jclass doubleArrCls = env->FindClass("[D");
    jobjectArray result = env->NewObjectArray((jsize)trajectories.size(), doubleArrCls, nullptr);

    for (size_t i = 0; i < trajectories.size(); ++i) {
        const auto& traj = trajectories[i];
        jdoubleArray pointArray = env->NewDoubleArray((jsize)traj.size() * 3);
        jdouble* raw = env->GetDoubleArrayElements(pointArray, nullptr);
        for (size_t j = 0; j < traj.size(); ++j) {
            raw[j * 3] = traj[j].x;
            raw[j * 3 + 1] = traj[j].y;
            raw[j * 3 + 2] = traj[j].z;
        }
        env->ReleaseDoubleArrayElements(pointArray, raw, 0);
        env->SetObjectArrayElement(result, i, pointArray);
        env->DeleteLocalRef(pointArray);
    }
    return result;
}