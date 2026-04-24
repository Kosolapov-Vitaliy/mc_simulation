package com.mcgui.javagui;

public class MCSimulationJNI {
    static{
        System.loadLibrary("mc_simulation_jni");
    }
    public native Object[][] runSimulation(double[] layersParams, long numPhotons,
                                           double sourceX, double sourceY, double sourceZ,
                                           double sourceDx, double sourceDy, double sourceDz,
                                           double sourceWeight);
}
