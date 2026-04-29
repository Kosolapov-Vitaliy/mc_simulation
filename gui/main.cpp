#define _USE_MATH_DEFINES
#include <vector>
#include <iostream>
#include <atomic>
#include <thread>
#include <mutex>
#include <cmath>
#include <future>
#include <utility>
#include <iterator>


#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include <GLFW/glfw3.h>
#include "implot.h"

#include "mc_method.h"
#include "biotissue.h"
#include "coordinate.h"
#include "rund_num_generate.h"

std::vector<Layer> userLayers;
std::vector<std::vector<Coordinate>> trajectories;   
std::vector<double> detected;
int photon_count = 1000;
bool ready = false;
bool simulation_running = false;
std::atomic<float> progress{ 0.0f };
std::mutex traj_mutex;

Biotissue buildTissueFromUI() {
    Biotissue t;
    for (auto& layer : userLayers) {
        layer.l = 1.0 / (layer.mu_a + layer.mu_s);
        t.AddLayer(layer);
    }
    return t;
}

void runSimulationWithProgress(const Biotissue& tissue, const Photon& init_photon,
    int num_photons, std::vector<std::vector<Coordinate>>& out_trajectories) {
    RNGenerate generator;
    out_trajectories.clear();
    out_trajectories.reserve(num_photons);

    for (int i = 0; i < num_photons; ++i) {
        std::vector<Coordinate> cur_path;
        Photon cur_photon = init_photon;
        RunOneIterMCM(tissue, cur_photon, generator, cur_path);
        {
            std::lock_guard<std::mutex> lock(traj_mutex);
            out_trajectories.push_back(std::move(cur_path));
        }
        progress = static_cast<float>(i + 1) / num_photons;
    }
}

void runSimulationDetectedWithProgress(const Biotissue& tissue, const Photon& init_photon,
    int num_photons, std::vector<std::vector<Coordinate>>& out_trajectories, std::vector<double>& detected) {
    const double max_distance = 10;
    const double step = 0.1;
    const int detector_count = static_cast<int>(max_distance / step);
    unsigned int threads = std::thread::hardware_concurrency();
    std::vector<std::future<std::pair<std::vector<std::vector<Coordinate>>, std::vector<double>>>> futures;

    int photons_per_thread = num_photons / threads;
    int remainder = num_photons % threads;
    int start = 0;
    std::atomic<int> processed{ 0 };

    for (int t = 0; t < threads; t++) {
        int photons_count = photons_per_thread + (t < remainder ? 1 : 0);
        futures.push_back(std::async(std::launch::async, [&, start, photons_count]() {
            RNGenerate local_gen;
            std::vector<std::vector<Coordinate>> local_traj;
            local_traj.reserve(photons_count);
            std::vector<double> local_detected(detector_count, 0.0);
            for (int i = 0; i < photons_count; i++) {
                std::vector<Coordinate> path;
                Photon photon = init_photon;
                double start_x = photon.x;
                double start_y = photon.y;
                double start_z = photon.z;
                RunOneIterMCM(tissue, photon, local_gen, path);
                local_traj.push_back(std::move(path));
                double last_x = photon.x;
                double last_y = photon.y;
                double last_z = photon.z;
                if (last_z == start_z) {
                    double dx = last_x - start_x, dy = last_y - start_y;
                    double r2 = dx * dx + dy * dy;
                    for (int j = 0; j < detector_count+1; j++) {
                        double rMin = j * step;
                        double rMax = (j + 1) * step;
                        if (r2 > rMin * rMin && r2 <= rMax * rMax) {
                            double area = M_PI * (rMax * rMax - rMin * rMin);
                            local_detected[j] += 1.0 / area;
                            break;
                        }
                    }
                }
                processed.fetch_add(1);
                progress = static_cast<float>(processed.load()) / num_photons;
            }
            
            return std::make_pair(std::move(local_traj), std::move(local_detected));
            }));
        start += photons_count;
    }
    out_trajectories.clear();
    detected.assign(detector_count, 0.0);
    for (auto& f : futures) {
        auto [traj, det] = f.get();
        out_trajectories.insert(out_trajectories.end(),
            std::make_move_iterator(traj.begin()),
            std::make_move_iterator(traj.end()));
        for (int j = 0; j < det.size(); j++) {
            detected[j] += det[j];
        }
    }
}

int main() {
    userLayers.emplace_back(10.0, 0.1, 0.9, 1.4, 3.0);

    if (!glfwInit())
        return -1;
    GLFWwindow* window = glfwCreateWindow(1200, 800, "MC Simulation", nullptr, nullptr);
    if (!window) {
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    ImGui::CreateContext();
    ImPlot::CreateContext();
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 130");

    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();

        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        ImGui::SetNextWindowPos(ImVec2(10, 10), ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(400, 500), ImGuiCond_FirstUseEver);
        ImGui::Begin("Control");
        bool correct = true;

        ImGui::InputInt("Photons", &photon_count);
        if (photon_count < 1) photon_count = 1;

        ImGui::SeparatorText("Tissue Layers");

        if (ImGui::Button("+ Add Layer")) {
            userLayers.emplace_back(10.0, 0.1, 0.9, 1.4, 3.0);
        }
        ImGui::SameLine();
        if (ImGui::Button("Clear All Layers") && !simulation_running) {
            userLayers.clear();
        }

        for (int i = 0; i < (int)userLayers.size(); ++i) {
            ImGui::PushID(i);
            ImGui::Text("Layer %d", i);
            auto& lay = userLayers[i];
            ImGui::InputDouble("mu_s", &lay.mu_s, 0.5, 1.0);
            ImGui::InputDouble("mu_a", &lay.mu_a, 0.05, 0.1);
            ImGui::InputDouble("g", &lay.g, 0.01, 0.1);
            ImGui::InputDouble("n", &lay.n, 0.02, 0.1);
            ImGui::InputDouble("thickness", &lay.thickness, 0.5, 1.0);
            if (lay.mu_s == 0 && lay.mu_a == 0) {
                correct = false;
            }
            if (ImGui::Button("Remove")) {
                userLayers.erase(userLayers.begin() + i);
                ImGui::PopID();
                break;
            }
            ImGui::Separator();
            ImGui::PopID();
        }

        if (ImGui::Button("Run simulation") && !simulation_running) {
            if (userLayers.empty()) {
                ImGui::TextColored(ImVec4(1, 0, 0, 1), "Error: no layers!");
            }
            else if (!correct) {
                ImGui::TextColored(ImVec4(1, 0, 0, 1), "Error: one of layer is incorrect!");
            }
            else {
                Biotissue tissue = buildTissueFromUI();
                Photon init_photon(0, 0, 0, 0, 0, 1, 1);
                ready = false;
                simulation_running = true;
                progress = 0.0f;
                std::thread sim_thread([tissue, init_photon]() {
                    runSimulationDetectedWithProgress(tissue, init_photon, photon_count, trajectories, detected);
                    ready = true;
                    simulation_running = false;
                    });
                sim_thread.detach();
            }
        }

        if (simulation_running) {
            ImGui::ProgressBar(progress.load(), ImVec2(-1, 0));
            ImGui::SameLine();
            ImGui::Text("%.1f%%", progress.load() * 100);
        }
        else if (ready) {
            ImGui::TextColored(ImVec4(0, 1, 0, 1), "Simulation finished. Paths: %zu", trajectories.size());
        }

        ImGui::End();

        ImGui::SetNextWindowPos(ImVec2(420, 10), ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(750, 600), ImGuiCond_FirstUseEver);
        ImGui::Begin("Trajectories XY");

        if (ready && ImPlot::BeginPlot("Monte Carlo paths (XY)", ImVec2(-1, -1))) {
            std::lock_guard<std::mutex> lock(traj_mutex);
            std::vector<double> end_x, end_y;
            for (const auto& path : trajectories) {
                if (path.empty()) continue;
                std::vector<double> x, y;
                x.reserve(path.size());
                y.reserve(path.size());
                for (const auto& p : path) {
                    x.push_back(p.x);
                    y.push_back(p.y);
                }
                ImPlot::PlotLine("path", x.data(), y.data(), (int)x.size());
                const auto& last = path.back();
                end_x.push_back(last.x);
                end_y.push_back(last.y);
            }
            if (!end_x.empty()) {
                ImPlot::PlotScatter("end points", end_x.data(), end_y.data(), (int)end_x.size());
            }            
            ImPlot::EndPlot();
        }
        else if (!ready && !simulation_running && trajectories.empty()) {
            ImGui::Text("Press 'Run simulation' to start.");
        }

        ImGui::End();

        ImGui::SetNextWindowPos(ImVec2(420, 620), ImGuiCond_FirstUseEver); 
        ImGui::SetNextWindowSize(ImVec2(750, 300), ImGuiCond_FirstUseEver);
        ImGui::Begin("Trajectories XZ");

        if (ready && ImPlot::BeginPlot("Monte Carlo paths (XZ)", ImVec2(-1, -1))) {
            std::lock_guard<std::mutex> lock(traj_mutex);
            std::vector<double> end_x, end_z;
            for (const auto& path : trajectories) {
                if (path.empty()) continue;
                std::vector<double> x, z;
                x.reserve(path.size());
                z.reserve(path.size());
                for (const auto& p : path) {
                    x.push_back(p.x);
                    z.push_back(p.z);
                }
                ImPlot::PlotLine("path", x.data(), z.data(), (int)x.size());
                const auto& last = path.back();
                end_x.push_back(last.x);
                end_z.push_back(last.z);
            }
            if (!end_x.empty()) {
                ImPlot::PlotScatter("end points", end_x.data(), end_z.data(), (int)end_x.size());
            }
            double prevThickness = 0;
            ImPlotSpec spec;
            spec.LineColor = ImVec4(1.0f, 0.0f, 0.0f, 1.0f);
            spec.Flags = ImPlotInfLinesFlags_Horizontal;
            for (const auto& layer : userLayers) {
                prevThickness += layer.thickness;
                ImPlot::PlotInfLines("layer_border", &prevThickness, 1, spec);
            }
            ImPlot::EndPlot();
        }
        else if (!ready && !simulation_running && trajectories.empty()) {
            ImGui::Text("Press 'Run simulation' to start.");
        }

        ImGui::End();

        ImGui::SetNextWindowPos(ImVec2(420, 620), ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(750, 300), ImGuiCond_FirstUseEver);
        ImGui::Begin("Denisty plot");
        double max_distance = 10;
        double step = 0.1;
        if (ready && ImPlot::BeginPlot("Monte Carlo Denisty plot", ImVec2(-1, -1))) {
            std::lock_guard<std::mutex> lock(traj_mutex);
            std::vector<double> x, y;
            x.reserve(detected.size());
            y.reserve(detected.size());
            for (int j = 0; j < (int)(max_distance / step);j++) {
                x.push_back(j*step);
                y.push_back(detected[j]);
                ImPlot::PlotLine("path", x.data(), y.data(), (int)x.size());
            }
            ImPlot::EndPlot();
        }
        else if (!ready && !simulation_running && trajectories.empty()) {
            ImGui::Text("Press 'Run simulation' to start.");
        }

        ImGui::End();

        ImGui::SetNextWindowPos(ImVec2(420, 620), ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(750, 300), ImGuiCond_FirstUseEver);
        ImGui::Begin("Denisty plot (logariphm)");
        if (ready && ImPlot::BeginPlot("Monte Carlo Denisty plot (logariphm)", ImVec2(-1, -1))) {
            std::lock_guard<std::mutex> lock(traj_mutex);
            std::vector<double> x, y;
            x.reserve(detected.size());
            y.reserve(detected.size());
            for (int j = 0; j < (int)(max_distance / step); j++) {
                x.push_back(j * step);
                y.push_back(std::log10(detected[j]));
                ImPlot::PlotLine("path", x.data(), y.data(), (int)x.size());
            }
            ImPlot::EndPlot();
        }
        else if (!ready && !simulation_running && trajectories.empty()) {
            ImGui::Text("Press 'Run simulation' to start.");
        }

        ImGui::End();

        ImGui::Render();
        int w, h;
        glfwGetFramebufferSize(window, &w, &h);
        glViewport(0, 0, w, h);
        glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        glfwSwapBuffers(window);
    }

    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImPlot::DestroyContext();
    ImGui::DestroyContext();
    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}