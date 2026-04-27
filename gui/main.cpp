// =============================== main.cpp ===============================
#include <vector>
#include <iostream>
#include <atomic>
#include <thread>
#include <mutex>

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include <GLFW/glfw3.h>
#include "implot.h"

#include "mc_method.h"
#include "biotissue.h"
#include "coordinate.h"
#include "rund_num_generate.h"

// ------------------- Глобальные данные для UI -------------------
std::vector<Layer> userLayers;                       // слои, редактируемые пользователем
std::vector<std::vector<Coordinate>> trajectories;   // траектории фотонов
int photon_count = 1000;
bool ready = false;
bool simulation_running = false;                     // флаг активной симуляции
std::atomic<float> progress{ 0.0f };                   // прогресс 0..1
std::mutex traj_mutex;                               // защита доступа к trajectories

// ------------------- Вспомогательная функция сборки ткани -------------------
Biotissue buildTissueFromUI() {
    Biotissue t;
    for (auto& layer : userLayers) {
        // Layer::l должен быть пересчитан (см. правки в layer.h)
        // Либо вычисляем l на лету, если поле l не используется.
        // В текущей реализации RunOneIterMCM использует cur_layer.l.
        // Поэтому мы обновляем l здесь.
        layer.l = 1.0 / (layer.mu_a + layer.mu_s);
        t.AddLayer(layer);
    }
    return t;
}

// ------------------- Симуляция с прогрессом (запускается в потоке) -------------------
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

// ------------------- Точка входа -------------------
int main() {
    // Начальный слой (как в оригинальной buildTissue)
    userLayers.emplace_back(10.0, 0.1, 0.9, 1.4, 3.0);

    // Инициализация GLFW и окна
    if (!glfwInit())
        return -1;
    GLFWwindow* window = glfwCreateWindow(1200, 800, "MC Simulation", nullptr, nullptr);
    if (!window) {
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    // ImGui / ImPlot
    ImGui::CreateContext();
    ImPlot::CreateContext();
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 130");

    // Главный цикл
    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();

        // Начало кадра ImGui
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        // ---------- Окно управления ----------
        // ---------- Окно управления ----------
        ImGui::SetNextWindowPos(ImVec2(10, 10), ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(400, 500), ImGuiCond_FirstUseEver);
        ImGui::Begin("Control");

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
            else {
                Biotissue tissue = buildTissueFromUI();
                Photon init_photon(0, 0, 0, 0, 0, 1, 1);
                ready = false;
                simulation_running = true;
                progress = 0.0f;
                // захват photon_count по значению
                std::thread sim_thread([tissue, init_photon]() {
                    runSimulationWithProgress(tissue, init_photon, photon_count, trajectories);
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

        // ---------- Окно траекторий ----------
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

        // ---------- Окно траекторий XZ ----------
        ImGui::SetNextWindowPos(ImVec2(420, 620), ImGuiCond_FirstUseEver); // под первым графиком
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

        // Отрисовка OpenGL
        ImGui::Render();
        int w, h;
        glfwGetFramebufferSize(window, &w, &h);
        glViewport(0, 0, w, h);
        glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        glfwSwapBuffers(window);
    }

    // Очистка
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImPlot::DestroyContext();
    ImGui::DestroyContext();
    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}