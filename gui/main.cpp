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
    #include <string>


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
    std::vector<std::pair<double,double>> max_deep_photons;
    std::vector<double> detected;
    int photon_count = 1000;
    bool ready = false;
    bool simulation_running = false;
    std::atomic<float> progress{ 0.0f };
    std::mutex traj_mutex;
    double n_external;
    double n_depth;
    bool saveDR = false;
    bool table_need_recalculation = false;

    double max_deep = 0.0;
    double max_distance = 10;
    double step = 0.1;
    std::vector<std::vector<double>> savedDenisty;
    int defer = 0;
    int size_deptdist_hm = 1000;

    static std::vector<std::vector<int>> table;


    Biotissue buildTissueFromUI(double& max_deep) {
        Biotissue t;
        max_deep = 0.0;
        for (auto& layer : userLayers) {
            layer.l = 1.0 / (layer.mu_a + layer.mu_s);
            max_deep += layer.thickness;
            t.AddLayer(layer);
        }
        return t;
    }


    void runSimulationDetectedWithProgress(const Biotissue& tissue, const Photon& init_photon,
        int num_photons, std::vector<std::vector<Coordinate>>& out_trajectories, std::vector<std::pair<double, double>>& photons_max_deep, std::vector<double>& detected,
        double n_external, double n_depth) {
        const double max_distance = 10;
        const double step = 0.1;
        const int detector_count = static_cast<int>(max_distance / step);
        unsigned int threads = std::thread::hardware_concurrency();
        std::vector<std::future<std::tuple<std::vector<std::vector<Coordinate>>, std::vector<double>, std::vector<std::pair<double, double>>>>> futures;

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
                std::vector<std::pair<double, double>> local_max_deep;
                local_max_deep.reserve(photons_count);
                std::vector<double> local_detected(detector_count, 0.0);
                for (int i = 0; i < photons_count; i++) {
                    std::vector<Coordinate> path;
                    double photon_max_deep;
                    Photon photon = init_photon;
                    double start_x = photon.x;
                    double start_y = photon.y;
                    double start_z = photon.z;
                    RunOneIterMCM(tissue, photon, local_gen, path, photon_max_deep, n_external, n_depth);
                    local_traj.push_back(std::move(path));                
                    double last_x = photon.x;
                    double last_y = photon.y;
                    double last_z = photon.z;
                    if (last_z == start_z) {
                        local_max_deep.push_back({ std::move(photon_max_deep), std::sqrt((last_x * last_x) + (last_y * last_y)) });
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
            
                return std::make_tuple(std::move(local_traj), std::move(local_detected), std::move(local_max_deep));
                }));
            start += photons_count;
        }
        out_trajectories.clear();
        photons_max_deep.clear();
        detected.assign(detector_count, 0.0);
        for (auto& f : futures) {
            auto [traj, det, max_d] = f.get();
            out_trajectories.insert(out_trajectories.end(),
                std::make_move_iterator(traj.begin()),
                std::make_move_iterator(traj.end()));
            for (int j = 0; j < det.size(); j++) {
                detected[j] += det[j];
            }
            for (int j = 0; j < max_d.size(); j++) {
                photons_max_deep.push_back(max_d[j]);
            }
        }
    }

    int main() {
        userLayers.emplace_back(10.0, 0.1, 0.9, 1.4, 3.0);
        table.assign(100, std::vector<int>(100, 0));

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
        ImPlotStyle& style = ImPlot::GetStyle();
        style.Colors[ImPlotCol_PlotBg] = ImVec4(1.0f, 1.0f, 1.0f, 1.0f);
        style.Colors[ImPlotCol_AxisGrid] = ImVec4(0.5f, 0.0f, 0.0f, 0.5f);
        style.Colors[ImPlotCol_AxisTick] = ImVec4(0.5f, 0.0f, 0.0f, 0.5f);        
        ImGui_ImplGlfw_InitForOpenGL(window, true);
        ImGui_ImplOpenGL3_Init("#version 130");
        ImFont* pMyBigFont = ImGui::GetIO().Fonts->AddFontFromFileTTF("C:\\Windows\\Fonts\\Arial.ttf", 24.0f);
        if (!pMyBigFont) {
            printf("Не удалось загрузить шрифт!\n");
            pMyBigFont = ImGui::GetIO().Fonts->AddFontDefault();
        }

        while (!glfwWindowShouldClose(window)) {
            glfwPollEvents();

            ImGui_ImplOpenGL3_NewFrame();
            ImGui_ImplGlfw_NewFrame();
            ImGui::NewFrame();
            ImGuiViewport* viewport = ImGui::GetMainViewport();

            ImGui::SetNextWindowPos(viewport->Pos);
            ImGui::SetNextWindowSize(viewport->Size);
        

            ImGuiWindowFlags flags =
                ImGuiWindowFlags_NoTitleBar
                | ImGuiWindowFlags_NoCollapse
                | ImGuiWindowFlags_NoResize
                | ImGuiWindowFlags_NoMove
                | ImGuiWindowFlags_NoSavedSettings;
            ImGui::Begin("Main", nullptr, flags);
            ImGui::BeginTabBar("MC-Simulation");
            if (ImGui::BeginTabItem("Control"))
            {
                bool correct = true;

                ImGui::InputInt("Photons", &photon_count);
                if (photon_count < 1) photon_count = 1;
                ImGui::InputInt("Size depth-distance map", &size_deptdist_hm);

                ImGui::SeparatorText("Tissue Layers");
                ImGui::Checkbox("Save denisty result", &saveDR);
                if (ImGui::Button("Clear saved denisty result")) {
                    savedDenisty.clear();
                }

                if (ImGui::Button("+ Add Layer")) {
                    userLayers.emplace_back(10.0, 0.1, 0.9, 1.4, 3.0);
                }
                ImGui::SameLine();
                if (ImGui::Button("Clear All Layers") && !simulation_running) {
                    userLayers.clear();
                }
                auto& n_ext = n_external;
                ImGui::InputDouble("n external enviroment", &n_ext, 0.5, 1.0);
                auto& n_dep = n_depth;
                ImGui::InputDouble("n depth", &n_dep, 0.5, 1.0);
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
                        Biotissue tissue = buildTissueFromUI(max_deep);
                        Photon init_photon(0, 0, 0, 0, 0, 1, 1);
                        ready = false;
                        simulation_running = true;
                        progress = 0.0f;
                        defer = 1;
                        std::thread sim_thread([tissue, init_photon]() {
                            runSimulationDetectedWithProgress(tissue, init_photon, photon_count, trajectories, max_deep_photons, detected, n_external, n_depth);
                            ready = true;
                            simulation_running = false;
                            table_need_recalculation = true;
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
                if (ready && table_need_recalculation) {
                    table.assign(size_deptdist_hm, std::vector<int>(size_deptdist_hm, 0));
                    for (auto& maxD : max_deep_photons)
                    {
                        double deep_didx = maxD.first / (max_deep / size_deptdist_hm);
                        double dist_didx = maxD.second / (max_distance / size_deptdist_hm);
                        int di = static_cast<int>(deep_didx);
                        int dj = static_cast<int>(dist_didx);
                        if (di >= 0 && di < size_deptdist_hm && dj >= 0 && dj < size_deptdist_hm)
                            table[di][dj] += 1;
                    }
                    table_need_recalculation = false;
                }
                ImGui::EndTabItem();
            }
            if (ImGui::BeginTabItem("Trajectories"))
            {
                enum class Surface { XY, XZ, YZ };
                static Surface sur = Surface::XZ;
                if (ImGui::Selectable("XY", sur == Surface::XY))  sur = Surface::XY;
                if (ImGui::Selectable("XZ", sur == Surface::XZ)) sur = Surface::XZ;
                if (ImGui::Selectable("YZ", sur == Surface::YZ))  sur = Surface::YZ;
                if (ready && (sur == Surface::XZ) &&ImPlot::BeginPlot("Monte Carlo paths (XZ)", ImVec2(-1, -1))) {
                    std::lock_guard<std::mutex> lock(traj_mutex);
                    ImPlot::SetupAxis(ImAxis_Y1, "axis Z (mm)");
                    ImPlot::SetupAxis(ImAxis_X1, "axis X (mm)");
                    std::vector<double> end_x, end_z;
                    for (int i = 0; i < std::min((int)trajectories.size(), 1000); i++) {
                        const auto& path = trajectories[i];
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
                        ImPlot::PlotScatter("end points", end_x.data(), end_z.data(), end_x.size());
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
                if (ready && sur == Surface::XY && ImPlot::BeginPlot("Monte Carlo paths (XY)", ImVec2(-1, -1))) {
                    std::lock_guard<std::mutex> lock(traj_mutex);
                    ImPlot::SetupAxis(ImAxis_Y1, "axis Y (mm)");
                    ImPlot::SetupAxis(ImAxis_X1, "axis X (mm)");
                    std::vector<double> end_x, end_y;
                    for (int i = 0; i < std::min((int)trajectories.size(), 1000); i++) {
                        const auto& path = trajectories[i];
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
                if (ready && (sur == Surface::YZ) && ImPlot::BeginPlot("Monte Carlo paths (YZ)", ImVec2(-1, -1))) {
                    std::lock_guard<std::mutex> lock(traj_mutex);
                    ImPlot::SetupAxis(ImAxis_Y1, "axis Z (mm)");
                    ImPlot::SetupAxis(ImAxis_X1, "axis Y (mm)");
                    std::vector<double> end_y, end_z;
                    for (int i = 0; i < std::min((int)trajectories.size(), 1000); i++) {
                        const auto& path = trajectories[i];
                        if (path.empty()) continue;
                        std::vector<double> y, z;
                        y.reserve(path.size());
                        z.reserve(path.size());
                        for (const auto& p : path) {
                            y.push_back(p.y);
                            z.push_back(p.z);
                        }
                        ImPlot::PlotLine("path", y.data(), z.data(), (int)y.size());
                        const auto& last = path.back();
                        end_y.push_back(last.y);
                        end_z.push_back(last.z);
                        i++;
                    }
                    if (!end_y.empty()) {
                        ImPlot::PlotScatter("end points", end_y.data(), end_z.data(), (int)end_y.size());
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
                ImGui::EndTabItem();
            }
            if (ImGui::BeginTabItem("Denisty Plot"))
            {
                enum class Scale { Linear, Logarithmical };
                static Scale sc = Scale::Linear;
                if (ImGui::Selectable("Linear", sc == Scale::Linear))  sc = Scale::Linear;
                if (ImGui::Selectable("Logarithmical", sc == Scale::Logarithmical)) sc = Scale::Logarithmical;
                ImGui::PushFont(pMyBigFont);
                if (ready && sc==Scale::Linear && ImPlot::BeginPlot("Monte Carlo Denisty plot", ImVec2(-1, -1))) {
                    std::lock_guard<std::mutex> lock(traj_mutex);
                    ImPlot::SetupAxis(ImAxis_Y1, "density");
                    ImPlot::SetupAxis(ImAxis_X1, "radius (mm)");
                    std::vector<double> x, y;
                    x.reserve(detected.size());
                    y.reserve(detected.size());
                    for (int j = 0; j < (int)(max_distance / step); j++) {
                        x.push_back(j * step);
                        y.push_back(detected[j]);

                    }
                    if (!savedDenisty.empty()) {
                        for (int k = 0; k < (int)savedDenisty.size(); k++) {
                            std::string namePath = std::to_string(k + 1) + "path";
                            ImPlot::PlotLine(namePath.c_str(), x.data(), savedDenisty[k].data(), (int)x.size());
                        }
                    }
                    if (saveDR && (defer > 0)) {
                        savedDenisty.push_back(y);
                        defer--;
                    }
                    ImPlot::PlotLine("newpath", x.data(), y.data(), (int)x.size());
                    ImPlot::EndPlot();
                }
                if (ready && ImPlot::BeginPlot("Monte Carlo Denisty plot (logariphm)", ImVec2(-1, -1))) {
                    std::lock_guard<std::mutex> lock(traj_mutex);

                    ImPlot::SetupAxis(ImAxis_Y1, "density");
                    ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);
                    ImPlot::SetupAxis(ImAxis_X1, "radius (mm)");

                    std::vector<double> x, y;
                    x.reserve(detected.size());
                    y.reserve(detected.size());
                    for (int j = 0; j < (int)(max_distance / step); j++) {
                        x.push_back(j * step);
                        y.push_back(detected[j]);

                    }
                    if (!savedDenisty.empty()) {
                        for (int k = 0; k < (int)savedDenisty.size(); k++) {
                            std::string namePath = std::to_string(k + 1) + "path";
                            ImPlot::PlotLine(namePath.c_str(), x.data(), savedDenisty[k].data(), (int)x.size());
                        }
                    }
                    if (saveDR && (defer > 0)) {
                        savedDenisty.push_back(y);
                        defer--;
                    }
                    ImPlot::PlotLine("newpath", x.data(), y.data(), (int)x.size());
                    ImPlot::EndPlot();
                }
                else if (!ready && !simulation_running && trajectories.empty()) {
                    ImGui::Text("Press 'Run simulation' to start.");
                }
                ImGui::PopFont();
                ImGui::EndTabItem();
            }
            if (ImGui::BeginTabItem("Dependence of the max depth")) {
                ImGuiTableFlags flags =
                    ImGuiTableFlags_Resizable
                    | ImGuiTableFlags_SizingFixedFit
                    | ImGuiTableFlags_ScrollX 
                    | ImGuiTableFlags_ScrollY;

                if (ImGui::BeginTable("table", size_deptdist_hm +1, flags)) {
                    ImGui::TableSetupScrollFreeze(1, 1);

                    ImGui::TableSetupColumn("Depth/Distance", 0, 128.0f);
                    for (int i = 0; i < size_deptdist_hm; i++) {
                        double dist = (static_cast<double>(i) + 0.5) * (max_distance / size_deptdist_hm);
                        ImGui::TableSetupColumn(std::to_string(dist).c_str(), 0, 40.0f);
                    }
                    ImGui::TableHeadersRow();
                    for (int i = 0; i < size_deptdist_hm; i++) {
                        ImGui::TableNextRow();
                        ImGui::TableSetColumnIndex(0);
                        double depth = (i + 0.5) * (max_deep / size_deptdist_hm);
                        ImGui::Text("%.2f", depth);
                        for (int j = 0; j < size_deptdist_hm; j++) {
                            ImGui::TableSetColumnIndex(j+1);
                            ImGui::Text("%d", table[i][j]);
                        }
                    }
                    ImGui::EndTable();
                }
                ImGui::EndTabItem();
            
            }
            if (ImGui::BeginTabItem("Dependence of the max depth(image)")) {
                ImPlot::PushColormap(ImPlotColormap_Jet);
                if (ImPlot::BeginPlot("Heatmap", ImVec2(-1,-1))) {
                    ImPlot::SetupAxes("X-Axis", "Y-Axis");
                    std::vector<double> flat;
                    flat.reserve(size_deptdist_hm* size_deptdist_hm);
                    for (int i = 0; i < size_deptdist_hm; ++i){
                        for (int j = 0; j < size_deptdist_hm; ++j){
                            flat.push_back(table[i][j]);
                        }
                    }
                    std::vector<double> sorted = flat;
                    std::sort(sorted.begin(), sorted.end());
                    int hi = sorted[(int)(sorted.size()*0.99)];
                    if (hi < 1) hi = 1;
                    ImPlot::PlotHeatmap("Matrix Data", flat.data(), size_deptdist_hm, size_deptdist_hm, 0.0f, hi, nullptr);
                    ImPlot::EndPlot();
                }
                ImGui::EndTabItem();
                ImPlot::PopColormap();
            }
            ImGui::EndTabBar();
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