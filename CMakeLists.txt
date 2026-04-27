cmake_minimum_required(VERSION 3.16)

project(mc_simulation LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ================= ENGINE =================
add_subdirectory(cpp/src)

# ================= UI =================
add_executable(mc_ui gui/main.cpp)

# ================= SYSTEM LIBS =================
find_package(OpenGL REQUIRED)
find_package(glfw3 CONFIG REQUIRED)

# ================= IMGUI (FULL MANUAL ONLY) =================
set(IMGUI_DIR ${CMAKE_SOURCE_DIR}/external/imgui)

set(IMGUI_SRC
    ${IMGUI_DIR}/imgui.cpp
    ${IMGUI_DIR}/imgui_draw.cpp
    ${IMGUI_DIR}/imgui_tables.cpp
    ${IMGUI_DIR}/imgui_widgets.cpp
    ${IMGUI_DIR}/imgui_demo.cpp

    ${IMGUI_DIR}/backends/imgui_impl_glfw.cpp
    ${IMGUI_DIR}/backends/imgui_impl_opengl3.cpp
)

# ================= IMPLOT =================
file(GLOB IMPLOT_SRC
    external/implot/*.cpp
)

target_sources(mc_ui PRIVATE
    ${IMGUI_SRC}
    ${IMPLOT_SRC}
)

target_include_directories(mc_ui PRIVATE
    cpp/include
    external/imgui
    external/imgui/backends
    external/implot
)

target_link_libraries(mc_ui PRIVATE
    mc_simulation
    glfw
    OpenGL::GL
)

set_target_properties(mc_ui PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)