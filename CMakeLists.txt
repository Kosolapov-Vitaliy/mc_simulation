cmake_minimum_required(VERSION 3.16)

project(mc_simulation LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ----- Автоматическая загрузка GLFW -----
include(FetchContent)
FetchContent_Declare(
  glfw
  GIT_REPOSITORY https://github.com/glfw/glfw.git
  GIT_TAG        3.3.8   # или использовать ветку/тег, например 3.3.8
)
# Отключаем ненужное, чтобы не тратить время
set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(glfw)
# -----------------------------------------

# ================= ENGINE =================
add_subdirectory(cpp/src)

# ================= UI =================
add_executable(mc_ui gui/main.cpp)

# ================= SYSTEM LIBS =================
find_package(OpenGL REQUIRED)  # Только OpenGL остаётся системной зависимостью

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
    glfw              # цель, созданная FetchContent
    OpenGL::GL
)

set_target_properties(mc_ui PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)