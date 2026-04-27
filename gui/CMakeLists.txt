cmake_minimum_required(VERSION 3.12)

project(mc_gui)

set(CMAKE_CXX_STANDARD 20)

find_package(OpenGL REQUIRED)

file(GLOB IMGUI_SRC
    imgui/*.cpp
    imgui/backends/imgui_impl_glfw.cpp
    imgui/backends/imgui_impl_opengl3.cpp
)

file(GLOB IMPLOT_SRC
    implot/*.cpp
)

add_executable(mc_ui
    main.cpp
    ${IMGUI_SRC}
    ${IMPLOT_SRC}
)

target_include_directories(mc_ui PRIVATE
    ../cpp/include
    imgui
    imgui/backends
    implot
)

target_link_libraries(mc_ui
    ${MCS_LIBRARY}   
    glfw
    OpenGL::GL
)