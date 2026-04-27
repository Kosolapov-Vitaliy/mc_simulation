set(MCS_LIBRARY mc_simulation)

file(GLOB SRC
    "*.cpp"
)

add_library(${MCS_LIBRARY} STATIC ${SRC})

target_include_directories(${MCS_LIBRARY} PUBLIC
    ${CMAKE_SOURCE_DIR}/cpp/include
)