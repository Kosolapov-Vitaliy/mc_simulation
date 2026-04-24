set(target ${MCS_LIBRARY})

file(GLOB hdrs "*.h*" "${MCS_INCLUDE}/*.h*")
file(GLOB srcs "*.cpp")
add_library(${target} STATIC ${srcs} ${hdrs})
target_include_directories(${target} PUBLIC ${MCS_INCLUDE})