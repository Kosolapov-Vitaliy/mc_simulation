cmake_minimum_required(VERSION 3.12)
project(mc_simulation)

set(MCS_LIBRARY "mc_simulation")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "Configs" FORCE)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Выходные директории
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE})

# Путь к заголовкам C++
set(MCS_INCLUDE "${CMAKE_CURRENT_SOURCE_DIR}/cpp/include")
include_directories("${MCS_INCLUDE}")

# ---- Поиск Java и JNI ----
find_package(Java REQUIRED)
find_package(JNI REQUIRED)
if(NOT JNI_FOUND)
    message(FATAL_ERROR "JNI не найден. Проверьте установку JDK.")
endif()
include_directories(${JNI_INCLUDE_DIRS})

# ---- Поиск JavaFX (через переменную JAVAFX_HOME) ----
if(DEFINED ENV{JAVAFX_HOME})
    set(JAVAFX_HOME $ENV{JAVAFX_HOME})
else()
    set(JAVAFX_HOME "" CACHE PATH "Path to JavaFX SDK (e.g., C:/javafx-sdk-21/lib)")
endif()

# Находим все нужные JAR-файлы JavaFX
set(JAVAFX_JARS "")
foreach(jar_name javafx.base javafx.controls javafx.fxml javafx.graphics)
    find_file(jar_file_${jar_name}
        NAMES ${jar_name}.jar ${jar_name}-controls.jar ${jar_name}.jar
        PATHS "${JAVAFX_HOME}"
        NO_DEFAULT_PATH
    )
    if(jar_file_${jar_name})
        list(APPEND JAVAFX_JARS ${jar_file_${jar_name}})
    else()
        message(WARNING "Missing JavaFX JAR: ${jar_name}.jar")
    endif()
endforeach()

if(JAVAFX_JARS)
    # Преобразуем список в строку с разделителем ; (для Windows)
    string(REPLACE ";" "\\;" JAVAFX_CLASSPATH "${JAVAFX_JARS}")
    message(STATUS "JavaFX classpath: ${JAVAFX_CLASSPATH}")
else()
    message(WARNING "No JavaFX JARs found. Set JAVAFX_HOME correctly.")
endif()

# ---- Java-исходники и ресурсы ----
set(JAVA_SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/java/src/main/java")
set(JAVA_RESOURCES_DIR "${CMAKE_CURRENT_SOURCE_DIR}/java/src/main/resources")
set(JAVA_CLASSES_DIR "${CMAKE_CURRENT_BINARY_DIR}/classes")
set(RESOURCES_DIR "${CMAKE_CURRENT_BINARY_DIR}/resources")
set(JNI_HEADER_DIR "${CMAKE_CURRENT_BINARY_DIR}/generated_headers")

# Копируем FXML-ресурс
file(COPY ${JAVA_RESOURCES_DIR}/com/mcgui/javagui/hello-view.fxml
    DESTINATION ${RESOURCES_DIR}/com/mcgui/javagui)

# ---- Генерация JNI-заголовка ----
add_custom_command(
    OUTPUT ${JNI_HEADER_DIR}/com_mcgui_javagui_MCSimulationJNI.h
    COMMAND ${CMAKE_COMMAND} -E make_directory ${JNI_HEADER_DIR}
    COMMAND ${Java_JAVAC_EXECUTABLE} -h ${JNI_HEADER_DIR}
        ${JAVA_SRC_DIR}/com/mcgui/javagui/MCSimulationJNI.java
    DEPENDS ${JAVA_SRC_DIR}/com/mcgui/javagui/MCSimulationJNI.java
    COMMENT "Generating JNI header"
)
add_custom_target(generate_jni_headers DEPENDS ${JNI_HEADER_DIR}/com_mcgui_javagui_MCSimulationJNI.h)

# ---- Компиляция Java с classpath JavaFX ----
set(JAVA_SOURCES
    ${JAVA_SRC_DIR}/com/mcgui/javagui/GUIController.java
    ${JAVA_SRC_DIR}/com/mcgui/javagui/HelloApplication.java
    ${JAVA_SRC_DIR}/com/mcgui/javagui/Launcher.java
    ${JAVA_SRC_DIR}/com/mcgui/javagui/LayerParams.java
    ${JAVA_SRC_DIR}/com/mcgui/javagui/MCSimulationJNI.java
)

add_custom_command(
    OUTPUT ${JAVA_CLASSES_DIR}/.timestamp
    COMMAND ${CMAKE_COMMAND} -E make_directory ${JAVA_CLASSES_DIR}
    COMMAND ${Java_JAVAC_EXECUTABLE}
        -cp ${JAVAFX_CLASSPATH}
        -d ${JAVA_CLASSES_DIR}
        ${JAVA_SOURCES}
    COMMAND ${CMAKE_COMMAND} -E touch ${JAVA_CLASSES_DIR}/.timestamp
    DEPENDS ${JAVA_SOURCES} ${JAVAFX_JARS} generate_jni_headers
    COMMENT "Compiling Java sources with JavaFX classpath"
)
add_custom_target(compile_java DEPENDS ${JAVA_CLASSES_DIR}/.timestamp)

# ---- Создание JAR с манифестом ----
set(MANIFEST_FILE "${CMAKE_CURRENT_BINARY_DIR}/MANIFEST.MF")
file(WRITE ${MANIFEST_FILE} "Main-Class: com.mcgui.javagui.Launcher\n")

add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/mc_gui-1.0.0.jar
    COMMAND ${Java_JAR_EXECUTABLE} cfm ${CMAKE_CURRENT_BINARY_DIR}/mc_gui-1.0.0.jar
        ${MANIFEST_FILE}
        -C ${JAVA_CLASSES_DIR} .
        -C ${RESOURCES_DIR} .
    DEPENDS compile_java
    COMMENT "Creating JAR file with manifest"
)
add_custom_target(mc_gui_jar DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/mc_gui-1.0.0.jar)

# ---- C++ статическая библиотека ----
add_subdirectory(cpp/src)

# ---- JNI-динамическая библиотека ----
add_library(mc_simulation_jni SHARED cpp/jni/jni_impl.cpp)

target_include_directories(mc_simulation_jni PRIVATE
    ${MCS_INCLUDE}
    ${JNI_HEADER_DIR}
    ${JNI_INCLUDE_DIRS}
)

target_link_libraries(mc_simulation_jni
    mc_simulation
    ${JNI_LIBRARIES}
)

if(WIN32)
    set_target_properties(mc_simulation_jni PROPERTIES PREFIX "")
endif()

add_dependencies(mc_simulation_jni generate_jni_headers)

# ---- Запуск приложения ----
# Формируем classpath для выполнения (наш JAR + все JavaFX JAR'ы)
set(RUNTIME_CLASSPATH "${JAVAFX_CLASSPATH};${CMAKE_CURRENT_BINARY_DIR}/mc_gui-1.0.0.jar")

add_custom_target(run
    COMMAND ${Java_JAVA_EXECUTABLE}
        -Djava.library.path=$<TARGET_FILE_DIR:mc_simulation_jni>
        -cp "${RUNTIME_CLASSPATH}"
        com.mcgui.javagui.Launcher
    DEPENDS mc_gui_jar mc_simulation_jni
    COMMENT "Running Monte Carlo GUI"
)

# ---- Отчёт ----
message(STATUS "")
message(STATUS "Configuration for ${PROJECT_NAME}")
message(STATUS "  Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "  JNI includes: ${JNI_INCLUDE_DIRS}")
message(STATUS "  JavaFX classpath: ${JAVAFX_CLASSPATH}")
message(STATUS "")