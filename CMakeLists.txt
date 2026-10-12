cmake_minimum_required(VERSION 3.20)
project(ism_swirl_dynamics CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Option to enable/disable SFML visualization (useful for headless servers)
option(ENABLE_VISUALIZATION "Build the interactive SFML visualizer target" ON)

# Find GSL (Required for stochastic forces and random generators)
find_package(GSL REQUIRED)

# 1. Core Physics Simulation Library (Headless-safe)
add_library(ism_core
    src/pvector.cpp
    src/bird.cpp
    src/flock.cpp
    src/simulation.cpp
    src/forcefield.cpp
)
target_include_directories(ism_core PUBLIC include)
target_link_libraries(ism_core PUBLIC GSL::gsl GSL::gslcblas)

# 2. Headless Batch Simulation Executable (No SFML required)
add_executable(ism_batch src/main.cpp)
target_link_libraries(ism_batch PRIVATE ism_core)

# 3. Optional Graphical Visualizer Executable (Requires SFML)
if(ENABLE_VISUALIZATION)
    find_package(SFML COMPONENTS graphics window system REQUIRED)
    
    add_executable(ism_visual
        src/visual_main.cpp
        src/manager.cpp
    )
    target_link_libraries(ism_visual PRIVATE 
        ism_core 
        sfml-graphics 
        sfml-window 
        sfml-system
    )
endif()