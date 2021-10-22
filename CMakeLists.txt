cmake_minimum_required(VERSION 3.0)
project(graphmatching)

# compiler and flags
#add_definitions("-Wall -O3")
add_definitions("-Wall -O0 -g")
# -fno-limit-debug-info")

set(CMAKE_CXX_STANDARD 17)

# include directory and the libraries we will build
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
add_subdirectory(lib)
add_subdirectory(test)

# build
add_executable(graphmatching main.cc "include/LpApproxSMFQ.h" "include/Exact_Exponential_SMFQ.h" 
"include/DirectApproachHR2LQ.h" "include/Node.h" "include/Edge.h" "include/Graph.h" 
"include/HApproxSMFQ.h" "include/Convert_HR_to_HR2LQ.h")
target_link_libraries(graphmatching graph)

install(TARGETS graphmatching RUNTIME DESTINATION bin)
