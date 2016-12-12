cmake_minimum_required(VERSION 2.8)
set(CMAKE_CXX_STANDARD 11) # C++11...
set(CMAKE_CXX_STANDARD_REQUIRED ON) #...is required...
set(CMAKE_CXX_EXTENSIONS ON) #...without compiler extensions like gnu++11
project(cuda_fphi CXX)

include_directories(${CMAKE_HOME_DIRECTORY}/src)

include_directories(${CMAKE_HOME_DIRECTORY}/include)
link_directories(${CMAKE_HOME_DIRECTORY}/lib)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
set(
    CMAKE_RUNTIME_OUTPUT_DIRECTORY
    ${CMAKE_HOME_DIRECTORY}/bin
    )
add_executable(cuda_fphi src/cuda_fphi_main.cc src/write_to_hdf5.cc)
target_link_libraries(cuda_fphi -lcuda_fphi -lcudart_static -lcublas_static -lculibos -lhdf5_cpp -lhdf5 -lz -lncurses -ltinfo -lgpm -pthread -ldl -lrt )



