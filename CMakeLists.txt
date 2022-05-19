cmake_minimum_required(VERSION 3.21)
project(HMM4GA_clion)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(HMM4GA_clion
        main.cpp
        Sequence.cpp Sequence.h SequenceGenerator.h SequenceGenerator.cpp StateTransitionMatrix.h StateTransitionMatrix.cpp
        PairHMM.cpp PairHMM.h ProbabilityMatrix.cpp ProbabilityMatrix.h DynamicMatrix.cpp DynamicMatrix.h)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(HMM4GA_clion PUBLIC OpenMP::OpenMP_CXX)
endif()


