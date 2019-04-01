cmake_minimum_required(VERSION 3.6)
project(MHC_code_OBA)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -std=c++14")

set(SOURCE_FILES
    src/Antigen.cpp
    src/Antigen.h
    src/DataHandler.cpp
    src/DataHandler.h
    src/Environment.cpp
    src/Environment.h
    src/Gene.cpp
    src/Gene.h
    src/H2Pinteraction.cpp
    src/H2Pinteraction.h
    src/Host.cpp
    src/Host.h
    main.cpp
    src/mainpage.h
    src/Pathogen.cpp
    src/Pathogen.h
    src/Random.cpp
    src/Random.h
    src/Tagging_system.cpp
    src/Tagging_system.h
    src/nlohmann/json.hpp)

add_executable(MHC_code_OBA ${SOURCE_FILES} main.cpp)