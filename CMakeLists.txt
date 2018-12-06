cmake_minimum_required(VERSION 3.6)
project(MHC_code_OBA)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")

set(SOURCE_FILES
    Antigen.cpp
    Antigen.h
    DataHandler.cpp
    DataHandler.h
    Environment.cpp
    Environment.h
    Gene.cpp
    Gene.h
    H2Pinteraction.cpp
    H2Pinteraction.h
    Host.cpp
    Host.h
    main.cpp
    mainpage.h
    Pathogen.cpp
    Pathogen.h
    Random.cpp
    Random.h
    Tagging_system.cpp
    Tagging_system.h mainTest.cpp)

add_executable(MHC_code_OBA ${SOURCE_FILES} mainTest.cpp)