
//
// Created by piotr on 04/12/18.
// compile: g++ Tagging_system.cpp Gene.cpp Antigen.cpp Random.cpp Host.cpp Pathogen.cpp H2Pinteraction.cpp Environment.cpp DataHandler.cpp mainTestEnv.cpp -fopenmp -std=c++14
//     run: ./a.out
//

#include <vector>
#include <iostream>
#include <thread>     // for reading the number of concurrent threads supported
#include "omp.h"
#include <time.h>

//#include "Random.h"
#include "Tagging_system.h"
#include "DataHandler.h"
//#include "Gene.h"
//#include "Host.h"
#include "Environment.h"



int main(int argc, char** argv)
{
    time_t begin_t, end_t;
    begin_t = time(nullptr);
    unsigned int numberOfThreads = 0;

    Tagging_system tag;
    DataHandler Data2file;  // Initialize the data harvesting mechanism
    Environment ENV(numberOfThreads);

    int hostPopSize = 1000;
    unsigned long mhcGeneLength = 16;
    unsigned long hostGeneNumbb = 3;

    ENV.setHostRandomPopulation(hostPopSize, mhcGeneLength, hostGeneNumbb, 0, tag);
    ENV.setPathoPopulatioDivSpecies(16000, 6000, 16, 16, 0, 0.33, tag);
    ENV.mutatePathogensWithRestric(0.05, 16, 0, tag);
    ENV.infectOneFromOneSpecHetero();
    ENV.selectAndReproducePathoFixedPopSizes();

    Data2file.setAllFilesAsFirtsTimers();
    Data2file.saveHostPopulToFile(ENV, 0);
    Data2file.savePathoPopulToFile(ENV, 0);
    Data2file.savePathoNoMuttList(ENV);

    ENV.clearPathoInfectionData();
    ENV.clearHostInfectionsData();

    end_t = time(nullptr);
    std::cout << "Czas obliczen: " << difftime(end_t, begin_t) << std::endl;

    return 0;
}
