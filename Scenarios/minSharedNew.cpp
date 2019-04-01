/*
 * File:   main.cpp
 * Author: Piotr Bentkowski : bentkowski.piotr@gmail.com
 *
 * Created on 15 March 2016, 19:50
 *
 *    This program is free software; you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation; either version 2 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *    MA 02110-1301, USA.
 */

#include <cstdlib>
#include <iostream>
#include <random>
#include <boost/lexical_cast.hpp>
#include <fstream>
#include <thread>     // for reading the number of concurrent threads supported

#include "omp.h"

#include "src/Tagging_system.h"
#include "src/Random.h"
#include "src/Environment.h"
#include "src/DataHandler.h"
#include "src/nlohmann/json.hpp"

using jsonf = nlohmann::json;

/**
 * @brief A handful of tips about the input parameters.
 */
void printTipsToRun(){
    std::cout << std::endl;
    std::cout << "This is the first sex scenario where most different MHC composition."
            " Parameters should be:" << std::endl;
    std::cout << " 1. The number of threads the program will use. Give 0 to use all the available CPU cores." << std::endl;
    std::cout << " 2. Number of bits in a MHC gene." << std::endl;
    std::cout << " 3. Number of bits in an antigen." << std::endl;
    std::cout << " 4. Host population size." << std::endl;
    std::cout << " 5. Pathogen population size." << std::endl;
    std::cout << " 6. Number of pathogen species." << std::endl;
    std::cout << " 7. Number of genes in one host chromosome (they have " <<
            "two chromosomes)." << std::endl;
    std::cout << " 8. Number of pathogen generations per one host generation. " <<
            std::endl;
    std::cout << " 9. Number of host generations (effective length of model run)." <<
            std::endl;
    std::cout << "10. Probability of mutation in hosts ([0,1] range)." << std::endl;
    std::cout << "11. Probability of mutation in pathogens ([0,1] range)." << std::endl;
    std::cout << "12. The heterozygote advantage / lack of advantage " <<
                "mode. It has to be 10 for heterozygote advantage or 11 for " <<
                "lack of thereof." << std::endl;
    std::cout << "13. Probability of deleting a gene in the host ([0,1] range)." <<
            std::endl;
    std::cout << "14. Probability of duplicating a gene in the host" <<
                " ([0,1] range)" << std::endl;
    std::cout << "15. Maximal number of genes permitted in one host chromosome." <<
            std::endl;
    std::cout << "16. Number of sexual partners an individual checks out before selecting one for mating." <<
            std::endl;
    std::cout << "17. Alpha factor for the host fitness function ([0,1] range)." << std::endl;
    std::cout << std::endl;

}


/**
 * @brief The main function. Things are happening here.
 *
 * Compile this program with:
 * g++ -O3 -o MHC_model main.cpp nlohmann/json.hpp Gene.cpp Antigen.cpp Host.cpp Pathogen.cpp H2Pinteraction.cpp \
  Random.cpp Tagging_system.cpp Environment.cpp DataHandler.cpp -fopenmp -std=c++14
 *
 * @param argc - number of arguments
 * @param argv - list of arguments
 * @return 0
 */
int main(int argc, char** argv) {
// === Check if the entered parameters make sense ===
    int numbOfArgs = 18; // how many arguments we need to run this model
    if (argc < numbOfArgs) {
        std::cout << std::endl;
        std::cout << "Not enough arguments. It has to be " <<
            "precisely " << numbOfArgs -1 << " of them and only " << argc - 1 <<
            " are provided." << std::endl;
        printTipsToRun();
        return 0;
    }
    if (argc > numbOfArgs) {
        std::cout << std::endl;
        std::cout << "Too many arguments. It has to be " <<
            "precisely " << numbOfArgs -1 << " of them but " << argc - 1 <<
            " are provided." << std::endl;
        printTipsToRun();
        return 0;
    }
    unsigned long maxGene, hostGeneNumbb, mhcGeneLength, antigenLength;
    unsigned int numberOfThreads;
    int hostPopSize, pathoPopSize, patho_sp, NumbPartners,
        patoPerHostGeneration, numOfHostGenerations, HeteroHomo;
    double hostMutationProb, pathoMutationProb, deletion, duplication, alpha;
    // Check if input params are numbers
    try {
        numberOfThreads = boost::lexical_cast<unsigned int>(argv[1]);
        mhcGeneLength = boost::lexical_cast<unsigned long>(argv[2]);
        antigenLength = boost::lexical_cast<unsigned long>(argv[3]);
        hostPopSize = boost::lexical_cast<int>(argv[4]);
        pathoPopSize = boost::lexical_cast<int>(argv[5]);
        patho_sp = boost::lexical_cast<int>(argv[6]);
        hostGeneNumbb = boost::lexical_cast<unsigned long>(argv[7]);
        patoPerHostGeneration = boost::lexical_cast<int>(argv[8]);
        numOfHostGenerations = boost::lexical_cast<int>(argv[9]);
        hostMutationProb = boost::lexical_cast<double>(argv[10]);
        pathoMutationProb = boost::lexical_cast<double>(argv[11]);
        HeteroHomo = boost::lexical_cast<int>(argv[12]);
        deletion = boost::lexical_cast<double>(argv[13]);
        duplication = boost::lexical_cast<double>(argv[14]);
        maxGene = boost::lexical_cast<unsigned long>(argv[15]);
        NumbPartners = boost::lexical_cast<int>(argv[16]);
        alpha = boost::lexical_cast<double>(argv[17]);
    }
    catch(boost::bad_lexical_cast& e) {
        std::cout << std::endl;
        std::cout << "Arguments from 1 to " << numbOfArgs-1 << " should be " <<
            "numbers. Not all are numbers. Check the params list!" << std::endl;
        printTipsToRun();
        return 0;
    }
    // Load the input params
    numberOfThreads = (unsigned int) strtol(argv[1], nullptr, 10);
    mhcGeneLength = (unsigned long) strtol(argv[2], nullptr, 10);
    antigenLength = (unsigned long) strtol(argv[3], nullptr, 10);
    hostPopSize = (int) strtol(argv[4], nullptr, 10);
    pathoPopSize = (int) strtol(argv[5], nullptr, 10);
    patho_sp = (int) strtol(argv[6], nullptr, 10);
    hostGeneNumbb = (unsigned long) strtol(argv[7], nullptr, 10);
    patoPerHostGeneration = (int) strtol(argv[8], nullptr, 10);
    numOfHostGenerations = (int) strtol(argv[9], nullptr, 10);
    hostMutationProb = strtod(argv[10], nullptr);
    pathoMutationProb = strtod(argv[11], nullptr);
    HeteroHomo = (int) strtol(argv[12], nullptr, 10);
    deletion = strtod(argv[13], nullptr);
    duplication = strtod(argv[14], nullptr);
    maxGene = (unsigned long) strtol(argv[15], nullptr, 10);
    NumbPartners = (int) strtol(argv[16], nullptr, 10);
    alpha = strtod(argv[17], nullptr);

    unsigned int threadsAvailable = std::thread::hardware_concurrency();
    // Initializing the multi-threaded environment
    if (numberOfThreads == 0)
    {
        numberOfThreads = std::thread::hardware_concurrency();
        if(numberOfThreads == 0) // if the value is not well defined or not computable, set at least 1 thread
            numberOfThreads = 1;
    } else if (numberOfThreads > threadsAvailable) {
        numberOfThreads = threadsAvailable;
    }
    std::cout << "We have " << numberOfThreads << " threads" << std::endl;
    omp_set_num_threads(numberOfThreads);


    DataHandler Data2file;  // Initialize the data harvesting mechanism

    // Check if input params are of any sense
    if (Data2file.checkParamsIfWrong(numberOfThreads, mhcGeneLength, antigenLength, hostPopSize,
            pathoPopSize, patho_sp, hostGeneNumbb, patoPerHostGeneration, numOfHostGenerations,
            hostMutationProb, pathoMutationProb, HeteroHomo, deletion, duplication,
            maxGene, alpha, NumbPartners)){
        std::cout << std::endl;
        std::cout << "Error in parameters on input. Check them." << std::endl;
        printTipsToRun();
        return 0;
    }
    std::cout << std::endl;
    std::cout << "Everything seems fine. Running the model." << std::endl;

    // Save input parameters to file
    Data2file.inputParamsToFile(numberOfThreads, mhcGeneLength, antigenLength, hostPopSize,
            pathoPopSize, patho_sp, hostGeneNumbb, patoPerHostGeneration, numOfHostGenerations, hostMutationProb,
            pathoMutationProb, HeteroHomo, deletion, duplication, maxGene,
            alpha, NumbPartners);

// === And now doing the calculations! ===

    // Initializing tagging system
    Tagging_system tag;
// ======================================

    Environment ENV(numberOfThreads); // Initialize the simulation environment
    Data2file.setAllFilesAsFirtsTimers();
    ENV.setHostRandomPopulation(hostPopSize, mhcGeneLength, hostGeneNumbb, 0, tag);
//    ENV.setHostClonalPopulation(hostPopSize, mhcGeneLength, hostGeneNumbb, 0);
    std::cout << "Host population all set!" << std::endl;
    ENV.setPathoPopulatioDivSpecies(pathoPopSize, antigenLength, patho_sp, mhcGeneLength, 0, 0, tag);
    std::cout << "Pathogen population all set!" << std::endl;
    hostMutationProb = ENV.MMtoPMscaling(hostMutationProb, mhcGeneLength);
    Data2file.savePathoNoMuttList(ENV);

    // Adding extra info about the parameters of this simulation
    std::ifstream inJson("InputParameters.json");
    jsonf jsonfile;
    inJson >> jsonfile;
    jsonfile["separated_species_genomes"] = "YES";
    // set "NO" when using ENV.setPathoPopulatioUniformGenome()
//    jsonfile["separated_species_genomes"] = "NO";
    jsonfile["point_mutation_in_host_is_used"] = hostMutationProb;
    std::string s = jsonfile.dump(4);
    std::ofstream InputParams;
    InputParams.open("InputParameters.json");
    InputParams << s;
    InputParams.close();
    
    std::cout << "Calculating...." << std::endl;
    // Heterozygote advantage
    if(HeteroHomo == 10){
//        Data2file.savePathoPopulToFile(ENV, 0);
        Data2file.saveHostPopulToFile(ENV, 0);
        Data2file.saveHostGeneticDivers(ENV, 0);
        Data2file.saveMhcNumbersBeforeMating(ENV, 0);
        Data2file.saveMhcNumbersWhenMating(ENV, 0);
        Data2file.saveMhcNumbersAfterMating(ENV, 0);
	Data2file.savePresentedPathos(ENV, 0);
        for(int i = 1; i <= numOfHostGenerations; ++i){
            for(int j = 0; j < patoPerHostGeneration; ++j){
                ENV.infectOneFromOneSpecHetero();
                ENV.selectAndReproducePathoFixedPopSizes();
                ENV.mutatePathogensWithRestric(pathoMutationProb, mhcGeneLength, i, tag);
                ENV.clearPathoInfectionData();
            }
            Data2file.savePresentedPathos(ENV, i);
            ENV.calculateHostsFitnessPlainPresent();
//            ENV.selectAndReprodHostsReplace();
            ENV.selectAndReprodHostsNoMating();  // changed for sexual reproduction
            Data2file.saveMhcNumbersBeforeMating(ENV, i);
            ENV.matingWithNoCommonMHCsmallSubset(NumbPartners); // changed for sexual reproduction
            Data2file.saveMhcNumbersWhenMating(ENV, i);
            Data2file.saveMhcNumbersAfterMating(ENV, i);
            ENV.mutateHostsWithDelDuplPointMuts(hostMutationProb, deletion, duplication, maxGene, i, tag);
            Data2file.saveHostGeneticDivers(ENV, i);
//            Data2file.saveHostGeneNumbers(ENV, i);
            ENV.clearHostInfectionsData();
//           std::cout << "Host loop " << i << " finished" << std::endl;
        }
        ENV.infectOneFromOneSpecHetero();
    } else {
       std::cout << "This instance of the model allows only heterozygote";
       std::cout << " advantage. Sorry :-(" << std::endl;
       return 0;
    }
//    Data2file.savePathoPopulToFile(ENV, numOfHostGenerations);
    Data2file.saveHostPopulToFile(ENV, numOfHostGenerations);

    std::cout << "Run finished. Check the output files for results." << std::endl;

    std::cout << std::endl;

    return 0;
}
