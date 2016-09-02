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

#include "RandomNumbs.h"
#include "Tagging_system.h"
#include "Environment.h"
#include "DataHandler.h"

/**
 * @brief A handful of tips about the input parameters.
 */
void printTipsToRun(){
    std::cout << std::endl;
    std::cout << "Parameters should be:" << std::endl;
    std::cout << " 1. Seed for the RNG (when set to < 0 the program will " <<
            "seed the RNG engine itself with a truly random number)." << std::endl;
    std::cout << " 2. Number of bits in a MHC gene." << std::endl;
    std::cout << " 3. Host population size." << std::endl;
    std::cout << " 4. Number of genes in one host chromosome (they have " <<
            "two chromosomes)." << std::endl;
    std::cout << " 5. Number of host generations (effective length of model run)." <<
            std::endl;
    std::cout << " 6. Probability of mutation in hosts ([0,1] range)." << std::endl;
    std::cout << " 7. The heterozygote advantage / lack of advantage " <<
                "mode. It has to be 10 for heterozygote advantage or 11 for " <<
                "lack of thereof." << std::endl;
    std::cout << " 8. Probability of deleting a gene in the host ([0,1] range)." <<
            std::endl;
    std::cout << " 9. Probability of duplicating a gene in the host" <<
                " ([0,1] range)" << std::endl;
    std::cout << "10. Maximal number of genes permitted in one host chromosome." <<
            std::endl;
    std::cout << "11. Number of potential sexual partners an individual courts " <<
            "before picking one to mate with." <<
            std::endl;
    std::cout << std::endl;
}


/**
 * @brief The main function. Things are happening here.
 *
 * Compile this program with:
 * g++ -O3 -o MHC_model main.cpp Gene.cpp Antigen.cpp Host.cpp Pathogen.cpp H2Pinteraction.cpp RandomNumbs.cpp Tagging_system.cpp Environment.cpp DataHandler.cpp -std=c++1y
 *
 * @param argc - number of arguments
 * @param argv - list of arguments
 * @return 0
 */
int main(int argc, char** argv) {
// === Check if the entered parameters make sense ===
    int numbOfArgs = 12; // how many arguments we need to run this model
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
    int rndSeed, mhcGeneLength, hostPopSize,  hostGeneNumbb,  numOfHostGenerations,
        HeteroHomo, maxGene, numOfMates;
    double hostMutationProb, deletion, duplication;
    // Check if input params are numbers
    try {
        rndSeed = boost::lexical_cast<int>(argv[1]);
        mhcGeneLength = boost::lexical_cast<int>(argv[2]);
        hostPopSize = boost::lexical_cast<int>(argv[3]);
        hostGeneNumbb = boost::lexical_cast<int>(argv[4]);
        numOfHostGenerations = boost::lexical_cast<int>(argv[5]);
        hostMutationProb = boost::lexical_cast<double>(argv[6]);
        HeteroHomo = boost::lexical_cast<int>(argv[7]);
        deletion = boost::lexical_cast<double>(argv[8]);
        duplication = boost::lexical_cast<double>(argv[9]);
        maxGene = boost::lexical_cast<int>(argv[10]);
        numOfMates = boost::lexical_cast<int>(argv[11]);
    }
    catch(boost::bad_lexical_cast& e) {
        std::cout << std::endl;
        std::cout << "Arguments from 1 to " << numbOfArgs-1 << " should be " <<
            "numbers. Not all are numbers. Check the params list!" << std::endl;
        printTipsToRun();
        return 0;
    }
    // Load the input params
    rndSeed = atoi(argv[1]);
    mhcGeneLength = atoi(argv[2]);
    hostPopSize = atoi(argv[3]);
    hostGeneNumbb = atoi(argv[4]);
    numOfHostGenerations = atoi(argv[5]);
    hostMutationProb = atof(argv[6]);
    HeteroHomo = atoi(argv[7]);
    deletion = atof(argv[8]);
    duplication = atof(argv[9]);
    maxGene = atoi(argv[10]);
    numOfMates = atoi(argv[11]);

    // When told so, fetching a truly random number to seed the RNG engine
    if (rndSeed < 0){
        std::random_device rd;
        std::uniform_int_distribution<int> dist(0, 99999);
        rndSeed = dist(rd);
    }

    DataHandler Data2file;  // Initialize the data harvesting mechanism

    // Check if input params are of any sense
    if (Data2file.checkParamsIfWrong(rndSeed, mhcGeneLength, hostPopSize,
            hostGeneNumbb, numOfHostGenerations, hostMutationProb,
            HeteroHomo, deletion, duplication, maxGene, numOfMates)){
        std::cout << std::endl;
        std::cout << "Error in parameters on input. Check them." << std::endl;
        printTipsToRun();
        return 0;
    }
    std::cout << std::endl;
    std::cout << "Everything seems fine. Running the model." << std::endl;

    // Save input parameters to file
    Data2file.inputParamsToFile(rndSeed, mhcGeneLength, hostPopSize,
            hostGeneNumbb, numOfHostGenerations, hostMutationProb,
            HeteroHomo, deletion, duplication, maxGene, numOfMates);

// === And now doing the calculations! ===

    // Initializing the random number generator engine and tagging system
    RandomNumbs* p_RandomNumbs = RandomNumbs::getInstance();
    p_RandomNumbs->SetSeed(rndSeed);
    Tagging_system* pTagging_system = Tagging_system::getInstance();
    pTagging_system->setValue(0);

// ======================================

    Environment ENV; // Initialize the simulation environment
    Data2file.setAllFilesAsFirtsTimers();
    ENV.setHostRandomPopulation(hostPopSize, mhcGeneLength, hostGeneNumbb, 0);
//    ENV.setHostClonalPopulation(hostPopSize, mhcGeneLength, hostGeneNumbb, 0);
    std::cout << "Host population all set!" << std::endl;
    hostMutationProb = ENV.MMtoPMscaling(hostMutationProb, mhcGeneLength);
    std::ofstream InputParams;
    InputParams.open("InputParameters.csv", std::ios::out | std::ios::ate | std::ios::app);
    InputParams << "# Other_information:" << std::endl;
    InputParams << "\tseparated_species_genomes = YES" << std::endl;
    // set "NO" when using ENV.setPathoPopulatioUniformGenome()
//    InputParams << "\tseparated_species_genomes = NO" << std::endl;
    InputParams << "\tpoint_mutation_in_host_is_used = " << hostMutationProb << std::endl;
    InputParams << std::endl;
    InputParams.close();
    std::cout << "Calculating...." << std::endl;
    // Heterozygote advantage
    if(HeteroHomo == 10){
        Data2file.saveHostPopulToFile(ENV, 0);
        Data2file.saveHostGeneticDivers(ENV, 0);
        Data2file.saveHostGeneNumbers(ENV, 0);
        for(int i = 1; i <= numOfHostGenerations; ++i){
//            ENV.matingWithOneDifferentMHC();
            ENV.matingWithNoCommonMHCsmallSubset(numOfMates);
            ENV.mutateHostsWithDelDupl(hostMutationProb, deletion, duplication,
                                       maxGene, i);
            Data2file.saveHostGeneticDivers(ENV, i);
            Data2file.saveHostGeneNumbers(ENV, i);
            ENV.clearHostInfectionsData();
//           std::cout << "Host loop " << i << " finished" << std::endl;
        }
    } else {
       std::cout << "This instance of the model allows only heterozygote";
       std::cout << " advantage. Sorry :-(" << std::endl;
       return 0;
    }
    Data2file.saveHostPopulToFile(ENV, numOfHostGenerations);

    std::cout << "Run finished. Check the output files for results." << std::endl;

    std::cout << std::endl;

    return 0;
}
