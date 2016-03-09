/* 
 * File:   main.cpp
 * Author: Piotr Bentkowski : bentkowski.piotr@gmail.com
 *
 * Created on 12 February 2015, 17:33
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
#include "DataHarvester.h"

/**
 * @brief A handful of tips about the input parameters.
 */
void printTipsToRun(){
    std::cout << std::endl;
    std::cout << "Parameters should be:" << std::endl;
    std::cout << " 1. Seed for the RNG (when set to < 0 the program will " <<
            "seed the RNG engine itself with a truly random number)." << std::endl;
    std::cout << " 2. Number of bits in a MHC gene." << std::endl;
    std::cout << " 3. Number of bits in an antigen." << std::endl;
    std::cout << " 4. Host population size." << std::endl;
    std::cout << " 5. Pathogen population size." << std::endl;
    std::cout << " 6. Number of pathogen species." << std::endl;
    std::cout << " 7. Number of genes in one host chromosome (they have " << 
            "two chromosomes)." << std::endl;
    std::cout << " 8. Number of antigens in a pathogen." << std::endl;
    std::cout << " 9. Number of pathogen generations per one host generation. " <<
            std::endl;
    std::cout << "10. Number of host generations (effective length of model run)." <<
            std::endl;
    std::cout << "11. Probability of mutation in hosts ([0,1] range)." << std::endl;
    std::cout << "12. Probability of mutation in pathogens ([0,1] range)." << std::endl;
    std::cout << "13. The heterozygote advantage / lack of advantage " <<
                "mode. It has to be 10 for heterozygote advantage or 11 for " <<
                "lack of thereof." << std::endl;
    std::cout << "14. Probability of deleting a gene in the host ([0,1] range)." << 
            std::endl;
    std::cout << "15. Probability of duplicating a gene in the host" <<
                " ([0,1] range)" << std::endl;
    std::cout << "16. Maximal number of genes permitted in one host chromosome." <<
            std::endl;
    std::cout << "17. Alpha factor for the host fitness function ([0,1] range)." <<
            std::endl;
    std::cout << "18. Fraction of antigen bits which get fixed and cannot" <<
            " mutate ([0,1] range)." << std::endl;
    std::cout << std::endl;
    
}


/**
 * @brief The main function. Things are happening here. 
 * 
 * Compile this program with:
 * g++ -O3 -o MHC_model main.cpp Gene.cpp Antigen.cpp Host.cpp Pathogen.cpp H2Pinteraction.cpp RandomNumbs.cpp Tagging_system.cpp Environment.cpp DataHarvester.cpp -std=c++1y
 * 
 * @param argc - number of arguments
 * @param argv - list of arguments
 * @return 0
 */
int main(int argc, char** argv) {
// === Check if the entered parameters make sense ===
    int numbOfArgs = 19; // how many arguments we need to run this model
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
    int rndSeed, mhcGeneLength, antigenLength, hostPopSize, pathoPopSize, patho_sp,
        hostGeneNumbb, pathoGeneNumb, patoPerHostGeneration, numOfHostGenerations,
        HeteroHomo, maxGene;
    double hostMutationProb, pathoMutationProb, deletion, duplication, alpha,
            fixedAntigPosit;
    // Check if input params are numbers
    try {
        rndSeed = boost::lexical_cast<int>(argv[1]);
        mhcGeneLength = boost::lexical_cast<int>(argv[2]);
        antigenLength = boost::lexical_cast<int>(argv[3]);
        hostPopSize = boost::lexical_cast<int>(argv[4]);
        pathoPopSize = boost::lexical_cast<int>(argv[5]);
        patho_sp = boost::lexical_cast<int>(argv[6]);
        hostGeneNumbb = boost::lexical_cast<int>(argv[7]);
        pathoGeneNumb = boost::lexical_cast<int>(argv[8]);
        patoPerHostGeneration = boost::lexical_cast<int>(argv[9]);
        numOfHostGenerations = boost::lexical_cast<int>(argv[10]);
        hostMutationProb = boost::lexical_cast<double>(argv[11]);
        pathoMutationProb = boost::lexical_cast<double>(argv[12]);
        HeteroHomo = boost::lexical_cast<int>(argv[13]);
        deletion = boost::lexical_cast<double>(argv[14]);
        duplication = boost::lexical_cast<double>(argv[15]);
        maxGene = boost::lexical_cast<int>(argv[16]);
        alpha = boost::lexical_cast<double>(argv[17]);
        fixedAntigPosit = boost::lexical_cast<double>(argv[18]);
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
    antigenLength = atoi(argv[3]);
    hostPopSize = atoi(argv[4]);
    pathoPopSize = atoi(argv[5]);
    patho_sp = atoi(argv[6]);
    hostGeneNumbb = atoi(argv[7]);
    pathoGeneNumb = atoi(argv[8]);
    patoPerHostGeneration = atoi(argv[9]);
    numOfHostGenerations = atoi(argv[10]);
    hostMutationProb = atof(argv[11]);
    pathoMutationProb = atof(argv[12]);
    HeteroHomo = atoi(argv[13]);
    deletion = atof(argv[14]);
    duplication = atof(argv[15]);
    maxGene = atoi(argv[16]);
    alpha = atof(argv[17]);
    fixedAntigPosit = atof(argv[18]);
    
    // When told so, fetching a truly random number to seed the RNG engine
    if (rndSeed < 0){
        std::random_device rd;
        std::uniform_int_distribution<int> dist(0, 99999);
        rndSeed = dist(rd);
    }
            
    DataHarvester Data2file;  // Initialize the data harvesting mechanism
    
    // Check if input params are of any sense
    if (Data2file.checkParamsIfWrong(rndSeed, mhcGeneLength, antigenLength, hostPopSize, 
            pathoPopSize, patho_sp, hostGeneNumbb, pathoGeneNumb,
            patoPerHostGeneration, numOfHostGenerations,
            hostMutationProb, pathoMutationProb, HeteroHomo, deletion, duplication,
            maxGene, alpha, fixedAntigPosit)){
        std::cout << std::endl;
        std::cout << "Error in parameters on input. Check them." << std::endl;
        printTipsToRun();
        return 0;
    }
    if (patho_sp % 4){
        std::cout << std::endl;
        std::cout << "This scenario needs the number of pathogen species " <<
                "to be multiplication of 4. Sorry." << std::endl;
        std::cout << std::endl;
        return 0;
    }
    std::cout << std::endl;
    std::cout << "Everything seems fine. Running the model." << std::endl;
    
    // Save input parameters to file
    Data2file.inputParamsToFile(rndSeed, mhcGeneLength, antigenLength, hostPopSize, 
            pathoPopSize, patho_sp, hostGeneNumbb, pathoGeneNumb,
            patoPerHostGeneration, numOfHostGenerations, hostMutationProb,
            pathoMutationProb, HeteroHomo, deletion, duplication, maxGene, 
            alpha, fixedAntigPosit);
    
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
    std::cout << "Host population all set!" << std::endl;
    ENV.setPathoPopulationFourClades(pathoPopSize, antigenLength, pathoGeneNumb,
                                     patho_sp, mhcGeneLength, 0, fixedAntigPosit);
    std::cout << "Pathogen population all set!" << std::endl;
    hostMutationProb = ENV.MMtoPMscaling(hostMutationProb, mhcGeneLength);
    Data2file.savePathoNoMuttList(ENV);
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
        Data2file.savePathoPopulToFile(ENV, 0);
        Data2file.saveHostPopulToFile(ENV, 0);
        Data2file.saveHostGeneticDivers(ENV, 0);
        Data2file.saveHostGeneNumbers(ENV, 0);
        for(int i = 1; i <= numOfHostGenerations; ++i){
            for(int j = 0; j < patoPerHostGeneration; ++j){
                ENV.infectOneFromOneSpecHetero();
                ENV.selectAndReproducePathoFixedPopSizes();
                ENV.mutatePathogensWithRestric(pathoMutationProb, mhcGeneLength, i);
                ENV.clearPathoInfectionData();
            }
            ENV.calculateHostsFitnessExpScalingUniqAlleles(alpha);
            ENV.selectAndReprodHostsReplace();
            ENV.mutateHostsWithDelDupl(hostMutationProb, deletion, duplication,
                    maxGene, i);
            Data2file.saveHostGeneticDivers(ENV, i);
            Data2file.saveHostGeneNumbers(ENV, i);
            ENV.clearHostInfectionsData();
//           std::cout << "Host loop " << i << " finished" << std::endl;
        }
        ENV.infectOneFromOneSpecHetero();
    } else {
       std::cout << "This instance of the model allows only heterozygote";
       std::cout << " advantage. Sorry :-(" << std::endl; 
       return 0;
    }

//    ENV.selectAndReproducePathoFixedPopSizes();
//    ENV.mutatePathogensWithRestric(pathoMutationProb, mhcGeneLength, numOfHostGenerations);
//    ENV.calculateHostsFitnessExpScalingUniqAlleles(alpha);
//    ENV.selectAndReprodHostsReplace();
//    ENV.mutateHostsWithDelDupl(hostMutationProb, deletion, duplication, 
//            maxGene, numOfHostGenerations);
//    Data2file.saveHostGeneticDivers(ENV, numOfHostGenerations);
//    Data2file.saveHostGeneNumbers(ENV, numOfHostGenerations);
    Data2file.savePathoPopulToFile(ENV, numOfHostGenerations);
    Data2file.saveHostPopulToFile(ENV, numOfHostGenerations);
       
    std::cout << "Run finished. Check the output files for results." << std::endl;
    
    std::cout << std::endl;
    
    return 0;
}
