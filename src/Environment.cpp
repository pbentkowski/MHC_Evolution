/*
 * File:   Environment.cpp
 * Author: Piotr Bentkowski : bentkowski.piotr@gmail.com
 *
 * Created on 18 February 2015, 16:46
 *
 *    This program is free software; you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation; either version 2 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *    MA 02110-1301, USA.
 */
#include <complex>
#include <vector>

#include <iostream>     // std::cout
#include <algorithm>    // std::shuffle
#include <chrono>       // std::chrono::system_clock
#include <iterator>
#include <thread>
#include <functional>   // std::bind

#include "Environment.h"
#include "H2Pinteraction.h"

typedef std::string sttr;
typedef boost::dynamic_bitset<> antigenstring;

/**
 * @brief Core method. Sets the number of threads used later for computing the model.
 *
 * @param numberOfThreads - number of threads you want to use for computing this mode.
 */
Environment::Environment(unsigned int numberOfThreads) {
    if(numberOfThreads == 0)
    {
        numberOfThreads = std::thread::hardware_concurrency();
        if(numberOfThreads == 0) //if the value is not well defined or not computable, set at least 1 thread
            numberOfThreads = 1;
    }
    omp_set_num_threads(numberOfThreads);
    mRandGenArrSize = numberOfThreads;
    mRandGenArr = new Random[mRandGenArrSize];
    seedEnvsRNG();
}

//Environment::Environment(const Environment& orig) {
//}

Environment::~Environment() = default;

void Environment::seedEnvsRNG() {
    for(unsigned int i = 0; i < mRandGenArrSize; ++i)
        mRandGenArr[i].reseed(std::random_device()());
}

/**
 * @brief Core method. It defines "no mutation sites" of the antigen for all
 * individual pathogen species in the simulation. It should be run only ones per
 * simulation.
 *
 * Creates a vector of <a href="http://www.cplusplus.com/reference/set/set/">
 * STD sets</a> containing indices of antigen bits in which  mutations are not
 * allowed. The length of the vector equals the number of pathogen species. Each
 * species has its own vector. Number of fixed sites is a user-defined parameter.
 *
 * @param numb_of_species - number of pathogen species
 * @param antigen_size - number of bits per antigen
 * @param fixedAntigenFrac - fraction of bits in antigens which need to be fixed
 */
void Environment::setNoMutsVector(int numb_of_species, unsigned long antigen_size,
        double fixedAntigenFrac)
{
    Random * rngGenPtr = mRandGenArr;
    std::set<unsigned long> NoMutSet;
//    #pragma omp parallel for default(none) shared(numb_of_species, rngGenPtr, antigen_size, fixedAntigenFrac) private(NoMutSet)
    for(int i = 0; i < numb_of_species; ++i){
        NoMutSet.clear();
        NoMutsVec.push_back(NoMutSet);
        for(unsigned long j = 0; j < antigen_size; ++j){
            if(rngGenPtr[omp_get_thread_num()].getUni() <= fixedAntigenFrac){
                NoMutsVec.back().insert(j);
            }
        }
    }
}


/**
 * @brief Core method. Initializes a vector containing host population.
 *
 * Sets and fills a vector containing the host population. Parameters like size
 * of the population, length of genes represented by bit-strings, number of
 * genes in a chromosome (remember that there are two chromosomes) are
 * user-defined, the rest is set at random (e.g. actual gene values).
 *
 * @param pop_size - number of host individuals in a simulation
 * @param gene_size - number of bits in bit-represented genes
 * @param chrom_size - number of genes in a chromosome
 * @param timeStamp - current time (number of the model iteration)
 * @param tag - pointer to the tagging system marking each gene variant
 */
void Environment::setHostRandomPopulation(int pop_size, unsigned long gene_size,
                                          unsigned long chrom_size, int timeStamp,
                                           Tagging_system &tag){
    for(int i = 0; i < pop_size; ++i) {
        HostPopulation.emplace_back(Host());
    }
    Random * rngGenPtr = mRandGenArr;
    #pragma omp parallel for default(none) shared(chrom_size, pop_size, rngGenPtr, gene_size, timeStamp, tag)
    for(int j = 0; j < pop_size; ++j){
        HostPopulation[j].setNewHost(chrom_size, gene_size, timeStamp, rngGenPtr[omp_get_thread_num()], tag);
    }
}

/**
 * @brief Core method. Initializes a vector containing host population.
 *
 * Sets and fills a vector containing the host population. Parameters like size
 * of the population and the length of genes represented by bit-strings are
 * user-defined. Whereas the range of the number of genes in a chromosome
 * (remember that there are two chromosomes) is selected randomly from a set
 * range of possible values, the rest of parameters is set totally at random
 * (e.g. actual gene values).
 *
 * @param pop_size - number of host individuals in a simulation
 * @param gene_size - number of bits in bit-represented genes
 * @param chrom_size_lower - lower limit of the number of genes in a chromosome
 * @param chrom_size_uper - upper limit of the number of genes in a chromosome
 * @param timeStamp - current time (number of the model iteration)
 * @param tag - pointer to the tagging system marking each gene variant
 */
void Environment::setHostRandomPopulation(int pop_size, unsigned long gene_size, unsigned long chrom_size_lower,
                                          unsigned long chrom_size_uper, int timeStamp,
                                           Tagging_system &tag){
    if(chrom_size_lower > chrom_size_uper){
        unsigned long tmp_size = chrom_size_lower;
        chrom_size_lower = chrom_size_uper;
        chrom_size_uper = tmp_size;
    }
    for(int i = 0; i < pop_size; ++i) {
        HostPopulation.emplace_back(Host());
    }
    Random * rngGenPtr = mRandGenArr;
    #pragma omp parallel for default(none) shared(pop_size, chrom_size_lower, chrom_size_uper, rngGenPtr, gene_size, timeStamp, tag)
    for(int j = 0; j < pop_size; ++j){
        HostPopulation[j].setNewHost(rngGenPtr[omp_get_thread_num()].getRandomFromUniform((unsigned int) chrom_size_lower,
                (unsigned int) chrom_size_uper), gene_size, timeStamp, rngGenPtr[omp_get_thread_num()], tag);
    }
}


/**
 * @brief Core method. Initializes a vector containing clonal host population.
 *
 * Sets and fills a vector containing the host population consisting of a single
 * homozygous clone. Parameters like size of the population, length of genes
 * represented by bit-strings, number of genes in a chromosome (remember that
 * there are two chromosomes) are user-defined, the rest is set at random
 * (e.g. actual gene values).
 *
 * @param pop_size - number of host individuals in a simulation
 * @param gene_size - number of bits in bit-represented genes
 * @param chrom_size - number of genes in a chromosome
 * @param timeStamp - current time (number of the model iteration)
 * @param tag - pointer to the tagging system marking each gene variant
 */
void Environment::setHostClonalPopulation(int pop_size, unsigned long gene_size,
                                          unsigned long chrom_size, int timeStamp,
                                           Tagging_system &tag){
    std::vector<Host> tmpPopulation;
    tmpPopulation.emplace_back(Host());
    Random * rngGenPtr = mRandGenArr;
    #pragma omp single
    {
        tmpPopulation.back().setNewHomozygHost(chrom_size, gene_size, timeStamp, rngGenPtr[omp_get_thread_num()], tag);
        //tmpPopulation.back().setNewHost(chrom_size, gene_size, timeStamp);
        for (int i = 0; i < pop_size; ++i) {
            HostPopulation.push_back(tmpPopulation.back());
        }
    }
}

/**
 * @brief Core method. Initializes the pathogen population.
 *
 * Given the number of individuals, number of bit per gene, desired number of
 * genes in a genome and desired number of pathogen species it generates
 * random population of pathogens. Number of individuals will be evenly
 * distributed between species and each species draws its genes from the same
 * pool of possible bit strings.
 *
 * @param pop_size - total number of individuals
 * @param antigenSize - number of bits per gene
 * @param chrom_size - number of genes per genome
 * @param numb_of_species - number of pathogen species
 * @param mhcSize - number of bits in MHC protein
 * @param timeStamp - current time (number of the model iteration)
 * @param fixedAntigenFrac - fraction of bits in antigens which need to be fixed
 * @param tag - pointer to the tagging system marking each gene variant
 *
void Environment::setPathoPopulatioUniformGenome(int pop_size, unsigned long antigenSize,
        int chrom_size, int numb_of_species, unsigned long mhcSize, int timeStamp,
        double fixedAntigenFrac, Tagging_system &tag){
    Environment::setNoMutsVector(numb_of_species, antigenSize, fixedAntigenFrac, rngGenPtr[omp_get_thread_num()]);
    if (numb_of_species > pop_size) numb_of_species = pop_size;
    int indiv_per_species = pop_size / numb_of_species;
    int indiv_left = pop_size % numb_of_species;
    std::vector<Pathogen> OneSpeciesVector;
    for (int i = 0; i < numb_of_species; ++i){
        for(int j = 0; j < indiv_per_species; ++j){
            OneSpeciesVector.push_back(Pathogen());
            if(i+1 != indiv_per_species){
                OneSpeciesVector.back().setNewPathogen(chrom_size, antigenSize,
                        mhcSize, i, timeStamp, randGen, tag);
            } else {
                OneSpeciesVector.back().setNewPathogen(chrom_size, antigenSize,
                        mhcSize, i, timeStamp, randGen, tag);
            }
            if(indiv_left){
                OneSpeciesVector.push_back(Pathogen());
                OneSpeciesVector.back().setNewPathogen(chrom_size, antigenSize,
                        mhcSize,  i, timeStamp, randGen, tag);
                indiv_left--;
            }
        }
        PathPopulation.push_back(OneSpeciesVector);
        OneSpeciesVector.clear();
    }
}
*/

/**
 * @brief Core method. Initializes the pathogen population.
 *
 * Given the number of individuals, number of bit per gene, desired number of
 * genes in a genome and desired number of pathogen species it generates
 * random population of pathogens. Number of individuals will be evenly
 * distributed between species and each species consists of identical clones.
 *
 * @param pop_size - total number of individuals
 * @param antigenSize - number of bits per gene
 * @param chrom_size - number of genes per genome
 * @param numb_of_species - number of pathogen species
 * @param mhcSize - number of bits in MHC protein
 * @param timeStamp - current time (number of the model iteration)
 * @param fixedAntigenFrac - fraction of bits in antigens which need to be fixed
 * @param tag - pointer to the tagging system marking each gene variant
 */
void Environment::setPathoPopulatioDivSpecies(int pop_size, unsigned long antigenSize,
        int numb_of_species, unsigned long mhcSize, int timeStamp,
        double fixedAntigenFrac, Tagging_system &tag)
{
    Environment::setNoMutsVector(numb_of_species, antigenSize, fixedAntigenFrac);
    if (numb_of_species > pop_size) numb_of_species = pop_size;
    int indiv_per_species = pop_size / numb_of_species;
    int indiv_left = pop_size % numb_of_species;
    std::vector<Pathogen> PathoSppTemplateVector;
    for(int kk = 0; kk < numb_of_species; ++kk) {
        PathoSppTemplateVector.emplace_back(Pathogen());
    }
    Random * rngGenPtr = mRandGenArr;
    #pragma omp parallel for default(none) \
        shared(numb_of_species, rngGenPtr, PathoSppTemplateVector, antigenSize, mhcSize, timeStamp, tag)
    for(int ll = 0; ll < numb_of_species; ++ll) {
        PathoSppTemplateVector[ll].setNewPathogen(antigenSize, mhcSize, ll, timeStamp,
                                                  rngGenPtr[omp_get_thread_num()], tag);
    }
    std::vector<Pathogen> OneSpeciesVector;
    for (int i = 0; i < numb_of_species; ++i){
        for(int j = 0; j < indiv_per_species; ++j){
            OneSpeciesVector.push_back(PathoSppTemplateVector[i]);
            if(indiv_left){
                OneSpeciesVector.push_back(PathoSppTemplateVector[i]);
                indiv_left--;
            }
        }
        PathPopulation.push_back(OneSpeciesVector);
        OneSpeciesVector.clear();
    }
}


/**
 * @brief Core method. Iterates through the host population and the parasite
 * population to "infect" the hosts with parasites. With heterozygote advantage
 * and it permits a pathogen species to infect a host only ONES.
 *
 * Each host is exposed to one, randomly selected individual from a pathogen
 * species. Exposition procedure for a single host is repeated for all
 * pathogen species. Maximum number of pathogens a host can contract in one go
 * equals to the number of pathogen species. Fitness is calculated with heterozygote
 * advantage added (antigen recognition just by one allele gives a full advantage).
 * One species can infect a host only ONES.
 *
 */
void Environment::infectOneFromOneSpecHetero(){
    H2Pinteraction H2P;
    unsigned long j;
    unsigned long HostPopulationSize = HostPopulation.size();
    Random * rngGenPtr = mRandGenArr;
    #pragma omp parallel for default(none) shared(rngGenPtr, HostPopulationSize) private(H2P, j)
    for(unsigned long i = 0; i < HostPopulationSize; ++i){
        unsigned long PathPopulationSize = PathPopulation.size();
        for(unsigned long sp = 0; sp < PathPopulationSize; ++sp){
            if(!PathPopulation[sp].empty()){
                j = rngGenPtr[omp_get_thread_num()].getRandomFromUniform(0, (unsigned int) PathPopulation[sp].size()-1);
                H2P.doesInfectedHeteroOnePerSpec(HostPopulation[i], PathPopulation[sp][j]);
            }
        }
    }
}

/**
 * @brief Core method. Iterates through the host population and calculates the
 * Fitness for each single individual by calling
 * Host::calculateFitnessAccChromSize(), which takes the number of MHC genes
 * under account.
 */
void Environment::calculateHostsFitnessPerGene(){
    unsigned long HostPopulationSize = HostPopulation.size();
    #pragma omp parallel for default(none) shared(HostPopulationSize)
    for(unsigned long i = 0; i < HostPopulationSize; ++i){
        HostPopulation[i].calculateFitnessAccChromSize();
    }
}

/**
 * @brief Core method. Iterates through the host population and calculates the
 * Fitness for each single individual by calling
 * Host::calculateHostsFitnessPlainInfect(), which is the plain-and-lame sum
 * of presented pathogens.
 */
void Environment::calculateHostsFitnessPlainPresent(){
    unsigned long HostPopulationSize = HostPopulation.size();
    #pragma omp parallel for default(none) shared(HostPopulationSize)
    for(unsigned long i = 0; i < HostPopulationSize; ++i){
        HostPopulation[i].calculateFitnessJustInfection();
    }
}

/**
 * @brief Core method. Iterates through the host population and calculates the
 * Fitness for each single individual by calling
 * Host::calculateFitnessForDrift(), which assigns "1" for each cell to make
 * the genetic driff work.
 */
void Environment::calculateHostsFitnessForDrift(){
    unsigned long HostPopulationSize = HostPopulation.size();
    #pragma omp parallel for default(none) shared(HostPopulationSize)
    for(unsigned long i = 0; i < HostPopulationSize; ++i){
        HostPopulation[i].calculateFitnessForDrift();
    }
}

/**
 * @brief Core method. Iterates through the host population and calculates the
 * Fitness for each single individual by calling
 * Host::calculateHostsFitnessAlphaXsqr(), which uses one over the square on
 * number of genes as a fitness cost.
 *
 * @param alpha - penalty for having too many MHC types factor for the host fitness function
 */
void Environment::calculateHostsFitnessAlphaXsqr(double alpha){
    unsigned long HostPopulationSize = HostPopulation.size();
    #pragma omp parallel for default(none) shared(HostPopulationSize, alpha)
    for(unsigned long i = 0; i < HostPopulationSize; ++i){
        HostPopulation[i].calculateFitnessAlphaXSqr(alpha);
    }
}

/**
 * @brief Core method. Iterates through the host population and calculates the
 * Fitness for each single individual by calling
 * Host::calculateFitnessExpFunc(), which uses a Gaussian function to accommodate
 * the costs of having lots of genes.
 *
 * @param alpha - penalty for having too many MHC types factor for the host fitness function
 */
void Environment::calculateHostsFitnessExpScaling(double alpha){
    unsigned long HostPopulationSize = HostPopulation.size();
    #pragma omp parallel for default(none) shared(HostPopulationSize, alpha)
    for(unsigned long i = 0; i < HostPopulationSize; ++i){
        HostPopulation[i].calculateFitnessExpFunc(alpha);
    }
}

/**
 * @brief Core method. Iterates through the host population and calculates the
 * Fitness for each single individual by calling
 * Host::calculateFitnessExpFuncUniqAlleles(), which uses a Gaussian function
 * to accommodate the costs of having lots of unique MHC alleles in chromosomes.
 *
 * @param alpha - penalty for having too many MHC types factor for the host fitness function
 */
void Environment::calculateHostsFitnessExpScalingUniqAlleles(double alpha){
    unsigned long HostPopulationSize = HostPopulation.size();
    #pragma omp parallel for default(none) shared(HostPopulationSize, alpha)
    for(unsigned long i = 0; i < HostPopulationSize; ++i){
        HostPopulation[i].calculateFitnessExpFuncUniqAlleles(alpha);
    }
}


/**
 * @brief Core method. Forms the next generation of hosts using the fitness
 * proportionate selection method. Replaces the old population with a new one.
 *
 * Iterates through the host population selecting the next generation
 * for reproduction using
 * <a href="http://en.wikipedia.org/wiki/Fitness_proportionate_selection">
 * fitness proportionate selection method</a> (also known as the roulette wheel
 * selection). Simulates random mating of hermaphrodites with no difference 
 * between sexes.
 */
void Environment::selectAndReprodHostsReplace(){
    std::vector<Host> NewHostsVec;
    NewHostsVec.clear();
    unsigned long pop_size = HostPopulation.size();
//    int tot_exposed = 0;
    double sum_of_fit = 0;
    double rnd;
    for(unsigned long i = 0; i < pop_size; ++i){
//        tot_exposed += HostPopulation[i].NumOfPathogesPresented;
//        HostPopulation[i].NumOfPathogesPresented += 1;
        sum_of_fit += HostPopulation[i].getFitness();
    }
    if(sum_of_fit == 0){
        return;
    }
    Random * rngGenPtr = mRandGenArr;
    #pragma omp single
    {
        int n = 0;
        aley_oop:
        while (n < pop_size) {
            rnd = rngGenPtr[omp_get_thread_num()].getRealDouble(0, sum_of_fit);
            unsigned long HostPopulationSize = HostPopulation.size();
            for (unsigned long k = 0; k < HostPopulationSize; ++k) {
                rnd = rnd - HostPopulation[k].getFitness();
                if (rnd <= 0) {
                    HostPopulation[k].SelectedForReproduction += 1;
                    NewHostsVec.push_back(HostPopulation[k]);
                    NewHostsVec.back().setMotherMhcNumber(HostPopulation[k].getUniqueMHCs().size());
                    n += 1;
                    goto second_parent;
                }
            }
            second_parent:
            rnd = rngGenPtr[omp_get_thread_num()].getRealDouble(0, sum_of_fit);
            HostPopulationSize = HostPopulation.size();
            for (unsigned long p = 0; p < HostPopulationSize; ++p) {
                rnd = rnd - HostPopulation[p].getFitness();
                if (rnd <= 0) {
                    HostPopulation[p].SelectedForReproduction += 1;
                    NewHostsVec.back().assignChromTwo(HostPopulation[p].getChromosomeTwo());
                    NewHostsVec.back().setFatherMhcNumber(HostPopulation[p].getUniqueMHCs().size());

                    // Randomly swaps places of chromosomes to avoid situation when
                    // they effectively become two separate populations.
                    NewHostsVec.back().swapChromosomes(rngGenPtr[omp_get_thread_num()]);
                    goto aley_oop;
                }
            }
        }
    }
    if (HostPopulation.size() == NewHostsVec.size()){
        HostPopulation.clear();
        HostPopulation = NewHostsVec;
    }else{
        std::cout << "Error in selectAndReprodHostsReplace(): Size mismatch " <<
                "between the new and the old population!" << std::endl;
        std::cout << "old pop: " << HostPopulation.size() <<
                " | new pop: " << NewHostsVec.size()  << std::endl;
    }
}


/**
 *@brief Core method. Forms the next generation of hosts using the fitness
 * proportionate selection method. Replaces the old population with a new one.
 *
 * Iterates through the host population selecting the next generation
 * for reproduction using
 * <a href="http://en.wikipedia.org/wiki/Fitness_proportionate_selection">
 * fitness proportionate selection method</a> (also known as the roulette wheel
 * selection). Successful individuals are simply cloned replacing the weak ones.
 */
void Environment::selectAndReprodHostsNoMating() {
    std::vector<Host> NewHostsVec;
    NewHostsVec.clear();
    unsigned long pop_size = HostPopulation.size();
    double sum_of_fit = 0;
    double rnd;
    for(unsigned long i = 0; i < pop_size; ++i){
        sum_of_fit += HostPopulation[i].getFitness();
    }
    if(sum_of_fit == 0){
        return;
    }
    Random * rngGenPtr = mRandGenArr;
    #pragma omp single
    {
        int n = 0;
        aley_oop:
        while (n < pop_size) {
            rnd = rngGenPtr[omp_get_thread_num()].getRealDouble(0, sum_of_fit);
            unsigned long HostPopulationSize = HostPopulation.size();
            for (unsigned long k = 0; k < HostPopulationSize; ++k) {
                rnd = rnd - HostPopulation[k].getFitness();
                if (rnd <= 0) {
                    HostPopulation[k].SelectedForReproduction += 1;
                    NewHostsVec.push_back(HostPopulation[k]);
                    n += 1;
                    goto aley_oop;
                }
            }
        }
    }
    if (HostPopulation.size() == NewHostsVec.size()){
        HostPopulation.clear();
        HostPopulation = NewHostsVec;
    }else{
        std::cout << "Error in selectAndReprodHostsNoMating(): Size mismatch " <<
                "between the new and the old population!" << std::endl;
        std::cout << "old pop: " << HostPopulation.size() <<
                " | new pop: " << NewHostsVec.size()  << std::endl;
    }
}


/**
 * @brief Core method. Forms the next generation of hosts using the fitness
 * proportionate selection method. It keeps population sizes of different
 * pathogens species at a fixed number.
 *
 * Iterates through the pathogen population selecting the next generation
 * for reproduction using fitness proportionate selection method (roulette
 * wheel selection).
 */
void Environment::selectAndReproducePathoFixedPopSizes(){
    std::vector<Pathogen> TmpPathVec;
    int rnd;
    unsigned long PopSizes[(int) PathPopulation.size()];
    int SpecTotInfected[(int) PathPopulation.size()];
    int total_ifected = 0;
    unsigned long PathPopulationSize = PathPopulation.size();
    for (unsigned long i = 0; i < PathPopulationSize; ++i){
        PopSizes[i] = PathPopulation[i].size();
        SpecTotInfected[i] = 0;
        for (int j = 0; j < PathPopulation[i].size(); ++j){
            SpecTotInfected[i] += PathPopulation[i][j].NumOfHostsInfected;
            total_ifected += PathPopulation[i][j].NumOfHostsInfected;
        }
    }
    if(total_ifected == 0)
        return;
    PathPopulationSize = PathPopulation.size();
    Random * rngGenPtr = mRandGenArr;
    #pragma omp parallel for ordered default(none) \
        shared(rngGenPtr, PathPopulationSize, PopSizes, SpecTotInfected, std::cout) private(rnd, TmpPathVec)
    for (int k = 0; k < PathPopulationSize; ++k){
        TmpPathVec.clear();
        int n = 0;
        aley_oop:
        while(n < PopSizes[k]){
           rnd = rngGenPtr[omp_get_thread_num()].getRandomFromUniform(0, (unsigned int) SpecTotInfected[k]);
           unsigned long PathPopulationKthSize = PathPopulation[k].size();
           for(unsigned long l = 0; l < PathPopulationKthSize; ++l){
               rnd = rnd - PathPopulation[k][l].NumOfHostsInfected;
               if(rnd <= 0){
                  PathPopulation[k][l].SelectedToReproduct += 1;
                  TmpPathVec.push_back(PathPopulation[k][l]);
                  n += 1;
                  goto aley_oop;
               }
           }
        }
        #pragma omp ordered
        {
            if (PathPopulation[k].size() == TmpPathVec.size()) {
                PathPopulation[k].clear();
                PathPopulation[k] = TmpPathVec;
            } else {
                std::cout << "Error in selectAndReproducePathoFixedPopSizes(): " <<
                          "Size mismatch between the new and the old population!" << std::endl;
                std::cout << "old pop: " << PathPopulation[k].size() <<
                          " | new pop: " << TmpPathVec.size() << std::endl;
            }
        }

    }
}

/**
 * @brief. Core method. Clears information about infection and fitness in the
 * whole pathogen population.
 */
void Environment::clearPathoInfectionData(){
    // Clear pathogens infection data
    unsigned long PathPopulationSize = PathPopulation.size();
    #pragma omp parallel for default(none) shared(PathPopulationSize)
    for (int i = 0; i < PathPopulationSize; ++i){
        unsigned long PathPopulationIthSize = PathPopulation[i].size();
        for (int j = 0; j < PathPopulationIthSize; ++j){
                PathPopulation[i][j].clearInfections();
            }
    }
}

/**
 * @brief. Core method. Clears information about infection and fitness in the
 * whole host population.
 */
void Environment::clearHostInfectionsData(){
    // Clear hosts infection data
    unsigned long HostPopulationSzie = HostPopulation.size();
    #pragma omp parallel for default(none) shared(HostPopulationSzie)
    for(unsigned long k = 0; k < HostPopulationSzie; ++k){
        HostPopulation[k].clearInfections();
    }
}


/**
 * @brief Core method. Iterates through the all genes of the host population
 * and performs POINT mutations in genes with a given probability. Also
 * deletions and duplications of genes.
 *
 * @param pm_mut_probabl - probability of a point mutation in a single gene.
 * @param del - mutation probability, probability a gene will be deleted
 * @param dupli - mutation probability, probability a gene will be duplicated
 * (and added at the end of the Chromosome vector)
 * @param timeStamp - current time (number of the model iteration)
 * @param maxGene - maximal allowed number of genes in a chromosome (user
 * defined parameter).
 * @param tag - pointer to the tagging system marking each gene variant
 */
void Environment::mutateHostsWithDelDuplPointMuts(double pm_mut_probabl,
        double del, double dupl, unsigned long maxGene, int timeStamp,
         Tagging_system &tag){
    unsigned long HostPopulationSzie = HostPopulation.size();
    Random * rngGenPtr = mRandGenArr;
    #pragma omp parallel for default(none) \
        shared(HostPopulationSzie, pm_mut_probabl, del, dupl, maxGene, timeStamp, rngGenPtr, tag)
    for(int k = 0; k < HostPopulationSzie; ++k){
        HostPopulation[k].chromoMutProcessWithDelDuplPointMuts(pm_mut_probabl,
                del, dupl, maxGene, timeStamp, rngGenPtr[omp_get_thread_num()], tag);
    }
}

/**
 *  @brief Core method. Iterates through the all genes of the host population
 * and performs all-MHC-mutation (writing a whole new bit string over an existing MHC)
 * of genes with a given probability. Also deletions and duplications of genes.
 *
 * @param mut_probabl - probability of mutation of the whole MHC flip.
 * @param del - mutation probability, probability a gene will be deleted
 * @param dupli - mutation probability, probability a gene will be duplicated
 * (and added at the end of the Chromosome vector)
 * @param timeStamp - current time (number of the model iteration)
 * @param maxGene - maximal allowed number of genes in a chromosome (user
 * defined parameter).
 * @param tag - pointer to the tagging system marking each gene variant
 */
void Environment::mutateHostWithDelDuplAllMHCchange(double mut_probabl, double del, double dupl, unsigned long maxGene,
                                                    int timeStamp, Tagging_system &tag) {
    unsigned long HostPopulationSzie = HostPopulation.size();
    Random * rngGenPtr = mRandGenArr;
    #pragma omp parallel for default(none) \
        shared(HostPopulationSzie, mut_probabl, del, dupl, maxGene, timeStamp, rngGenPtr, tag)
    for(int k = 0; k < HostPopulationSzie; ++k){
        HostPopulation[k].chromoMutProcessWithDelDupl(mut_probabl, del, dupl, maxGene,
                timeStamp, rngGenPtr[omp_get_thread_num()], tag);
    }
}


/**
 * @brief Core method. When given micro-recombination mutation probability it
 * returns point-mutation probability calculated the way that it the average
 * number of new MHC match in both mutation scenarios.
 *
 * It implements the way in Ejsmond MJ, Radwan J. 2015. *Red Queen drives
 * positive selection on Major Histocompatibility Complex genes (MHC)*
 * micro-recombination mutation is transformed to point-mutation probability:
 *
 *  \f$ p_{t} = 1 - \left[1 - p_{m} (1 - 0.5^{b})\right]^{1/b} \f$
 *
 * where \f$b\f$ is the number of bits, \f$p_{m}\f$ is the old-style mutation
 * probability. Only to avoid problems with numerical stability when calculating
 * the power of \f$1/b\f$ when \f$b\f$ is large and the base is small we
 * re-phrased it into:
 *
 * \f$ p_{t} = 1 - \exp\left( 1/b \cdot \ln [1 - p_{m} (1 - 0.5^{b})]\right) \f$
 *
 * @param MM_prob_mut - micro-recombination mutation ('old-style' mutation)
 * @param geneLength - number of bits in a gene representation
 * @return - calculated corresponding point mutation
 */
double Environment::MMtoPMscaling(double MM_prob_mut, unsigned long geneLength){
    double bitt = (double) geneLength;
    return 1.0 - std::exp((1.0 / bitt) * std::log(1.0 -  MM_prob_mut
            * (1.0 - std::pow(0.5, bitt))));
}

/**
 * @brief Core method. Iterates through the all genes of the pathogen population
 * and performs mutations in genes with a given probability.
 *
 * @param mut_probabl - probability of a mutation in a single gene.
 * @param mhcSize - number of bits in MHC protein
 * @param timeStamp - current time (number of the model iteration)
 *
void Environment::mutatePathogens(double mut_probabl, unsigned long mhcSize, int timeStamp){
    unsigned long PathPopulationSize = PathPopulation.size();
    for (int i = 0; i < PathPopulationSize; ++i){
        unsigned long PathPopulationIthSize = PathPopulation[i].size();
        for (int j = 0; j < PathPopulationIthSize; ++j){
            PathPopulation[i][j].chromoMutProcess(mut_probabl, mhcSize, timeStamp);
        }
    }
}
*/

/**
 * @brief Core method. Iterates through the all genes of the pathogen population
 * and performs mutations in genes with a given probability. But some positions
 * in the bit-string (gene) are not allowed to change.
 *
 * The "No-Mutation Vector" is defined within the Environment object.
 *
 * @param mut_probabl - probability of a mutation in a single gene.
 * @param mhcSize - number of bits in MHC protein
 * @param timeStamp - current time (number of the model iteration)
 * @param tag - pointer to the tagging system marking each gene variant
 */
void Environment::mutatePathogensWithRestric(double mut_probabl, unsigned long mhcSize,
        int timeStamp, Tagging_system &tag){
    if (PathPopulation.size() == NoMutsVec.size()){
        Random * rngGenPtr = mRandGenArr;
//    #pragma omp parallel for default(none) shared(
        unsigned long PathPopulationSize = PathPopulation.size();
        for(unsigned long i = 0; i < PathPopulationSize; ++i){
            unsigned long PathPopulationIthSize = PathPopulation[i].size();
            #pragma omp parallel for default(none) shared(rngGenPtr, tag, mut_probabl, mhcSize, timeStamp, i, PathPopulationIthSize)
            for(unsigned long j = 0; j < PathPopulationIthSize; ++j){
                PathPopulation[i][j].chromoMutProcessWithRestric(mut_probabl,
                        mhcSize, timeStamp, NoMutsVec[i], rngGenPtr[omp_get_thread_num()], tag);
            }
        }
    } else {
        std::cout << "Error in Environment::mutatePathogensWithRestric(): " <<
                "unequal number of species between PathPopulation and " <<
                "NoMutsVec" << std::endl;
    }
}

/**
 * @brief Core method. Gets the number of species of pathogens.
 *
 * @return number of pathogen species
 */
unsigned long Environment::getPathoNumOfSpecies(){
    return PathPopulation.size();
}

/**
 * @brief Core method. Gets a number of individuals in a selected species
 * of pathogen.
 *
 * @param spec_numb - number of selected species.
 * @return number of individuals of a selected species.
 */
unsigned long Environment::getPathoSpeciesPopSize(unsigned long spec_numb){
    return PathPopulation[spec_numb].size();
}

/**
 * @brief Core method. Gets the host population size.
 *
 * @return host population size
 */
unsigned long Environment::getHostsPopSize(){
    return HostPopulation.size();
}


/**
 * @brief Core method. Creates a new generation of hosts by sexual reproduction
 *  with strong negative preference towards MHC similarity between mates. Takes 
 * given number of possible sexual partners.
 * 
 * There are no sexes as the host species is assumed a hermaphrodite. Each 
 * individual checks out a user defined N number of random individuals from the 
 * population and mates with the individual that has the lowest number of similar
 * MHC genes. The process of random mating and selection is repeated until
 * the algorithm will recreate a population of the same size as the original one.
 * 
 * @param matingPartnerNumber - number of randomly selected partners an individual
 * will checks out eventually selecting one best to mate with.
 */
void Environment::matingWithNoCommonMHCsmallSubset(unsigned long matingPartnerNumber){
    unsigned long popSize = HostPopulation.size();
    std::vector<Host> NewHostsVec;
    NewHostsVec.clear();
    unsigned long i, maxGenomeSize, highScore, geneIndCount;
    int theBestMatch;
    maxGenomeSize = 0;
    for(auto indvidual : HostPopulation){
        if(indvidual.getGenomeSize() > maxGenomeSize){
            maxGenomeSize = indvidual.getUniqueMhcSize();
        }
    }
    // First create an instance of an random engine.
//    std::random_device rnd_device;
    // Specify the size of the mates vector
    std::vector<int> matesVec(matingPartnerNumber);
    // the mating procedure
    Random * rngGenPtr = mRandGenArr;
    #pragma omp single
    {
        while (NewHostsVec.size() < popSize) {
            i = rngGenPtr[omp_get_thread_num()].getRandomFromUniform(0, (unsigned int) popSize - 1);
            // Specify the engine and create unique vector listing mates from population
//            std::mt19937 mersenne_engine(rnd_device());
            std::mt19937 mersenne_engine = rngGenPtr[omp_get_thread_num()].returnEngene();
            std::uniform_int_distribution<unsigned long> dist(0, popSize - 1);
            auto genn = std::bind(dist, mersenne_engine);
            generate(begin(matesVec), end(matesVec), genn);
            // find the best mate out of N randomly chosen
            highScore = maxGenomeSize;
            theBestMatch = matesVec[0];
            for (auto mate : matesVec) {
//            std::cout << mate << " ";
                geneIndCount = 0;
                for (unsigned long l = 0; l < HostPopulation[i].getUniqueMhcSize(); ++l) {
                    for (unsigned long m = 0; m < HostPopulation[mate].getUniqueMhcSize(); ++m) {
                        if (HostPopulation[i].getOneGeneFromUniqVect(l) ==
                            HostPopulation[mate].getOneGeneFromUniqVect(m)) {
                            geneIndCount++;
                        }
                    }
                }
                if (geneIndCount < highScore) {
                    highScore = geneIndCount;
                    theBestMatch = mate;
                }
            }
            if (int(NewHostsVec.size()) < popSize) {
                // mate two hosts, set a new individual
                NewHostsVec.push_back(HostPopulation[i]);
                NewHostsVec.back().setMotherMhcNumber(HostPopulation[i].getUniqueMHCs().size());
                NewHostsVec.back().assignChromTwo(HostPopulation[theBestMatch].getChromosomeTwo());
                NewHostsVec.back().setFatherMhcNumber(HostPopulation[theBestMatch].getUniqueMHCs().size());
                NewHostsVec.back().swapChromosomes(rngGenPtr[omp_get_thread_num()]);
            }
        }
    }
    if (HostPopulation.size() == NewHostsVec.size()) {
        HostPopulation.clear();
        HostPopulation = NewHostsVec;
    } else {
        std::cout << "Error in matingWithNoCommonMHCsmallSubset(): Size mismatch " <<
                  "between the new and the old population!" << std::endl;
        std::cout << "old pop: " << HostPopulation.size() <<
                  " | new pop: " << NewHostsVec.size() << std::endl;
    }

//    // **** printing for testing  ****
//    for(auto indvidual : HostPopulation){
//        std::cout << "(" << indvidual.getMotherMhcNumber() << ", "
//                      << indvidual.getFatherMhcNumber() << ") ";
//    }
//    std::cout << std::endl;
}


/**
 * @brief Core method. Creates a new generation of hosts by sexual reproduction
 * with weak negative preference towards MHC similarity between mates.
 *
 * There are no sexes as the host species is assumed a hermaphrodite. Each
 * individual checks out a user defined N number of random individuals from the
 * population and mates with the first of the N individual that has at least one different
 * MHC gene not present in the choosing party genome. The process of random mating and
 * selection is repeated until the algorithm will recreate a population of the
 * same size as the original one.
 *
 * @param matingPartnerNumber - number of randomly selected partners an individual
 * will checks out eventually selecting one best to mate with.
 */
void  Environment::matingWithOneDifferentMHCsmallSubset(int matingPartnerNumber) {
    unsigned long popSize = HostPopulation.size();
    std::vector<Host> NewHostsVec;
    NewHostsVec.clear();
    unsigned long i, maxGenomeSize;
    int theBestMatch;
    bool differentGene;
    maxGenomeSize = 0;
    for(auto indvidual : HostPopulation){
        if(indvidual.getGenomeSize() > maxGenomeSize){
            maxGenomeSize = indvidual.getUniqueMhcSize();
        }
    }
    // First create an instance of an random engine.
//    std::random_device rnd_device;
    // Specify the size of the mates vector
    std::vector<int> matesVec((unsigned int) matingPartnerNumber);
    // the mating procedure
    Random * rngGenPtr = mRandGenArr;
    #pragma omp single
    {
        while (NewHostsVec.size() < popSize) {
            i = rngGenPtr[omp_get_thread_num()].getRandomFromUniform(0, (unsigned int) popSize - 1);
            // Specify the engine and create unique vector listing mates from population
//            std::mt19937 mersenne_engine(rnd_device());
            std::mt19937 mersenne_engine = rngGenPtr[omp_get_thread_num()].returnEngene();
            std::uniform_int_distribution<unsigned long> dist(0, popSize - 1);
            auto genn = std::bind(dist, mersenne_engine);
            generate(begin(matesVec), end(matesVec), genn);
            // find the best mate out of N randomly chosen
            theBestMatch = matesVec[0];
            for (auto mate : matesVec) {
//            std::cout << mate << " ";
                for (unsigned long l = 0; l < HostPopulation[i].getUniqueMhcSize(); ++l) {
                    differentGene = true;
                    for (unsigned long m = 0; m < HostPopulation[mate].getUniqueMhcSize(); ++m) {
                        if (HostPopulation[i].getOneGeneFromUniqVect(l) ==
                            HostPopulation[mate].getOneGeneFromUniqVect(m)) {
                            differentGene = false;
                            break;
                        }
                    }
                }
                if (differentGene) {
                    theBestMatch = mate;
                }
            }
            if (int(NewHostsVec.size()) < popSize) {
                // mate two hosts, set a new individual
                NewHostsVec.push_back(HostPopulation[i]);
                NewHostsVec.back().setMotherMhcNumber(HostPopulation[i].getUniqueMHCs().size());
                NewHostsVec.back().assignChromTwo(HostPopulation[theBestMatch].getChromosomeTwo());
                NewHostsVec.back().setFatherMhcNumber(HostPopulation[theBestMatch].getUniqueMHCs().size());
                NewHostsVec.back().swapChromosomes(rngGenPtr[omp_get_thread_num()]);
            }
//        std::cout << std::endl;
        }
    }
    if (HostPopulation.size() == NewHostsVec.size()){
        HostPopulation.clear();
        HostPopulation = NewHostsVec;
    }else{
        std::cout << "Error in matingWithOneDifferentMHCsmallSubset(): Size mismatch " <<
                  "between the new and the old population!" << std::endl;
        std::cout << "old pop: " << HostPopulation.size() <<
                  " | new pop: " << NewHostsVec.size()  << std::endl;
    }
}

/**
 * @brief Core method. Creates a new generation of hosts by sexual reproduction
 * picking a mate which ensures that the offspring will have the number of unique MHC
 * greater then their parents. Takes given number of possible sexual partners.
 *
 * @param matingPartnerNumber - number of randomly selected partners an individual
 * will checks out eventually selecting one best to mate with.
 */
void Environment::matingMeanOptimalNumberMHCsmallSubset(int matingPartnerNumber) {
    unsigned long popSize = HostPopulation.size();
    std::vector<Host> NewHostsVec;
    NewHostsVec.clear();
    unsigned long i;
    int theBestMatch;
    double sameGeneCount, highScore, uniqueMHCcount, score;
// First create an instance of an random engine.
//    std::random_device rnd_device;
    // Specify the size of the mates vector
    std::vector<int> matesVec(matingPartnerNumber);
    // the mating procedure
    Random * rngGenPtr = mRandGenArr;
    #pragma omp single
    {
        while (NewHostsVec.size() < popSize) {
            i = rngGenPtr[omp_get_thread_num()].getRandomFromUniform(0, popSize - 1);
            // Specify the engine and create unique vector listing mates from population
//        std::mt19937 mersenne_engine(rnd_device());
            std::mt19937 mersenne_engine = rngGenPtr[omp_get_thread_num()].returnEngene();
            std::uniform_int_distribution<unsigned long> dist(0, popSize - 1);
            auto genn = std::bind(dist, mersenne_engine);
            generate(begin(matesVec), end(matesVec), genn);
            // find the best mate out of N randomly chosen
            highScore = 0.0;
            theBestMatch = matesVec[0];
            for (auto mate : matesVec) {
//            std::cout << mate << " ";
                sameGeneCount = 0;
                uniqueMHCcount = (double) ( HostPopulation[i].getUniqueMhcSize()
                                            + HostPopulation[mate].getUniqueMhcSize());
                for (unsigned long l = 0; l < HostPopulation[i].getUniqueMhcSize(); ++l) {
                    for (unsigned long m = 0; m < HostPopulation[mate].getUniqueMhcSize(); ++m) {
                        if (HostPopulation[i].getOneGeneFromUniqVect(l) ==
                            HostPopulation[mate].getOneGeneFromUniqVect(m)) {
                            sameGeneCount++;
                        }
                    }
                }
                score = ( uniqueMHCcount - sameGeneCount ) / sameGeneCount ;
                if (score > highScore) {
                    highScore = score;
                    theBestMatch = mate;
                }
            }
            if (int(NewHostsVec.size()) < popSize) {
                // mate two hosts, set a new individual
                NewHostsVec.push_back(HostPopulation[i]);
                NewHostsVec.back().setMotherMhcNumber(HostPopulation[i].getUniqueMHCs().size());
                NewHostsVec.back().assignChromTwo(HostPopulation[theBestMatch].getChromosomeTwo());
                NewHostsVec.back().setFatherMhcNumber(HostPopulation[theBestMatch].getUniqueMHCs().size());
                NewHostsVec.back().swapChromosomes(rngGenPtr[omp_get_thread_num()]);
            }
//        std::cout << std::endl;
        }
    }
    if (HostPopulation.size() == NewHostsVec.size()){
        HostPopulation.clear();
        HostPopulation = NewHostsVec;
    }else{
        std::cout << "Error in matingMeanOptimalNumberMHCsmallSubset(): Size mismatch " <<
                "between the new and the old population!" << std::endl;
        std::cout << "old pop: " << HostPopulation.size() <<
                " | new pop: " << NewHostsVec.size()  << std::endl;
    }
}


/**
 *@brief Core method. Creates a new generation of hosts by sexual reproduction
 * picking a mate which has the most different MHCs then the selecting partner.
 * Takes a given number of possible sexual partners that are selected.
 *
 * @param matingPartnerNumber - number of randomly selected partners an individual
 * will checks out eventually selecting one that is the best to mate with.
 */
void Environment::matingMaxDifferentNumber(int matingPartnerNumber) {
    unsigned long popSize = HostPopulation.size();
    std::vector<Host> NewHostsVec;
    NewHostsVec.clear();
    unsigned long i;
    int theBestMatch;
    double sameGeneCount, highScore, score;
// First create an instance of an random engine.
//    std::random_device rnd_device;
    // Specify the size of the mates vector
    std::vector<int> matesVec((unsigned int) matingPartnerNumber);
    // the mating procedure
    Random * rngGenPtr = mRandGenArr;
    #pragma omp single
    {
        while (NewHostsVec.size() < popSize) {
            i = rngGenPtr[omp_get_thread_num()].getRandomFromUniform(0, (unsigned int) popSize - 1);
            // Specify the engine and create unique vector listing mates from population
//            std::mt19937 mersenne_engine(rnd_device());
            std::mt19937 mersenne_engine = rngGenPtr[omp_get_thread_num()].returnEngene();
            std::uniform_int_distribution<unsigned long> dist(0, popSize - 1);
            auto genn = std::bind(dist, mersenne_engine);
            generate(begin(matesVec), end(matesVec), genn);
            // find the best mate out of N randomly chosen
            highScore = 0.0;
            theBestMatch = matesVec[0];
            for (auto mate : matesVec) {
//            std::cout << mate << " ";
                sameGeneCount = 0;
                for (unsigned long l = 0; l < HostPopulation[i].getUniqueMhcSize(); ++l) {
                    for (unsigned long m = 0; m < HostPopulation[mate].getUniqueMhcSize(); ++m) {
                        if (HostPopulation[i].getOneGeneFromUniqVect(l) ==
                            HostPopulation[mate].getOneGeneFromUniqVect(m)) {
                            sameGeneCount++;
                        }
                    }
                }
                score = (double) HostPopulation[mate].getUniqueMhcSize() - sameGeneCount;
                if (score > highScore) {
                    highScore = score;
                    theBestMatch = mate;
                }
            }
            if (int(NewHostsVec.size()) < popSize) {
                // mate two hosts, set a new individual
                NewHostsVec.push_back(HostPopulation[i]);
                NewHostsVec.back().setMotherMhcNumber(HostPopulation[i].getUniqueMHCs().size());
                NewHostsVec.back().assignChromTwo(HostPopulation[theBestMatch].getChromosomeTwo());
                NewHostsVec.back().setFatherMhcNumber(HostPopulation[theBestMatch].getUniqueMHCs().size());
                NewHostsVec.back().swapChromosomes(rngGenPtr[omp_get_thread_num()]);
            }
//        std::cout << std::endl;
        }
    }
    if (HostPopulation.size() == NewHostsVec.size()){
        HostPopulation.clear();
        HostPopulation = NewHostsVec;
    }else{
        std::cout << "Error in matingMaxDifferentNumber(): Size mismatch " <<
                  "between the new and the old population!" << std::endl;
        std::cout << "old pop: " << HostPopulation.size() <<
                  " | new pop: " << NewHostsVec.size()  << std::endl;
    }

}


/**
 * @brief Core method. Creates a new generation of hosts by sexual reproduction
 * picking mates at random.
 *
 * Mates are randomly selected. It only makes sure that an individual does not mate with itself.
 *
 */
void Environment::matingRandom() {
    unsigned long popSize = HostPopulation.size();
    std::vector<Host> NewHostsVec;
    NewHostsVec.clear();
    unsigned long i;
    int theMatch;
    // the random mating procedure
    Random * rngGenPtr = mRandGenArr;
    #pragma omp single
    {
        while (NewHostsVec.size() < popSize) {
            i = 0;
            theMatch = 0;
            while (i == theMatch) {
                i = rngGenPtr[omp_get_thread_num()].getRandomFromUniform(0, (unsigned int) popSize - 1);
                theMatch = rngGenPtr[omp_get_thread_num()].getRandomFromUniform(0, (unsigned int) popSize - 1);
            }
            if (int(NewHostsVec.size()) < popSize) {
                // mate two hosts, set a new individual
                NewHostsVec.push_back(HostPopulation[i]);
                NewHostsVec.back().setMotherMhcNumber(HostPopulation[i].getUniqueMHCs().size());
                NewHostsVec.back().assignChromTwo(HostPopulation[theMatch].getChromosomeTwo());
                NewHostsVec.back().setFatherMhcNumber(HostPopulation[theMatch].getUniqueMHCs().size());
                NewHostsVec.back().swapChromosomes(rngGenPtr[omp_get_thread_num()]);
            }
        }
    }
    if (HostPopulation.size() == NewHostsVec.size()){
        HostPopulation.clear();
        HostPopulation = NewHostsVec;
    }else{
        std::cout << "Error in matingRandom(): Size mismatch " <<
                  "between the new and the old population!" << std::endl;
        std::cout << "old pop: " << HostPopulation.size() <<
                  " | new pop: " << NewHostsVec.size()  << std::endl;
    }
}

//==============================================================//

/**
 * @brief Data harvesting method. Gets the genome of selected pathogen in
 * human-readable format.
 *
 * @param i pathogen species number
 * @param j individual's index within the pathogen species vector
 * @return a string of the pathogen's chromosome in a human-readable format.
 */
std::string Environment::getPathoGenesToString(unsigned long i, unsigned long j){
    return PathPopulation[i][j].stringGenesFromGenome();
}

/**
 * @brief Data harvesting method. Gets the genome of selected host in
 * human-readable format.
 *
 * @param i individual's index within the host vector
 * @return  a string of the host's chromosome in a human-readable format.
 */
std::string Environment::getHostGenesToString(int i){
    return HostPopulation[i].stringChromosomes();
}


std::string Environment::getHostUniqMHCtoString(int i) {
    return HostPopulation[i].stringUniqMHCs();
}

/**
 * @brief Data harvesting method. Prepares a string with the list of fixed
 * "no mutation sites" in antigens used by the DataHarvester to write these data
 * to file.
 *
 * @return a string listing all the fixed sites in all pathogen species
 */
std::string Environment::getFixedBitsInAntigens(){
    sttr fixedMutStr;
    for (auto &i : NoMutsVec) {
        for (unsigned long const& possit : i){
            fixedMutStr += std::to_string(possit) + sttr(" ");
        }
        fixedMutStr += sttr("\n");
    }
    return fixedMutStr;
}

/**
 * @brief Data harvesting method. Prepares a string with the list of how many pathogen
 * species each host had presented.
 *
 * @return a string listing how many pathogens hosts presented
 */
std::string Environment::getNumbersOfPathogensPresented() {
    sttr presentedPatho;
    for(auto indvidual : HostPopulation){
            presentedPatho += sttr(" ") + std::to_string(indvidual.getNumberOfPresentedPatho());
    }
    presentedPatho +=  sttr("\n");
    return presentedPatho;
}

/**
 * @brief Data harvesting method. Prepares a string with the list of how many unique MHC
 * a "mother" had (a selecting party in mating procedures)
 *
 * @return a string listing how many MHCs a "mother" had
 */
std::string Environment::getNumbersOfMhcInMother() {
    sttr mhcsInMother;
    for(auto indvidual : HostPopulation){
        mhcsInMother += sttr(" ") + std::to_string(indvidual.getMotherMhcNumber());
    }
    mhcsInMother +=  sttr("\n");
    return mhcsInMother;
}

/**
 * @brief Data harvesting method. Prepares a string with the list of how many unique MHC
 * a "father" had (the party that is selected in mating procedures)
 *
 * @return a string listing how many MHCs a "father" had
 */
std::string Environment::getNumbersOfMhcInFather()  {
    sttr mhcsInFather;
    for(auto indvidual : HostPopulation){
        mhcsInFather += sttr(" ") + std::to_string(indvidual.getFatherMhcNumber());
    }
    mhcsInFather +=  sttr("\n");
    return mhcsInFather;
}


std::string Environment::getNumbersOfUniqueMHCs() {
    sttr uniqueMHCs;
    for(auto indvidual : HostPopulation){
        uniqueMHCs += sttr(" ") + std::to_string(indvidual.getNumbOfUniqMHCgenes());
    }
    uniqueMHCs +=  sttr("\n");
    return uniqueMHCs;
}

unsigned long Environment::getSingleHostGenomeSize(unsigned long indx){
    return HostPopulation[indx].getGenomeSize();
}


unsigned long Environment::getSingleHostChromoOneSize(unsigned long indx){
    return HostPopulation[indx].getChromoOneSize();
}


unsigned long Environment::getSingleHostChromoTwoSize(unsigned long indx){
    return HostPopulation[indx].getChromoTwoSize();
}


unsigned long Environment::getSingleHostRealGeneOne(unsigned long i, unsigned long j){
    return HostPopulation[i].getOneGeneFromOne(j);
}

unsigned long Environment::getSingleHostRealGeneTwo(unsigned long i, unsigned long j){
    return HostPopulation[i].getOneGeneFromTwo(j);
}

double Environment::getHostFitness(unsigned long indx){
    return (double) HostPopulation[indx].NumOfPathogesPresented;
}

