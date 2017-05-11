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

#include "RandomNumbs.h"
#include "Environment.h"
#include "H2Pinteraction.h"

typedef std::string sttr;
typedef boost::dynamic_bitset<> antigenstring;

Environment::Environment() {
}

//Environment::Environment(const Environment& orig) {
//}

Environment::~Environment() {
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
        double fixedAntigenFrac){
    std::set<unsigned long> NoMutSet;
    RandomNumbs * p_RandomNumbs = RandomNumbs::getInstance();
    for(int i = 0; i < numb_of_species; ++i){
        NoMutSet.clear();
        NoMutsVec.push_back(NoMutSet);
        for(unsigned long j = 0; j < antigen_size; ++j){
            if(p_RandomNumbs->NextReal(0.0, 1.0) <= fixedAntigenFrac){
                NoMutsVec.back().insert(j);
            }
        }
    }
}

/**
 * @brief Core method. It defines "no mutation sites" in 4-bit-long packages of
 * the antigen for all  individual pathogen species in the simulation. It should
 * be run only ones per simulation.
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
void Environment::setNoMutsVecInFours(int numb_of_species, int antigen_size,
        double fixedAntigenFrac){
    std::set<unsigned long> NoMutSet;
    RandomNumbs * p_RandomNumbs = RandomNumbs::getInstance();
    fixedAntigenFrac = fixedAntigenFrac / 4.0;
    for(int i = 0; i < numb_of_species; ++i){
        NoMutSet.clear();
        NoMutsVec.push_back(NoMutSet);
        unsigned long indexCounter = 0;
        while (indexCounter < antigen_size) {
            if(p_RandomNumbs->NextReal(0.0, 1.0) <= fixedAntigenFrac){
                NoMutsVec.back().insert(indexCounter);
                NoMutsVec.back().insert(indexCounter + 1);
                NoMutsVec.back().insert(indexCounter + 2);
                NoMutsVec.back().insert(indexCounter + 3);
                indexCounter = indexCounter + 3;
            } else {
                indexCounter++;
            }
        }
    }
}


/**
 *  @brief Core method. It defines "no mutation sites" of the antigen for all
 * individual pathogen species in the simulation. It should be run only ones per
 * simulation.
 *
 * Creates a vector of <a href="http://www.cplusplus.com/reference/set/set/">
 * STD sets</a> containing indices of antigen bits in which  mutations are not
 * allowed. The length of the vector equals the number of pathogen species. Each
 * species has its own vector. Number of fixed sites is a user-defined parameter.
 * This version generates only four unique "no-mutation" vectors (all vectors are
 * just copies of one of the four) to differentiate "clads" of pathogen species.
 *
 * @param numb_of_species - number of pathogen species
 * @param antigen_size - number of bits per antigen
 * @param fixedAntigenFrac - fraction of bits in antigens which need to be fixed
 */
void Environment::setNoMutsVecFourClads(int numb_of_species, unsigned long antigen_size,
        double fixedAntigenFrac){
    std::set<unsigned long> NoMutSet;
    std::vector<std::set<unsigned long>> TmpNoMutsVec;
    RandomNumbs * p_RandomNumbs = RandomNumbs::getInstance();
    for(int i = 0; i < 4; ++i){
        NoMutSet.clear();
        TmpNoMutsVec.push_back(NoMutSet);
        for(unsigned long j = 0; j < antigen_size; ++j){
            if(p_RandomNumbs->NextReal(0.0, 1.0) <= fixedAntigenFrac){
                TmpNoMutsVec.back().insert(j);
            }
        }
    }
    int spp_count;
    for(int i = 0; i < numb_of_species; ++i){
        spp_count = i % 4;
        NoMutsVec.push_back(TmpNoMutsVec[spp_count]);
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
 */
void Environment::setHostRandomPopulation(int pop_size, unsigned long gene_size,
                                          unsigned long chrom_size, int timeStamp){
    for(int i = 0; i < pop_size; ++i){
        HostPopulation.push_back(Host());
        HostPopulation.back().setNewHost(chrom_size, gene_size, timeStamp);
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
 */
void Environment::setHostRandomPopulation(int pop_size, unsigned long gene_size,
                                          unsigned long chrom_size_lower,
                                          unsigned long chrom_size_uper, int timeStamp){
    if(chrom_size_lower > chrom_size_uper){
        unsigned long tmp_size = chrom_size_lower;
        chrom_size_lower = chrom_size_uper;
        chrom_size_uper = tmp_size;
    }
    RandomNumbs * p_RandomNumbs = RandomNumbs::getInstance();
    for(int i = 0; i < pop_size; ++i){
        HostPopulation.push_back(Host());
        HostPopulation.back().setNewHost(p_RandomNumbs->NextLongInt(chrom_size_lower,
                chrom_size_uper), gene_size, timeStamp);
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
 */
void Environment::setHostClonalPopulation(int pop_size, unsigned long gene_size,
                                          unsigned long chrom_size, int timeStamp){
    std::vector<Host> tmpPopulation;
    tmpPopulation.push_back(Host());
    tmpPopulation.back().setNewHomozygHost(chrom_size, gene_size, timeStamp);
    //tmpPopulation.back().setNewHost(chrom_size, gene_size, timeStamp);
    for(int i = 0; i < pop_size; ++i){
        HostPopulation.push_back(tmpPopulation.back());
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
 */
void Environment::setPathoPopulatioUniformGenome(int pop_size, unsigned long antigenSize,
        int chrom_size, int numb_of_species, unsigned long mhcSize, int timeStamp,
        double fixedAntigenFrac){
    Environment::setNoMutsVector(numb_of_species, antigenSize, fixedAntigenFrac);
    if (numb_of_species > pop_size) numb_of_species = pop_size;
    int indiv_per_species = pop_size / numb_of_species;
    int indiv_left = pop_size % numb_of_species;
    std::vector<Pathogen> OneSpeciesVector;
    for (int i = 0; i < numb_of_species; ++i){
        for(int j = 0; j < indiv_per_species; ++j){
            OneSpeciesVector.push_back(Pathogen());
            if(i+1 != indiv_per_species){
                OneSpeciesVector.back().setNewPathogen(chrom_size, antigenSize,
                        mhcSize, i, timeStamp);
            } else {
                OneSpeciesVector.back().setNewPathogen(chrom_size, antigenSize,
                        mhcSize, i, timeStamp);
            }
            if(indiv_left){
                OneSpeciesVector.push_back(Pathogen());
                OneSpeciesVector.back().setNewPathogen(chrom_size, antigenSize,
                        mhcSize,  i, timeStamp);
                indiv_left--;
            }
        }
        PathPopulation.push_back(OneSpeciesVector);
        OneSpeciesVector.clear();
    }
}

/**
 * @brief Core method. Initializes the pathogen population.
 *
 * Given the number of individuals, number of bit per gene, desired number of
 * genes in a genome and desired number of pathogen species it generates
 * random population of pathogens. Number of individuals will be evenly
 * distributed between species and each species consists of identical clones.
 * Antigens in species are not assigned at random, but are spread evenly in all
 * possible bit strings space.
 *
 * @@param pop_size - total number of individuals
 * @param antigenSize - number of bits per gene
 * @param chrom_size - number of genes per genome
 * @param numb_of_species - number of pathogen species
 * @param mhcSize - number of bits in MHC protein
 * @param timeStamp - current time (number of the model iteration)
 * @param fixedAntigenFrac - fraction of bits in antigens which need to be fixed
 */
void Environment::setPathoPopulatioDistincSpp(int pop_size, unsigned long antigenSize,
        int chrom_size, int numb_of_species, unsigned long mhcSize, int timeStamp,
        double fixedAntigenFrac){
    Environment::setNoMutsVector(numb_of_species, antigenSize, fixedAntigenFrac);
    if (numb_of_species > pop_size) numb_of_species = pop_size;
    int indiv_per_species = pop_size / numb_of_species;
    int indiv_left = pop_size % numb_of_species;
    std::vector<Pathogen> PathoSppTemplateVector;
    PathoSppTemplateVector.push_back(Pathogen());
    antigenstring atniStrr;
    PathoSppTemplateVector.back().setNewPathogen(chrom_size, antigenSize,
                                                 mhcSize, 0, timeStamp);
    for(int kk = 1; kk < numb_of_species; ++kk){
        if(kk % 2){
            atniStrr = PathoSppTemplateVector.back().getSingleAntigen(0);
            PathoSppTemplateVector.push_back(Pathogen());
            PathoSppTemplateVector.back().setNewPathogenNthSwap(chrom_size, atniStrr,
            mhcSize, kk, timeStamp, kk);
        } else {
            atniStrr = PathoSppTemplateVector.back().getSingleAntigen(0);
            PathoSppTemplateVector.push_back(Pathogen());
            PathoSppTemplateVector.back().setNewPathogenNthSwap(chrom_size, atniStrr,
            mhcSize, kk, timeStamp, 1);
        }
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
 */
void Environment::setPathoPopulatioDivSpecies(int pop_size, unsigned long antigenSize,
        int chrom_size, int numb_of_species, unsigned long mhcSize, int timeStamp,
        double fixedAntigenFrac){
    Environment::setNoMutsVector(numb_of_species, antigenSize, fixedAntigenFrac);
    if (numb_of_species > pop_size) numb_of_species = pop_size;
    int indiv_per_species = pop_size / numb_of_species;
    int indiv_left = pop_size % numb_of_species;
    std::vector<Pathogen> PathoSppTemplateVector;
    for(int kk = 0; kk < numb_of_species; ++kk){
        PathoSppTemplateVector.push_back(Pathogen());
        PathoSppTemplateVector.back().setNewPathogen(chrom_size, antigenSize,
                        mhcSize, kk, timeStamp);
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
 * @brief Core method. Initializes the pathogen population.
 *
 * Given the number of individuals, number of bit per gene, desired number of
 * genes in a genome and desired number of pathogen species it generates
 * random population of pathogens. Number of individuals will be evenly
 * distributed between species, each species draws its genes from the same
 * pool of possible bit strings and there are only 4 possible "clads" - groups
 * of species initialized with the same antigens and same mutations restrictions.
 *
 * @param pop_size - total number of individuals
 * @param antigenSize - number of bits per gene
 * @param chrom_size - number of genes per genome
 * @param numb_of_species - number of pathogen species
 * @param mhcSize - number of bits in MHC protein
 * @param timeStamp - current time (number of the model iteration)
 * @param fixedAntigenFrac - fraction of bits in antigens which need to be fixed
 */
void Environment::setPathoPopulationFourClades(int pop_size, unsigned long antigenSize,
        int chrom_size, int numb_of_species, unsigned long mhcSize, int timeStamp,
        double fixedAntigenFrac){
    Environment::setNoMutsVecFourClads(numb_of_species, antigenSize, fixedAntigenFrac);
    if (numb_of_species > pop_size) numb_of_species = pop_size;
    int indiv_per_species = pop_size / numb_of_species;
    int indiv_left = pop_size % numb_of_species;
    std::vector<Pathogen> PathoSppTemplateVector;
    PathoSppTemplateVector.push_back(Pathogen());
    antigenstring atniStrr;
    PathoSppTemplateVector.back().setNewPathogen(chrom_size, antigenSize,
                                                 mhcSize, 0, timeStamp);
    for(int kk = 1; kk < 4; ++kk){
        if(kk % 2){
            atniStrr = PathoSppTemplateVector.back().getSingleAntigen(0);
            PathoSppTemplateVector.push_back(Pathogen());
            PathoSppTemplateVector.back().setNewPathogenNthSwap(chrom_size, atniStrr,
            mhcSize, kk, timeStamp, 2);
        } else {
            atniStrr = PathoSppTemplateVector.back().getSingleAntigen(0);
            PathoSppTemplateVector.push_back(Pathogen());
            PathoSppTemplateVector.back().setNewPathogenNthSwap(chrom_size, atniStrr,
            mhcSize, kk, timeStamp, 1);
        }
    }
    int spp_count;
    std::vector<Pathogen> OneSpeciesVector;
    for (int i = 0; i < numb_of_species; ++i){
        spp_count = i % 4;
        for(int j = 0; j < indiv_per_species; ++j){
            OneSpeciesVector.push_back(PathoSppTemplateVector[spp_count]);
            OneSpeciesVector.back().setNewSpeciesNumber(i);
            if(indiv_left){
                OneSpeciesVector.push_back(PathoSppTemplateVector[spp_count]);
                OneSpeciesVector.back().setNewSpeciesNumber(i);
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
 * @param simil_mesure - number of bits which have to be similar, to expose
 * a pathogen. It's passed to H2Pinteraction::doesInfected() method.
 */
void Environment::infectOneFromOneSpecHetero(){
    H2Pinteraction H2P;
    unsigned long j;
        RandomNumbs * p_RandomNumbs = RandomNumbs::getInstance();
    unsigned long HostPopulationSize = HostPopulation.size();
    for(unsigned long i = 0; i < HostPopulationSize; ++i){
        unsigned long PathPopulationSize = PathPopulation.size();
        for(unsigned long sp = 0; sp < PathPopulationSize; ++sp){
            if(PathPopulation[sp].size()){
                j = p_RandomNumbs->NextLongInt(0, PathPopulation[sp].size()-1);
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
    for(unsigned long i = 0; i < HostPopulationSize; ++i){
        HostPopulation[i].calculateFitnessForDrift();
    }
}

/**
 * @brief Core method. Iterates through the host population and calculates the
 * Fitness for each single individual by calling
 * Host::calculateHostsFitnessAlphaXsqr(), which uses one over the square on
 * number of genes as a fitness cost.
 */
void Environment::calculateHostsFitnessAlphaXsqr(double alpha){
    unsigned long HostPopulationSize = HostPopulation.size();
    for(unsigned long i = 0; i < HostPopulationSize; ++i){
        HostPopulation[i].calculateFitnessAlphaXSqr(alpha);
    }
}

/**
 * @brief Core method. Iterates through the host population and calculates the
 * Fitness for each single individual by calling
 * Host::calculateFitnessExpFunc(), which uses a Gaussian function to accommodate
 * the costs of having lots of genes.
 */
void Environment::calculateHostsFitnessExpScaling(double alpha){
    unsigned long HostPopulationSize = HostPopulation.size();
    for(unsigned long i = 0; i < HostPopulationSize; ++i){
        HostPopulation[i].calculateFitnessExpFunc(alpha);
    }
}

/**
 * @brief Core method. Iterates through the host population and calculates the
 * Fitness for each single individual by calling
 * Host::calculateFitnessExpFuncUniqAlleles(), which uses a Gaussian function
 * to accommodate the costs of having lots of unique MHC alleles in chromosomes.
 */
void Environment::calculateHostsFitnessExpScalingUniqAlleles(double alpha){
    unsigned long HostPopulationSize = HostPopulation.size();
    for(unsigned long i = 0; i < HostPopulationSize; ++i){
        HostPopulation[i].calculateFitnessExpFuncUniqAlleles(alpha);
    }
}

/**
 * @brief Core method. Forms the next generation of hosts using the fitness
 * proportionate selection method.
 *
 * The new generation takes place of individuals which got removed from
 * population due to too low fitness.
 *
 * Iterates through the host population selecting the next generation
 * for reproduction using
 * <a href="http://en.wikipedia.org/wiki/Fitness_proportionate_selection">
 * fitness proportionate selection method</a> (also known as the roulette wheel
 * selection).
 */
void Environment::selectAndReprodHostsAddOffspring(){
    unsigned long pop_size = HostPopulation.size();
    double sum_of_fit = 0;
    double rnd;
    for(unsigned i = 0; i < pop_size; ++i){
        sum_of_fit += HostPopulation[i].getFitness();
    }
    if(sum_of_fit == 0){
        return;
    }
    RandomNumbs * p_RandomNumbs = RandomNumbs::getInstance();
    int n = 0;
    aley_oop:
    while(n < pop_size){
        rnd = p_RandomNumbs->NextReal(0, sum_of_fit);
        for(unsigned k = 0; k < pop_size; ++k) {
            rnd = rnd - HostPopulation[k].getFitness();
            if(rnd <= 0){
               HostPopulation[k].SelectedForReproduction += 1;
               n += 1;
               goto aley_oop;
            }
        }
    }
    // Remove the unfit ones
    unsigned long PopSize = HostPopulation.size();
    for (unsigned long i = PopSize - 1; i >= 0; --i){
        if (HostPopulation[i].SelectedForReproduction == 0){
            HostPopulation.erase(HostPopulation.begin() + i);
        }
    }
    // Reproduce the fit ones
//    PopSize = HostPopulation.size();
//    for (int i = PopSize - 1; i >= 0; --i){
//        if (HostPopulation[i].SelectedForReproduction){
//            for(int k = 1; k < HostPopulation[i].SelectedForReproduction; ++k){
//                HostPopulation.push_back(HostPopulation[k]);
//            }
//        }
//    }
    // Random mating
    PopSize = HostPopulation.size();
    unsigned long emptySlots = pop_size - PopSize;
    unsigned long indx_father, indx_mother;
    for (unsigned long p = 0; p < emptySlots; ++p){
        indx_father = p_RandomNumbs->NextLongInt(0, PopSize);
        indx_mother = p_RandomNumbs->NextLongInt(0, PopSize);
        while (indx_father == indx_mother){
            indx_mother = p_RandomNumbs->NextLongInt(0, PopSize);
        }
        HostPopulation.push_back(HostPopulation[indx_mother]);
        HostPopulation.back().assignChromTwo(HostPopulation[indx_father].getChromosomeTwo());
        // Randomly swaps places of chromosomes to avoid situation when they
        // effectively become two separate populations.
        HostPopulation.back().swapChromosomes();
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
    RandomNumbs * p_RandomNumbs = RandomNumbs::getInstance();
    int n = 0;
    aley_oop:
    while(n < pop_size){
        rnd = p_RandomNumbs->NextReal(0, sum_of_fit);
        unsigned long HostPopulationSize = HostPopulation.size();
        for(unsigned long k = 0; k < HostPopulationSize; ++k) {
            rnd = rnd - HostPopulation[k].getFitness();
            if(rnd <= 0){
               HostPopulation[k].SelectedForReproduction += 1;
               NewHostsVec.push_back(HostPopulation[k]);
               n += 1;
               goto second_parent;
            }
        }
        second_parent:
        rnd = p_RandomNumbs->NextReal(0, sum_of_fit);
        HostPopulationSize = HostPopulation.size();
        for(unsigned long p = 0; p < HostPopulationSize; ++p) {
            rnd = rnd - HostPopulation[p].getFitness();
            if(rnd <= 0){
               HostPopulation[p].SelectedForReproduction += 1;
               NewHostsVec.back().assignChromTwo(HostPopulation[p].getChromosomeTwo());
               // Randomly swaps places of chromosomes to avoid situation when
               // they effectively become two separate populations.
               NewHostsVec.back().swapChromosomes();
               goto aley_oop;
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
 * @brief Core method. Forms the next generation of hosts using the fitness
 * proportionate selection method. It can adjust species population sizes in
 * proportion to its individuals fitness values keeping the total number of
 * pathogens fixed.
 *
 * Iterates through the pathogen population selecting the next generation
 * for reproduction using  fitness proportionate selection method (roulette
 * wheel selection).
 */
void Environment::selectAndReproducePathoFlexPopSizes(){
    int rnd;
    int tot_patho_pop_size = 0;
    int total_ifected = 0;
    int sum_of_fit = 0;
    unsigned long PathPopulationSize = PathPopulation.size();
    for (unsigned long i = 0; i < PathPopulationSize; ++i){
        tot_patho_pop_size = tot_patho_pop_size + (int) PathPopulation[i].size();
        for (int j = 0; j < PathPopulation[i].size(); ++j){
            total_ifected += PathPopulation[i][j].NumOfHostsInfected;
//            PathPopulation[i][j].NumOfHostsInfected += 1;
            // A trick to have more survivors
            sum_of_fit += PathPopulation[i][j].NumOfHostsInfected;
        }
    }
//    std::cout << "Total infected: " << total_ifected << std::endl;
    if (total_ifected == 0){
        return;
    }
    RandomNumbs * p_RandomNumbs = RandomNumbs::getInstance();
    int n = 0;
    aley_oop:
    while(n < tot_patho_pop_size){
        rnd = p_RandomNumbs->NextInt(0, sum_of_fit);
        PathPopulationSize = PathPopulation.size();
        for (int i = 0; i < PathPopulationSize; ++i){
            unsigned long PathPopulationIthSize = PathPopulation[i].size();
            for (int j = 0; j < PathPopulationIthSize; ++j){
                rnd = rnd - PathPopulation[i][j].NumOfHostsInfected;
                if(rnd <= 0){
                    PathPopulation[i][j].SelectedToReproduct += 1;
                    n += 1;
                    goto aley_oop;
                }
            }
        }
    }
    unsigned long PopSize;
    // Elimination of the unfit
    unsigned long PopOfPopsSize = PathPopulation.size();
    for (unsigned long i = PopOfPopsSize - 1; i >= 0; --i){
        PopSize = PathPopulation[i].size();
        for (unsigned long j = PopSize - 1; j >= 0; --j){
            if(PathPopulation[i][j].SelectedToReproduct  == 0){
                PathPopulation[i].erase(PathPopulation[i].begin() + j);
            }
        }
    }
    // Reproduction of the fit ones
    PopOfPopsSize = PathPopulation.size();
    for (unsigned long i = PopOfPopsSize - 1; i >= 0; --i){
        PopSize = PathPopulation[i].size();
        for (unsigned long j = PopSize - 1; j >= 0; j--){
            if(PathPopulation[i][j].SelectedToReproduct){
                for(int k = 1; k < PathPopulation[i][j].SelectedToReproduct; ++k){
                    PathPopulation[i].push_back(PathPopulation[i][j]);
                }
            }
        }
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
    RandomNumbs * p_RandomNumbs = RandomNumbs::getInstance();
    PathPopulationSize = PathPopulation.size();
    for (int k = 0; k < PathPopulationSize; ++k){
        TmpPathVec.clear();
        int n = 0;
        aley_oop:
        while(n < PopSizes[k]){
           rnd = p_RandomNumbs->NextInt(0, SpecTotInfected[k]);
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
        if(PathPopulation[k].size() == TmpPathVec.size()){
            PathPopulation[k].clear();
            PathPopulation[k] = TmpPathVec;
        }else{
            std::cout << "Error in selectAndReproducePathoFixedPopSizes(): " <<
                "Size mismatch between the new and the old population!" << std::endl;
            std::cout << "old pop: " << PathPopulation[k].size() <<
                " | new pop: " << TmpPathVec.size() << std::endl;
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
    for(unsigned long k = 0; k < HostPopulationSzie; ++k){
        HostPopulation[k].clearInfections();
    }
}

/**
 * @brief Core method. Iterates through the all genes of the host population
 * and performs mutations in genes with a given probability.
 *
 * @param mut_probabl - probability of a mutation in a single gene.
 * @param timeStamp - current time (number of the model iteration)
 */
void Environment::mutateHosts(double mut_probabl, int timeStamp){
    unsigned long HostPopulationSzie = HostPopulation.size();
    for(unsigned long k = 0; k < HostPopulationSzie; ++k){
        HostPopulation[k].chromoMutProcess(mut_probabl, timeStamp);
    }
}

/**
 * @brief Core method. Iterates through the all genes of the host population
 * and performs mutations in genes with a given probability. Also deletions and
 * duplications of genes.
 *
 * @param mut_probabl - probability of a mutation in a single gene.
 * @param del - mutation probability, probability a gene will be deleted
 * @param dupli - mutation probability, probability a gene will be duplicated
 * (and added at the end of the Chromosome vector)
 * @param maxGene - maximal allowed number of genes in a chromosome
 * @param timeStamp - current time (number of the model iteration)
 */
void Environment::mutateHostsWithDelDupl(double mut_probabl, double del,
        double dupl, unsigned long maxGene, int timeStamp){
    unsigned long HostPopulationSzie = HostPopulation.size();
    for(unsigned long k = 0; k < HostPopulationSzie; ++k){
        HostPopulation[k].chromoMutProcessWithDelDupl(mut_probabl, del, dupl,
                maxGene, timeStamp);
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
 */
void Environment::mutateHostsWithDelDuplPointMuts(double pm_mut_probabl,
        double del, double dupl, unsigned long maxGene, int timeStamp){
    unsigned long HostPopulationSzie = HostPopulation.size();
    for(int k = 0; k < HostPopulationSzie; ++k){
        HostPopulation[k].chromoMutProcessWithDelDuplPointMuts(pm_mut_probabl,
                del, dupl, maxGene, timeStamp);
    }
}

/**
 * @brief Core method. When given micro-recombination mutation probability it
 * returns point-mutation probability calculated the way that it the average
 * number of new MHC match in both mutation scenarios.
 *
 * It implements the way in Ejsmond MJ, Radwan J. 2015. *Red Queen drives
 * positive selection on Major Histocompatibility Complex genes (MHC)*
 * [not published yet] micro-recombination mutation is transformed to
 * point-mutation probability:
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
 */
void Environment::mutatePathogens(double mut_probabl, unsigned long mhcSize, int timeStamp){
    unsigned long PathPopulationSize = PathPopulation.size();
    for (int i = 0; i < PathPopulationSize; ++i){
        unsigned long PathPopulationIthSize = PathPopulation[i].size();
        for (int j = 0; j < PathPopulationIthSize; ++j){
            PathPopulation[i][j].chromoMutProcess(mut_probabl, mhcSize, timeStamp);
        }
    }
}

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
 */
void Environment::mutatePathogensWithRestric(double mut_probabl, unsigned long mhcSize,
        int timeStamp){
    if (PathPopulation.size() == NoMutsVec.size()){
        unsigned long PathPopulationSize = PathPopulation.size();
        for(unsigned long i = 0; i < PathPopulationSize; ++i){
            unsigned long PathPopulationIthSize = PathPopulation[i].size();
            for(unsigned long j = 0; j < PathPopulationIthSize; ++j){
                PathPopulation[i][j].chromoMutProcessWithRestric(mut_probabl,
                        mhcSize, timeStamp, NoMutsVec[i]);
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
 *  with strong negative preference towards MHC similarity between mates.
 * 
 * There are no sexes as the host species is assumed a hermaphrodite. Each 
 * individual checks out all other individuals in the population in a random
 * order and mates with the first individual that has 0 same MHCs. If
 * 0-similarity-rule will not provide enough offspring to match the old population
 * size, then 1 same MHC is accepted, after that 2 same MHC are OK, etc. until
 * the algorithm will recreate a population of the same size as the original one.
 * 
 */
void Environment::matingWithNoCommonMHCwholePop(){
    unsigned long popSize = HostPopulation.size();
    // Generating a vector of shuffled indices 
    std::vector<int> indxVec;
    for (int &&i = 0; i < popSize; ++i){
        indxVec.push_back((int) i);
    }
    long seed = std::chrono::system_clock::now().time_since_epoch().count();
    shuffle (indxVec.begin(), indxVec.end(), std::default_random_engine(seed));
    
    std::vector<Host> NewHostsVec;
    NewHostsVec.clear();
    // proper mating
    int similCount, geneIndCount;
    similCount = 0;
    while(NewHostsVec.size() < popSize){
        for (unsigned long j = 0; j < HostPopulation.size(); ++j){
            for (unsigned long k = 0; k < popSize; ++k){
                geneIndCount = 0;
                for (unsigned long l = 0; l < HostPopulation[j].getChromoOneSize(); ++l){
                    for (unsigned long m = 0; m < HostPopulation[indxVec[k]].getChromoOneSize(); ++m){
                        if(HostPopulation[j].getOneGeneFromOne(l) == 
                           HostPopulation[indxVec[k]].getOneGeneFromOne(m)){
                            geneIndCount++;
                        }
                    }
                }
                for (unsigned long l = 0; l < HostPopulation[j].getChromoOneSize(); ++l){
                    for (unsigned long m = 0; m < HostPopulation[indxVec[k]].getChromoTwoSize(); ++m){
                        if(HostPopulation[j].getOneGeneFromOne(l) == 
                           HostPopulation[indxVec[k]].getOneGeneFromTwo(m)){
                            geneIndCount++;
                        }
                    }
                }
                for (unsigned long l = 0; l < HostPopulation[j].getChromoTwoSize(); ++l){
                    for (unsigned long m = 0; m < HostPopulation[indxVec[k]].getChromoOneSize(); ++m){
                        if(HostPopulation[j].getOneGeneFromTwo(l) == 
                           HostPopulation[indxVec[k]].getOneGeneFromOne(m)){
                            geneIndCount++;
                        }
                    }
                }
                for (unsigned long l = 0; l < HostPopulation[j].getChromoTwoSize(); ++l){
                    for (unsigned long m = 0; m < HostPopulation[indxVec[k]].getChromoTwoSize(); ++m){
                        if(HostPopulation[j].getOneGeneFromTwo(l) == 
                           HostPopulation[indxVec[k]].getOneGeneFromTwo(m)){
                            geneIndCount++;
                        }
                    }
                }
                if(geneIndCount == similCount and NewHostsVec.size() < popSize){
                    // mate two hosts, set a new individual
                    NewHostsVec.push_back(HostPopulation[j]);
                    NewHostsVec.back().assignChromTwo(HostPopulation[indxVec[k]].getChromosomeTwo());
                    NewHostsVec.back().swapChromosomes();
                    break;
                }
            }
        }
        similCount++;
    }
    if (HostPopulation.size() == NewHostsVec.size()){
        HostPopulation.clear();
        HostPopulation = NewHostsVec;
    }else{
        std::cout << "Error in matingWithNoCommonMHC(): Size mismatch " <<
                "between the new and the old population!" << std::endl;
        std::cout << "old pop: " << HostPopulation.size() <<
                " | new pop: " << NewHostsVec.size()  << std::endl;
    }
}

/**
 * @brief Core method. Creates a new generation of hosts by sexual reproduction
 *  with weak negative preference towards MHC similarity between mates.
 * 
 * There are no sexes as the host species is assumed a hermaphrodite. Each 
 * individual checks out all other individuals in the population in a random
 * order and mates with the first individual that has at least one different
 * MHC. If one-similarity-rule will not provide enough offspring to match the
 * old population size, then partners are matched at random until the algorithm
 * will recreate a population of the same size as the original one.
 * 
 */
void Environment::matingWithOneDifferentMHCwholePop(){
    unsigned long popSize = HostPopulation.size();
    // Generating a vector of shuffled indices 
    std::vector<int> indxVec;
    for (unsigned long i = 0; i < popSize; ++i){
        indxVec.push_back((int) i);
    }
    long seed = std::chrono::system_clock::now().time_since_epoch().count();
    shuffle (indxVec.begin(), indxVec.end(), std::default_random_engine(seed));
    RandomNumbs * p_RandomNumbs = RandomNumbs::getInstance();
    
    std::vector<Host> NewHostsVec;
    NewHostsVec.clear();
    bool isThereMatch;
    // proper mating
    while(NewHostsVec.size() < popSize){
        for (unsigned long j = 0; j < HostPopulation.size(); ++j){
            for (unsigned long k = 0; k < popSize; ++k){
                isThereMatch = false;
                for (unsigned long l = 0; l < HostPopulation[j].getChromoOneSize(); ++l){
                    for (unsigned long m = 0; m < HostPopulation[indxVec[k]].getChromoOneSize(); ++m){
                        if(HostPopulation[j].getOneGeneFromOne(l) == 
                           HostPopulation[indxVec[k]].getOneGeneFromOne(m)){
                            isThereMatch = true;
                            goto got_a_match;
                        }
                    }
                }
                for (unsigned long l = 0; l < HostPopulation[j].getChromoOneSize(); ++l){
                    for (unsigned long m = 0; m < HostPopulation[indxVec[k]].getChromoTwoSize(); ++m){
                        if(HostPopulation[j].getOneGeneFromOne(l) == 
                           HostPopulation[indxVec[k]].getOneGeneFromTwo(m)){
                            isThereMatch = true;
                            goto got_a_match;
                        }
                    }
                }
                for (unsigned long l = 0; l < HostPopulation[j].getChromoTwoSize(); ++l){
                    for (unsigned long m = 0; m < HostPopulation[indxVec[k]].getChromoOneSize(); ++m){
                        if(HostPopulation[j].getOneGeneFromTwo(l) == 
                           HostPopulation[indxVec[k]].getOneGeneFromOne(m)){
                            isThereMatch = true;
                            goto got_a_match;
                        }
                    }
                }
                for (unsigned long l = 0; l < HostPopulation[j].getChromoTwoSize(); ++l){
                    for (unsigned long m = 0; m < HostPopulation[indxVec[k]].getChromoTwoSize(); ++m){
                        if(HostPopulation[j].getOneGeneFromTwo(l) == 
                           HostPopulation[indxVec[k]].getOneGeneFromTwo(m)){
                            isThereMatch = true;
                            goto got_a_match;
                        }
                    }
                }
                got_a_match:
                if(isThereMatch and NewHostsVec.size() < popSize){
                    // mate two hosts, set a new individual
                    NewHostsVec.push_back(HostPopulation[j]);
                    NewHostsVec.back().assignChromTwo(HostPopulation[indxVec[k]].getChromosomeTwo());
                    NewHostsVec.back().swapChromosomes();
                    break;
                }
            }
            if (!isThereMatch and NewHostsVec.size() < popSize){
                NewHostsVec.push_back(HostPopulation[j]);
                NewHostsVec.back().
                  assignChromTwo(HostPopulation[p_RandomNumbs->NextLongInt(0, popSize-1)].
                  getChromosomeTwo());
                NewHostsVec.back().swapChromosomes();
            }
        }
    }
    if (HostPopulation.size() == NewHostsVec.size()){
        HostPopulation.clear();
        HostPopulation = NewHostsVec;
    }else{
        std::cout << "Error in matingWithOneDifferentMHC(): Size mismatch " <<
                "between the new and the old population!" << std::endl;
        std::cout << "old pop: " << HostPopulation.size() <<
                " | new pop: " << NewHostsVec.size()  << std::endl;
    }
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
    RandomNumbs * p_RandomNumbs = RandomNumbs::getInstance();
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
    std::random_device rnd_device;
    // Specify the size of the mates vector
    std::vector<int> matesVec(matingPartnerNumber);
    // the mating procedure
    while(NewHostsVec.size() < popSize){
        i = p_RandomNumbs->NextLongInt(0, popSize-1);
        // Specify the engine and create unique vector listing mates from population
        std::mt19937 mersenne_engine(rnd_device());
        std::uniform_int_distribution<unsigned long> dist(0, popSize-1);
        auto genn = std::bind(dist, mersenne_engine);
        generate(begin(matesVec), end(matesVec), genn);
        // find the best mate out of N randomly chosen
        highScore = maxGenomeSize;
        theBestMatch = matesVec[0];
        for (auto mate : matesVec) {
//            std::cout << mate << " ";
            geneIndCount = 0;
            for (unsigned long l = 0; l < HostPopulation[i].getUniqueMhcSize(); ++l){
                for (unsigned long m = 0; m < HostPopulation[mate].getUniqueMhcSize(); ++m){
                    if(HostPopulation[i].getOneGeneFromUniqVect(l) ==
                       HostPopulation[mate].getOneGeneFromUniqVect(m)){
                        geneIndCount++;
                    }
                }
            }
            if(geneIndCount < highScore){
                highScore = geneIndCount;
                theBestMatch = mate;
            }
        }
        if (int(NewHostsVec.size()) < popSize) {
            // mate two hosts, set a new individual
            NewHostsVec.push_back(HostPopulation[i]);
            NewHostsVec.back().assignChromTwo(HostPopulation[theBestMatch].getChromosomeTwo());
            NewHostsVec.back().swapChromosomes();
        }
//        std::cout << std::endl;
    }
    if (HostPopulation.size() == NewHostsVec.size()){
        HostPopulation.clear();
        HostPopulation = NewHostsVec;
    }else{
        std::cout << "Error in matingWithNoCommonMHC(): Size mismatch " <<
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
    for(int i = 0; i < NoMutsVec.size(); ++i){
        for (unsigned long const& possit : NoMutsVec[i]){
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
 * @brief  Data harvesting method. Creates a vector of unique genes with all their stats.
 */
void Environment::setUniqueGenes() {
    UniqueGenes.clear();
    std::vector<Gene> AllGenes;
    AllGenes.clear();
    // Finding unique genes
    if(HostPopulation.size()){
        for(unsigned long i = 0; i < HostPopulation.size(); ++i){
//            std::cout << HostPopulation[i].getChromoOneSize() << " ";
            for(unsigned long j = 0; j < HostPopulation[i].getChromoOneSize(); ++j){
                AllGenes.push_back(HostPopulation[i].ChromosomeOne[j]);
            }
            for(unsigned long l = 0; l < HostPopulation[i].getChromoTwoSize(); ++l){
                AllGenes.push_back(HostPopulation[i].ChromosomeTwo[l]);
            }
        }
        unsigned long GeneCounter = AllGenes.size();
        bool IfCountedLyst[GeneCounter];
        for (unsigned long w = 0; w < GeneCounter; ++w) {
            IfCountedLyst[w] = false;
        }
        for (unsigned long i = 0; i < GeneCounter; ++i) {
            if (!IfCountedLyst[i]) {
                UniqueGenes.push_back(AllGenes[i]);
                for (unsigned long j = i + 1; j < GeneCounter; ++j) {
                    if (i != j && !IfCountedLyst[j] && AllGenes[i].GenesTag == AllGenes[j].GenesTag){
                        IfCountedLyst[j] = true;
                    }
                }
            }
        }
        for(unsigned long l = 0; l < UniqueGenes.size(); ++l){
            UniqueGenes[l].presentedPathos = 0;
            UniqueGenes[l].occurence = 0;
            UniqueGenes[l].numbOfHostsBearing = 0;
        }
        for (unsigned long k = 0; k < UniqueGenes.size(); ++k) {
            for(unsigned long i = 0; i < GeneCounter; ++i){
                if(UniqueGenes[k].getTheRealGene() == AllGenes[i].getTheRealGene()){
                    UniqueGenes[k].presentedPathos += AllGenes[i].presentedPathos;
                    UniqueGenes[k].occurence +=  1;
                }
            }
        }
        for(unsigned long l = 0; l < UniqueGenes.size(); ++l){
            for(unsigned long i = 0; i < HostPopulation.size(); ++i){
                for(unsigned long j = 0; j < HostPopulation[i].getChromoOneSize(); ++j){
                    if(UniqueGenes[l].GenesTag ==  HostPopulation[i].ChromosomeOne[j].GenesTag){
                        UniqueGenes[l].numbOfHostsBearing += 1;
                        goto geneCheckDone;
                    }
                }
                for(unsigned long k = 0; k < HostPopulation[i].getChromoTwoSize(); ++k){
                    if(UniqueGenes[l].GenesTag ==  HostPopulation[i].ChromosomeTwo[k].GenesTag){
                        UniqueGenes[l].numbOfHostsBearing += 1;
                        goto geneCheckDone;
                    }
                }
                geneCheckDone:
                continue;
            }
        }
    }
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

