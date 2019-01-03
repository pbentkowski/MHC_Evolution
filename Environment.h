/* 
 * File:   Environment.h
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
#ifndef ENVIRONMENT_H
#define	ENVIRONMENT_H

#include <vector>
#include <list>
#include <string>

#include "Random.h"
#include "Tagging_system.h"
#include "Host.h"
#include "Pathogen.h"

/**
 * @brief Core class. Stores and handles the environment object that is the 
 * place where all the things are happening. It contains host and pathogen
 * populations stored in separate STL Vectors. Handles all the population-wide
 * things (setting populations up, calculating fitness in all population, matching
 * pathogens and hosts for infections, mating hosts etc.). Also has methods for
 * fetching some of the stats.
 */
class Environment {
private:
    std::vector<Host> HostPopulation;
    std::vector<std::vector<Pathogen> > PathPopulation;
    std::vector<std::set<unsigned long>> NoMutsVec;
    Random* mRandGenArr;      //array of random generators: one for each thread
    unsigned int mRandGenArrSize;
public:
    // === Core methods ===
    explicit Environment(unsigned int numberOfThreads);
//    Environment(const Environment& orig);
    virtual ~Environment();
    void seedEnvsRNG();
    void setNoMutsVector(int numb_of_species, unsigned long antigen_size, double fixedAntigenFrac);
//    void setNoMutsVecInFours(int numb_of_species, int antigen_size, double fixedAntigenFrac);
//    void setNoMutsVecFourClads(int numb_of_species, unsigned long antigen_size, double fixedAntigenFrac);
    void setHostRandomPopulation(int pop_size, unsigned long gene_size, unsigned long chrom_size, int timeStamp,
                                 Tagging_system &tag);
    void setHostRandomPopulation(int pop_size, unsigned long gene_size, unsigned long chrom_size_lower,
                                 unsigned long chrom_size_uper, int timeStamp, Tagging_system &tag);
    void setHostClonalPopulation(int pop_size, unsigned long gene_size, unsigned long chrom_size, int timeStamp,
                                 Tagging_system &tag);
    void setPathoPopulatioUniformGenome(int pop_size, unsigned long gene_size,
        int chrom_size, int numb_of_species, unsigned long mhcSize, int timeStamp,
        double fixedAntigenFrac, Tagging_system &tag);
    void setPathoPopulatioDivSpecies(int pop_size, unsigned long gene_size,
        int numb_of_species, unsigned long mhcSize, int timeStamp, double fixedAntigenFrac,
         Tagging_system &tag);
    void infectOneFromOneSpecHetero();
    //void infectEveryOne(int simil_mesure);
    void calculateHostsFitnessPerGene();
    void calculateHostsFitnessPlainPresent();
    void calculateHostsFitnessForDrift();
    void calculateHostsFitnessAlphaXsqr(double alpha);
    void calculateHostsFitnessExpScaling(double alpha);
    void calculateHostsFitnessExpScalingUniqAlleles(double alpha);
    void selectAndReprodHostsReplace();
    void selectAndReprodHostsNoMating();
    void selectAndReproducePathoFixedPopSizes();
    void clearHostInfectionsData();
    void clearPathoInfectionData();
    void mutatePathogens(double mut_probabl, unsigned long mhcSize, int timeStamp);
    void mutatePathogensWithRestric(double mut_probabl,  unsigned long mhcSize, int timeStamp,
                                     Tagging_system &tag);
    void mutateHostsWithDelDuplPointMuts(double mut_probabl, double del, 
        double dupl, unsigned long maxGene, int timeStamp, Tagging_system &tag);
    double MMtoPMscaling(double MM_prob_mut, unsigned long geneLength);

    void matingWithNoCommonMHCsmallSubset(unsigned long matingPartnerNumber);
    void matingWithOneDifferentMHCsmallSubset(int matingPartnerNumber);
    void matingMeanOptimalNumberMHCsmallSubset(int matingPartnerNumber);
    void matingMaxDifferentMHCs(int matingPartnerNumber);
    void matingRandom();

    // === Data harvesting methods ===
    unsigned long getPathoNumOfSpecies();
    unsigned long getPathoSpeciesPopSize(unsigned long spec_numb);
    unsigned long getHostsPopSize();
    std::string getHostsTags();
    std::string getPathoGenesToString(unsigned long i, unsigned long j);
    std::string getHostGenesToString(int i);
    std::string getHostUniqMHCtoString(int i);
    std::string getFixedBitsInAntigens();
    std::string getNumbersOfPathogensPresented();
    std::string getNumbersOfMhcInMother();
    std::string getNumbersOfMhcInFather();
    std::string getNumbersOfUniqueMHCs();
    unsigned long getSingleHostGenomeSize(unsigned long indx);
    unsigned long getSingleHostChromoOneSize(unsigned long indx);
    unsigned long getSingleHostChromoTwoSize(unsigned long indx);
    unsigned long getSingleHostRealGeneOne(unsigned long i, unsigned long j);
    unsigned long getSingleHostRealGeneTwo(unsigned long i, unsigned long j);
    double getHostFitness(unsigned long indx);
    
};

#endif	/* ENVIRONMENT_H */

