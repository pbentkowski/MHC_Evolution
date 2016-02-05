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

#include "Host.h"
#include "Pathogen.h"

/**
 * @brief Core class. Stores and handles the environment object that is the 
 * place where all the things are happening. It contains host and pathogen
 * populations stored in separate STL Vectors.
 */
class Environment {
private:
    std::vector<Host> HostPopulation;
    std::vector<std::vector<Pathogen> > PathPopulation;
    std::vector<std::set<int>> NoMutsVec;
public:
    // === Core methods ===
    Environment();
//    Environment(const Environment& orig);
    virtual ~Environment();
    void setNoMutsVector(int numb_of_species, int antigen_size, double fixedAntigenFrac);
    void setNoMutsVecInFours(int numb_of_species, int antigen_size, double fixedAntigenFrac);
    void setHostPopulation(int pop_size, int gene_size, int chrom_size, int timeStamp);
    void setHostPopulation(int pop_size, int gene_size, int chrom_size_lower,
        int chrom_size_uper, int timeStamp);
    void setPathoPopulatioUniformGenome(int pop_size, int gene_size, 
        int chrom_size, int numb_of_species, int mhcSize, int timeStamp);
    void setPathoPopulatioDivSpecies(int pop_size, int gene_size, 
        int chrom_size, int numb_of_species, int mhcSize, int timeStamp);
    void infectOneFromOneSpecHetero();
    void infectEveryOne(int simil_mesure);
    void calculateHostsFitnessPerGene();
    void calculateHostsFitnessPlainPresent();
    void calculateHostsFitnessForDrift();
    void calculateHostsFitnessAlphaXsqr(double alpha);
    void calculateHostsFitnessExpScaling(double alpha);
    void calculateHostsFitnessExpScalingUniqAlleles(double alpha);
    void selectAndReprodHostsAddOffspring();
    void selectAndReprodHostsReplace();
    void selectAndReproducePathoFlexPopSizes();
    void selectAndReproducePathoFixedPopSizes();
    void clearHostInfectionsData();
    void clearPathoInfectionData();
    void mutatePathogens(double mut_probabl, int mhcSize, int timeStamp);
    void mutatePathogensWithRestric(double mut_probabl,  int mhcSize, int timeStamp);
    void mutateHosts(double mut_probabl, int timeStamp);
    void mutateHostsWithDelDupl(double mut_probabl, double del, double dupl, 
        unsigned int maxGene, int timeStamp);
    void mutateHostsWithDelDuplPointMuts(double mut_probabl, double del, 
        double dupl, unsigned int maxGene, int timeStamp);
    double MMtoPMscaling(double MM_prob_mut, int geneLength);

    // === Data harvesting methods ===
    unsigned getPathoNumOfSpecies();
    unsigned getPathoSpeciesPopSize(unsigned spec_numb);
    unsigned getHostsPopSize();
    std::string getHostsTags();
    std::string getPathoGenesToString(int i, int j);
    std::string getHostGenesToString(int i);
    std::string getFixedBitsInAntigens();
    int getSingleHostGenomeSize(int indx);
    int getSingleHostChromoOneSize(int indx);
    int getSingleHostChromoTwoSize(int indx);
    int getSingleHostRealGeneOne(int i, int j);
    int getSingleHostRealGeneTwo(int i, int j);
    double getHostFitness(int indx);
    
};

#endif	/* ENVIRONMENT_H */

