/*
 * File:   Host.h
 * Author: Piotr Bentkowski : bentkowski.piotr@gmail.com
 *
 * Created on 16 February 2015, 16:13
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
#ifndef HOST_H
#define	HOST_H

#include <cstdlib>
#include <string>
#include <vector>
#include "boost/dynamic_bitset.hpp"

#include "Gene.h"

typedef boost::dynamic_bitset<> genestring;
typedef std::vector<Gene>  chromovector;

/**
 * @brief Core class. Stores and handles a single host object. Each host
 * has multiple instances of Gene class objects and stores them in 2 vectors
 * called ChromosomeOne and ChromosomeTwo which simulate a diploid genome.
 */
class Host {
public:
    Host();
    //Host(const Host& orig);
    virtual ~Host();
    // === Core methods ===
    std::vector<int> PathoSpecInfecting;
    std::vector<int> PathogesPresented;
    unsigned NumOfPathogesInfecting;
    unsigned NumOfPathogesPresented;
    unsigned NumOfMhcAlleles;
    int SelectedForReproduction;
    int TimeOfRecombin;
    void setNewHost(int num_of_loci, int gene_size, int timeStamp);
    void setNewHomozygHost(int num_of_loci, int gene_size, int timeStamp);
    void chromoMutProcess(double mut_probabl, int timeStamp);
    void chromoMutProcessWithDelDupl(double mut_probabl, double del,
        double dupli, int maxGene, int timeStamp);
    void chromoMutProcessWithDelDuplPointMuts(double mut_probabl, double del,
        double dupli, int maxGene, int timeStamp);
    void chromoRecombination(double recomb_prob, int timeStamp);
    void clearInfections();
    chromovector doCrossAndMeiosis(double corssing_prob);
    chromovector getChromosomeOne();
    chromovector getChromosomeTwo();
    chromovector mergeChromosomes();
    unsigned getGenomeSize();
    unsigned getChromoOneSize();
    unsigned getChromoTwoSize();
    double getChromoOneUniqAlleles();
    double getChromoTwoUniqAlleles();
    void assignChromOne(chromovector One);
    void assignChromTwo(chromovector Two);
    genestring getSingleGeneFromOne(int indx);
    genestring getSingleGeneFromTwo(int indx);
    unsigned long getHostIndvTag();
    unsigned long getHostMotherTag();
    void setHostIndvTag(unsigned long theTag);
    void setHostMotherTag(unsigned long theTag);
    void swapChromosomes();
    void calculateFitnessJustInfection();
    void calculateFitnessAccChromSize();
    void calculateFitnessForDrift();
    void calculateFitnessAlphaXSqr(double alpha);
    void calculateFitnessExpFunc(double alpha);
    void calculateFitnessExpFuncUniqAlleles(double alpha);
    double getFitness();
    // === Data harvesting methods ===
    std::string stringChromosomes();
    unsigned long int getOneGeneFromOne(int indx);
    unsigned long int getOneGeneFromTwo(int indx);
private:
    // === Very core methods ===
    std::vector<Gene> ChromosomeOne;
    std::vector<Gene> ChromosomeTwo;
    double Fitness;
};

#endif	/* HOST_H */

