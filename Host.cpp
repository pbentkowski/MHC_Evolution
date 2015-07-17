/* 
 * File:   Host.cpp
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


#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include "boost/dynamic_bitset.hpp"

#include "Host.h"
#include "Gene.h"
#include "RandomNumbs.h"

typedef boost::dynamic_bitset<> genestring;
typedef std::vector<Gene> chromovector;
typedef std::string sttr;

Host::Host() {
}

//Host::Host(const Host& orig) {
//}

Host::~Host() {
}

/**
 * @brief Core method. Sets a new Host object and assigns a user-defined
 * number of genes of user-defined number of bits (as genes are bit-strings).
 * 
 * A host has two chromosome of equal length. 
 * 
 * @param num_of_loci - number of gene loci in a chromosome
 * @param gene_size - the length of the bit-string representing a gene
 */
void Host::setNewHost(int num_of_loci, int gene_size, int timeStamp){
    NumOfPathogesInfecting = 0;
    NumOfPathogesPresented = 0;
    SelectedForReproduction = 0;
    Fitness = 0.0;
    for(int i = 0; i < num_of_loci; ++i){
        ChromosomeOne.push_back(Gene());
        ChromosomeOne.back().setNewGene(gene_size, timeStamp);
        ChromosomeTwo.push_back(Gene());
        ChromosomeTwo.back().setNewGene(gene_size, timeStamp);
    }
}

/**
 * @brief Core method. Decides (on a random basis) if there will be any mutations
 * in the genome.
 * 
 * Iterates through the both chromosomes and calls gene mutation function. Genes 
 * mutate at random depending on the probability which was user-defined.
 * 
 * @param mut_probabl - mutation probability, a probability a gene will be
 * replaced by a new one
 * @param timeStamp - current time (number of the model iteration)
 */
void Host::chromoMutProcess(double mut_probabl, int timeStamp){
    for(int i = 0; i < ChromosomeOne.size(); ++i){
        ChromosomeOne[i].mutateGeneWhole(mut_probabl, timeStamp);
    }
    for(int i = 0; i < ChromosomeTwo.size(); ++i){
        ChromosomeTwo[i].mutateGeneWhole(mut_probabl, timeStamp);
    }
}

/**
 * @brief Core method. Decides (on a random basis) if there will be any mutations
 * in the genome.
 * 
 * Iterates through the both chromosomes and calls gene mutation function. Genes 
 * mutate at random depending on the probability which was user-defined. Also at
 * random a gene can be deleted or duplicated.
 * 
 * @param mut_probabl - mutation probability, a probability a gene will be
 * replaced by a new one
 * @param del - mutation probability, probability a gene will be deleted
 * @param dupli - mutation probability, probability a gene will be duplicated 
 * (and added at the end of the Chromosome vector)
 * @param timeStamp - current time (current number of the model iteration)
 */
void Host::chromoMutProcessWithDelDupl(double mut_probabl, double del, 
        double dupli, int maxGene, int timeStamp){
    RandomNumbs * p_RandomNumbs = RandomNumbs::getInstance();
    if(ChromosomeOne.size()){
        for(int i = ChromosomeOne.size() - 1; i >= 0; --i){
            ChromosomeOne[i].mutateGeneWhole(mut_probabl, timeStamp);
            if(ChromosomeOne.size() and ChromosomeOne.size() < maxGene 
                    and p_RandomNumbs->NextReal(0.0, 1.0) < dupli){
                ChromosomeOne.push_back(ChromosomeOne[i]);
            }
            if(ChromosomeOne.size() and p_RandomNumbs->NextReal(0.0, 1.0) < del){
                ChromosomeOne.erase(ChromosomeOne.begin() + i);
            }
        }
    }
    if(ChromosomeTwo.size()){
        for(int i = ChromosomeTwo.size() - 1; i >= 0; --i){
            ChromosomeTwo[i].mutateGeneWhole(mut_probabl, timeStamp);
            if(ChromosomeTwo.size() and ChromosomeTwo.size() < maxGene 
                    and p_RandomNumbs->NextReal(0.0, 1.0) < dupli){
                ChromosomeTwo.push_back(ChromosomeTwo[i]);
            }
            if(ChromosomeTwo.size() and p_RandomNumbs->NextReal(0.0, 1.0) < del){
                ChromosomeTwo.erase(ChromosomeTwo.begin() + i);
            }
        }
    }
}

/**
 * @brief Core method. Decides (on a random basis) if there will be any mutations
 * in the genome.
 * 
 * Iterates through the both chromosomes and calls gene mutation function. Point
 * sites in genes are mutated at random depending on the probability which was 
 * user-defined. Also at random a gene can be deleted or duplicated.
 * 
 * @param pm_mut_probabl - point mutation probability, a probability a gene 
 * will be replaced by a new one
 * @param del - mutation probability, probability a gene will be deleted
 * @param dupli - mutation probability, probability a gene will be duplicated 
 * (and added at the end of the Chromosome vector)
 * @param timeStamp - current time (current number of the model iteration)
 */
void Host::chromoMutProcessWithDelDuplPointMuts(double pm_mut_probabl,
        double del, double dupli, int maxGene, int timeStamp){
    RandomNumbs * p_RandomNumbs = RandomNumbs::getInstance();
    if(ChromosomeOne.size()){
        for(int i = ChromosomeOne.size() - 1; i >= 0; --i){
            ChromosomeOne[i].mutateGeneBitByBit(pm_mut_probabl, timeStamp);
            if(ChromosomeOne.size() and ChromosomeOne.size() < maxGene 
                    and p_RandomNumbs->NextReal(0.0, 1.0) < dupli){
                ChromosomeOne.push_back(ChromosomeOne[i]);
            }
            if(ChromosomeOne.size() and p_RandomNumbs->NextReal(0.0, 1.0) < del){
                ChromosomeOne.erase(ChromosomeOne.begin() + i);
            }
        }
    }
    if(ChromosomeTwo.size()){
        for(int i = ChromosomeTwo.size() - 1; i >= 0; --i){
            ChromosomeTwo[i].mutateGeneBitByBit(pm_mut_probabl, timeStamp);
            if(ChromosomeTwo.size() and ChromosomeTwo.size() < maxGene 
                    and p_RandomNumbs->NextReal(0.0, 1.0) < dupli){
                ChromosomeTwo.push_back(ChromosomeTwo[i]);
            }
            if(ChromosomeTwo.size() and p_RandomNumbs->NextReal(0.0, 1.0) < del){
                ChromosomeTwo.erase(ChromosomeTwo.begin() + i);
            }
        }
    }
}

/**
 * @brief Core method. Performing crossing over and chromosome selection before
 * mating.
 * 
 * It mixes up genes from two host chromosomes to form one child chromosome 
 * later pushed for mating. One can regulate which chromosome is preferred as
 * a source of genes by changing crossing-over probability parameter.
 * 
 * @param corssing_prob -probability of a gene crossing over; p=0.5 means that 
 * each chromosome gives half of its genes to the resulting meiotic chromosome.
 * p lesser then 0.5 means ChromosomeOne will give majority of genes; p larger
 * then 0.5 will favour ChromosomeTwo as a donor.
 * 
 * @return a Chromosome vector.
 */
 chromovector Host::doCrossAndMeiosis(double corssing_prob){
    RandomNumbs * p_RandomNumbs = RandomNumbs::getInstance();
    if(ChromosomeOne.size() == ChromosomeTwo.size()){
        chromovector new_chromos;
        if(p_RandomNumbs->NextReal(0.0, 1.0) < 0.5){
            corssing_prob = 1.0 - corssing_prob;
        }
        Gene tmpGene;
        for(int i = 0; i < ChromosomeOne.size(); ++i){
            if(p_RandomNumbs->NextReal(0.0, 1.0) < corssing_prob){
                new_chromos.push_back(ChromosomeOne[i]);
            } else {
               new_chromos.push_back(ChromosomeTwo[i]);
            }
        }
        return new_chromos;
    } else {
        std::cout << "Error in Host::doesItMutate(): chromosomes size"\
                " mismatch!" << std::endl;
        if(p_RandomNumbs->NextReal(0.0, 1.0) < 0.5){
            return ChromosomeTwo;
        }else{
            return ChromosomeOne;
        }
    }
}

/**
 * @brief Core method. Fetches a gene (as integer) with a given position on 
 * the Chromosome One.
 * 
 * @param indx - index of the gene you wanna get.
 * @return an integer representation of a gene with a given index.
 */
int Host::getOneGeneFromOne(int indx){
    if(ChromosomeOne.size()){
        if(indx < ChromosomeOne.size()){
            return ChromosomeOne[indx].getTheRealGene();
        }else{
            std::cout << "Error in Host::getOneGeneFromOne(): Index out of the "\
                    "range of the chromosome size. Fetching the last gene." << std::endl;
            std::cout << "Index: " << indx << ", Chromosome size: " << ChromosomeOne.size() << std::endl;
            return ChromosomeOne.back().getTheRealGene();
        }
    }else{
        std::cout << "Error in Host::getOneGeneFromOne(): Chromosome One "
                "has no genes! Returning 0." << std::endl;
        return -1;
    }
}

/**
 * @brief Core method. Fetches a gene (as integer) with a given position on 
 * the Chromosome Two.
 * 
 * @param indx - index of the gene you wanna get.
 * @return an integer representation of a gene with a given index.
 */
int Host::getOneGeneFromTwo(int indx){
    if(ChromosomeTwo.size()){
        if(indx < ChromosomeTwo.size()){
            return ChromosomeTwo[indx].getTheRealGene();
        }else{
            std::cout << "Error in Host::getOneGeneFromTwo(): Index out of the "\
                    "range of the chromosome size. Fetching the last gene." << std::endl;
            std::cout << "Index: " << indx << ", Chromosome size: " << ChromosomeTwo.size() << std::endl;
            return ChromosomeTwo.back().getTheRealGene();
        }
    }else{
        std::cout << "Error in Host::getOneGeneFromTwo(): Chromosome Two "
                "has no genes! Returning 0." << std::endl;
        return -1;
    }
}

/**
 * @brief Returns the total number of genes in the genome.
 * 
 * @return number of all MHC genes in the host
 */
unsigned Host::getGenomeSize(){
    return ChromosomeOne.size() + ChromosomeTwo.size();
}
 
/**
 * @brief Returns the number of genes in the Chromosome One.
 * 
 * @return number of MHC genes in the Chromosome One.
 */
unsigned Host::getChromoOneSize(){
    return ChromosomeOne.size();
}

/**
 * @brief Returns the number of genes in the Chromosome Two.
 * 
 * @return number of MHC genes in the Chromosome Two.
 */
unsigned Host::getChromoTwoSize(){
    return ChromosomeTwo.size();
}

/**
 * @brief Core method. Fetches a gene (bit string) with a given position on 
 * the chromosome One.
 * 
 * @param indx - index of the gene you wanna get.
 * @return a bitstring representation of a gene with a given index.
 */
genestring Host::getSingleGeneFromOne(int indx){
    if(indx < ChromosomeOne.size()){
        return ChromosomeOne[indx].getBitGene();
    }else{
        std::cout << "Error in Host::getSingleGene(): Index out of the range of"\
                " the chromosome size. Fetching the last gene." << std::endl;
        return ChromosomeOne.back().getBitGene();
    }
}

/**
 * @brief Core method. Fetches a gene (bit string) with a given position on 
 * the chromosome Two.
 * 
 * @param indx - index of the gene you wanna get.
 * @return a bitstring representation of a gene with a given index.
 */
genestring Host::getSingleGeneFromTwo(int indx){
    if(indx < ChromosomeTwo.size()){
        return ChromosomeTwo[indx].getBitGene();
    }else{
        std::cout << "Error in Host::getSingleGene(): Index out of the range of"\
                " the chromosome size. Fetching the last gene." << std::endl;
        return ChromosomeTwo.back().getBitGene();
    }
}

/**
 * @brief Core method. Returns the first chromosome.
 * 
 * @return Chromosome One
 */
chromovector Host::getChromosomeOne(){
    return ChromosomeOne;
}

/**
 * @brief Core method.  Returns the first chromosome.
 * 
 * @return Chromosome Two
 */
chromovector Host::getChromosomeTwo(){
    return ChromosomeTwo;
}

/**
 * @brief Data harvesting method. Merges both chromosomes to create an object
 * easier to handle is some situations.
 * 
 * @return Chromosome vector containing genes from both, Chromosome One and Two
 */
chromovector Host::mergeChromosomes(){
    std::vector<Gene> ChromoMerged;
    ChromoMerged.reserve(ChromosomeOne.size() + ChromosomeTwo.size() );
    ChromoMerged.insert(ChromoMerged.end(), ChromosomeOne.begin(), ChromosomeOne.end());
    ChromoMerged.insert(ChromoMerged.end(), ChromosomeTwo.begin(), ChromosomeTwo.end());
    return ChromoMerged;
}

/**
 * @brief Core method. Assigns a new chromosome to host's Chromosome ONE.
 * 
 * @param One - a STL vector of genes 
 */
void Host::assignChromOne(chromovector One){
    ChromosomeOne = One;
}

/**
 * @brief Core method. Assigns a new chromosome to host's Chromosome TWO.
 * 
 * @param One - a STL vector of genes 
 */
void Host::assignChromTwo(chromovector Two){
    ChromosomeTwo = Two;
}

/**
 * @brief Core method. Randomly swaps places of Chromosome One and Chromosome Two
 * to avoid situation when they effectively become two separate populations. 
 */
void Host::swapChromosomes(){
    RandomNumbs * p_RandomNumbs = RandomNumbs::getInstance();
    if(p_RandomNumbs->NextReal(0.0, 1.0) < 0.5){
        chromovector tmpChrome;
        tmpChrome = ChromosomeOne;
        ChromosomeOne = ChromosomeTwo;
        ChromosomeTwo = tmpChrome;
    }
}

/**
 * @brief Core method. Calculates host individual fitness as the number 
 * of exposed pathogens divided by the genome size.
 * 
 */
void Host::calculateFitnessAccChromSize(){
    if (ChromosomeOne.size() + ChromosomeTwo.size()){
        Fitness = (double) NumOfPathogesPresented / (double) (ChromosomeOne.size()
                + ChromosomeTwo.size() );
    } else {
        Fitness = 0.0;
    }
}

/**
 * @brief Core method. Calculates host individual fitness simply as the number 
 * of exposed pathogens.
 * 
 */
void Host::calculateFitnessJustInfection(){
    Fitness = (double) NumOfPathogesPresented;
}

/**
 * @brief Core method. Sets fitness as 1, a fixed value to simulate genetic drift.
 * 
 */
void Host::calculateFitnessForDrift(){
    Fitness = 1.0;
}

/**
 * @brief Core method. Calculates host individual fitness in proportion to one
 * over the square of the number o genes scaled to by factor \f$ \alpha \f$:
 * 
 * \f$ F = \frac{P}{(\alpha \cdot N)^{2}} \f$
 * 
 * where \f$0 < \alpha < 1 \f$, \f$ P \f$ is the number of pathogens exposed and
 * \f$ N \f$ in the sum of number of genes in both chromosomes.
 * 
 * @param alpha - scaling parameter
 */
void Host::calculateFitnessAlphaXSqr(double alpha){
    double NN = (double) (ChromosomeOne.size() + ChromosomeTwo.size());
    if (ChromosomeOne.size() + ChromosomeTwo.size()){
        Fitness = (double) NumOfPathogesPresented / std::pow(alpha * NN, 2.0);
    } else {
        Fitness = 0.0;
    }
}

/**
 * @brief Core method. Returns current value of host's fitness.
 * 
 * @return current value of host's fitness
 */
double Host::getFitness(){
    return Fitness;
}

/**
 * @brief Core method. Clears data regarding infections and fitness.
 */
void Host::clearInfections(){
//    PathogesInfecting.clear();
//    PathogesPresented.clear();
    NumOfPathogesInfecting = 0;
    NumOfPathogesPresented = 0;
    SelectedForReproduction = 0;
    Fitness = 0.0;
}

/**
 * @brief Data harvesting method. Gives a host's genome in a human-readable 
 * format. With all the gene specs.
 * 
 * @return a STL string containing the host's genome and its annotations in
 *  a human-readable format.
 */
 std::string Host::stringChromosomes(){
    std::string outString;
    outString = sttr(" === Host has ") +  std::to_string(NumOfPathogesInfecting) +
        sttr(" parasites and presented ") + std::to_string(NumOfPathogesPresented) +
        sttr(" of them ===\n");
    sttr g1;
    sttr g2;
    for(int i = 0; i < ChromosomeOne.size(); ++i){
        boost::to_string(ChromosomeOne[i].getBitGene(), g1);
        outString += sttr(g1) + sttr("\tch_one\t") 
                   + std::to_string(ChromosomeOne[i].timeOfOrigin) + sttr("\t")
                   + std::to_string(ChromosomeOne[i].GenesTag);
        if (ChromosomeOne[i].ParentTags.size()){
            for (int j = 0; j < ChromosomeOne[i].ParentTags.size(); ++j){
                outString += sttr("\t") + std::to_string(ChromosomeOne[i].MutationTime[j])
                        + sttr("\t") + std::to_string(ChromosomeOne[i].ParentTags[j]);
            }
            outString += sttr("\n");
        } else {
            outString += sttr("\t-1\n");
        }
    }
//    outString += sttr("-----\n");
    for(int k = 0; k < ChromosomeTwo.size(); ++k){
        boost::to_string(ChromosomeTwo[k].getBitGene(), g2);
        outString += sttr(g2) + sttr("\tch_two\t") 
                   + std::to_string(ChromosomeTwo[k].timeOfOrigin) + sttr("\t")
                   + std::to_string(ChromosomeTwo[k].GenesTag);
        if (ChromosomeTwo[k].ParentTags.size()){
            for (int l = 0; l < ChromosomeTwo[k].ParentTags.size(); ++l){
                outString += sttr("\t") + std::to_string(ChromosomeTwo[k].MutationTime[l])
                        + sttr("\t") + std::to_string(ChromosomeTwo[k].ParentTags[l]);
            }
            outString += sttr("\n");
        } else {
            outString += sttr("\t-1\n");
        }       
    } 
    return outString;
}
 