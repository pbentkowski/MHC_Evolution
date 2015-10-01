/*
 * File:   Pathogen.cpp
 * Author: Piotr Bentkowski : bentkowski.piotr@gmail.com
 *
 * Created on 18 February 2015, 13:11
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

#include <complex>
#include <vector>
#include <string>
#include "boost/dynamic_bitset.hpp"

#include "Pathogen.h"

typedef boost::dynamic_bitset<> genestring;
typedef std::vector<Gene> chromovector;
typedef std::string sttr;

Pathogen::Pathogen() {
}

//Pathogen::Pathogen(const Pathogen& orig) {
//}

Pathogen::~Pathogen() {
}

/**
 * @brief Core method. Sets a new Pathogen object.
 *
 * Sets a new Pathogen object and assigns a user-defined number of genes of
 * user-defined number of bits (as genes are bit-strings). Value of a gene
 * is randomly assigned from an interval between 0 and 2^N-1, where N is the
 * number of bits in a gene.
 *
 * @param num_of_loci - number of gene loci in a chromosome
 * @param gene_size - the length of the bit-string representing a gene
 * @param species - user-defined number of species
 * @param timeStamp - current time (current number of the model iteration)
 */
void Pathogen::setNewPathogen(int num_of_loci, int gene_size, int species, int timeStamp){
    Species = species;
    NumOfHostsInfected = 0;
    SelectedToReproduct = 0;
    LowerGeneValue = 0;
    UpperGeneValue = std::pow(2, gene_size) - 1;
    for(int i = 0; i < num_of_loci; ++i){
        ChromosomePat.push_back(Gene());
        ChromosomePat.back().setNewGene(gene_size, timeStamp);
    }
}

/**
 * @brief Core method. Sets a new Pathogen object.
 *
 * Sets a new Pathogen object and assigns a user-defined number of genes
 * of user-defined number of bits (as genes are bit-strings). Value of a gene
 * is randomly assigned from a user-specified interval.
 *
 * @param num_of_loci - number of gene loci in a chromosome
 * @param gene_size - the length of the bit-string representing a gene
 * @param species - user-defined number of species
 * @param low_lim - lower limit of the range of permitted bit-string values
 * @param up_lim - upper limit of the range of permitted bit-string values
 * @param timeStamp - current time (current number of the model iteration)
 */
void Pathogen::setNewPathogen(int num_of_loci, int gene_size, int species,
        int low_lim, int up_lim, int timeStamp){
    Species = species;
    NumOfHostsInfected = 0;
    SelectedToReproduct = 0;
    LowerGeneValue = low_lim;
    UpperGeneValue = up_lim;
    for(int i = 0; i < num_of_loci; ++i){
        ChromosomePat.push_back(Gene());
        ChromosomePat.back().setNewGene(gene_size, low_lim, up_lim, timeStamp);
//        // remove line below
//        ChromosomePat.back().printGeneToScreen();
    }
}

/**
 * @brief Core method. Decides (on a random basis) if there will be any mutations
 * in the genome.
 *
 * Iterates through the chromosome and calls the gene mutation function. Genes
 * mutate at random depending on the probability which was user-defined.
 *
 * @param mut_probabl - mutation probability, a probability a gene will be
 * replaced by a new one
 * @param timeStamp - current time (current number of the model iteration)
 */
void Pathogen::chromoMutProcess(double mut_probabl, int timeStamp){
    for(int i = 0; i < ChromosomePat.size(); ++i){
//        ChromosomePat[i].mutateGeneWhole(mut_probabl);
//        ChromosomePat[i].mutateGeneWhole(mut_probabl, LowerGeneValue,
//                UpperGeneValue, timeStamp);
        ChromosomePat[i].mutateGeneBitByBit(mut_probabl, timeStamp);
    }
}

/**
 * @brief Core method. Decides (on a random basis) if there will be any mutations
 * in the genome. But some positions in the bit-string are not allowed to change. 
 * 
 * @param mut_probabl - mutation probability, a probability a gene will be
 * replaced by a new one
 * @param timeStamp - current time (current number of the model iteration)
 * @param noMutts - a std::set containing indices of residues of the bit-string 
 * that are not allowed to change
 */
void Pathogen::chromoMutProcessWithRestric(double mut_probabl, int timeStamp,
        std::set<int>& noMutts){
    for(int i = 0; i < ChromosomePat.size(); ++i){
        ChromosomePat[i].mutateBitByBitWithRestric(mut_probabl, timeStamp, noMutts);
    }
}

/**
 * @brief Core method. Fetches a single gene from a genome in a bit-string format.
 *
 * @param indx - index number of gene in a chromosome
 * @return - a Boost's dynamic bit string object
 */
genestring Pathogen::getSingleGene(int indx){
    if(indx < ChromosomePat.size()){
        return ChromosomePat[indx].getBitGene();
    }else{
        std::cout << "Error in Host::getSingleGene(): Index out of the range of"\
                " the chromosome size. Fetching the last gene." << std::endl;
        return ChromosomePat.back().getBitGene();
    }
}

/**
 * @brief Core method. Fetches the pathogene's chromosome.
 *
 * @return - vector of genes object (a chromosome).
 */
chromovector Pathogen::getChomosome(){
    return ChromosomePat;
}


/**
 * @brief Core method. Fetches a species tag.
 *
 * @return - integer being the species tag.
 */
int Pathogen::getSpeciesTag(){
    return Species;
}

/**
 * @brief Core method. Zeroing infection and reproduction flags.
 */
void Pathogen::clearInfections(){
//    HostsInfected.clear();
    NumOfHostsInfected = 0;
    SelectedToReproduct = 0;
}

/**
 * @brief Data harvesting method. Returns a sting containing the whole pathogen
 * genome in a human-readable format.
 *
 * @return STD string being in human readable format.
 */
std::string Pathogen::stringGenesFromGenome(){
    sttr genomeString;
    genomeString = sttr(" === Patho. sp. No. ") + std::to_string(Species) +
            sttr(" has infected ") + std::to_string(NumOfHostsInfected) +
            sttr(" hosts ===\n");
    sttr bitGene;
    for(int i = 0; i < ChromosomePat.size(); i++){
        boost::to_string(ChromosomePat[i].getBitGene(), bitGene);
        genomeString += sttr(bitGene) + sttr("\tch_pat\t")
                   + std::to_string(ChromosomePat[i].timeOfOrigin) + sttr("\t")
                   + std::to_string(ChromosomePat[i].GenesTag);
        for (int j = 0; j < ChromosomePat[i].ParentTags.size(); ++j){
            genomeString += sttr("\t") + std::to_string(ChromosomePat[i].ParentTags[j]);
        }
        genomeString += sttr("\n");
    }
    return genomeString;
}

/**
 * @brief Auxiliary method useful for debugging. Prints genome to screen.
 */
void Pathogen::printGenesFromGenome(){
    std::cout << " === Patho. sp. No. " <<  Species << " has infected " \
              << NumOfHostsInfected << " hosts ===" << std::endl;
    for(int i = 0; i < ChromosomePat.size(); i++){
        std::cout << ChromosomePat[i].getBitGene() << std::endl;
    }
}
