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
#include "boost/dynamic_bitset.hpp"

#include "Pathogen.h"

typedef boost::dynamic_bitset<> antigenstring;
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
 * user-defined number of bits (as antigens are bit-strings). Value of a antigens
 * is randomly assigned from an interval between 0 and 2^N-1, where N is the
 * number of bits in a gene.
 *
 * @param num_of_loci - number of gene loci in a chromosome
 * @param gene_size - the length of the bit-string representing a gene
 * @param species - user-defined number of species
 * @param timeStamp - current time (current number of the model iteration)
 */
void Pathogen::setNewPathogen(unsigned long antigen_size, unsigned long mhcSize, int species, int timeStamp,
                              Random& randGen, Tagging_system& tag){
    Species = species;
    NumOfHostsInfected = 0;
    SelectedToReproduct = 0;
    PathoProtein.setNewAntigen(antigen_size, mhcSize, timeStamp, randGen, tag);
}

/**
 * @brief Core method. Sets a new Pathogen object with pre-defined antigen.
 * 
 * Sets a new Pathogen object and assigns a user-defined number of antigens of
 * user-defined number of bits (as antigens are bit-strings). Value of an 
 * antigen is loaded as a bit-string. Later on every N-th bit in the antigen
 * is flipped to the opposite value.
 * 
 * @param antigen - bit string containing pre-defined antigen
 * @param mhcSize - number of bits in MHC protein
 * @param species- user-defined number of species
 * @param timeStamp - current time (current number of the model iteration)
 * @param Nth - step at each a bit should be flipped
 */
void Pathogen::setNewPathogenNthSwap(anigenstring antigen, unsigned long int Tag, unsigned long mhcSize,
                                     int species, int timeStamp, int Nth){
    Species = species;
    NumOfHostsInfected = 0;
    SelectedToReproduct = 0;
    PathoProtein.setAntigenFlipedPositions(antigen, Tag, Nth, mhcSize, timeStamp);
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
void Pathogen::chromoMutProcess(double mut_probabl, unsigned long mhcSize, int timeStamp,
        Random& randGen, Tagging_system& tag){
    PathoProtein.mutateAntigenBitByBit(mut_probabl, mhcSize, timeStamp, randGen, tag);
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
void Pathogen::chromoMutProcessWithRestric(double mut_probabl, unsigned long mhcSize,
        int timeStamp, std::set<unsigned long>& noMutts, Random& randGen, Tagging_system& tag){
    PathoProtein.mutateAntgBitByBitWithRes(mut_probabl, mhcSize, timeStamp, noMutts, randGen, tag);
}



/**
 * @brief Changes pathogen's species tag. Don't over use it. 
 * 
 * @param new_spp_num
 */
void Pathogen::setNewSpeciesNumber(int new_spp_num){
   Species = new_spp_num;
}

/**
 * @brief Core method. Fetches the pathogene's antigens.
 *
 * @return - vector of genes object (a chromosome).
 */
Antigen Pathogen::getAntigenProt(){
    return PathoProtein;
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
 * @brief Data harvesting method. Returns a string containing the whole pathogen
 * genome in a human-readable format.
 *
 * @return STD string being in human readable format.
 */
std::string Pathogen::stringGenesFromGenome(){
    sttr genomeString;
    genomeString = sttr(" === Patho. sp. No. ") + std::to_string(Species) +
            sttr(" has infected ") + std::to_string(NumOfHostsInfected) +
            sttr(" hosts ===\n");
    sttr bitAntigen;
    boost::to_string(PathoProtein.getBitAntigen(), bitAntigen);
    genomeString += sttr(bitAntigen) + sttr("\tch_pat\t")
               + std::to_string(PathoProtein.timeOfOrigin) + sttr("\t")
               + std::to_string(PathoProtein.AntigenTag);
    unsigned long PathoProtsIthParentTagsSize = PathoProtein.ParentTags.size();
    for (unsigned long j = 0; j < PathoProtsIthParentTagsSize; ++j){
        genomeString += sttr("\t") + std::to_string(PathoProtein.ParentTags[j]);
    }
    genomeString += sttr("\n");

    return genomeString;
}

/**
 * @brief Auxiliary method useful for debugging. Prints genome to screen.
 */
void Pathogen::printGenesFromGenome(){
    std::cout << " === Patho. sp. No. " <<  Species << " has infected " \
              << NumOfHostsInfected << " hosts ===" << std::endl;
    std::cout << PathoProtein.getBitAntigen() << std::endl;
}
