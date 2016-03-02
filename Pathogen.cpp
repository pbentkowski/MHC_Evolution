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
#include "Tagging_system.h"

typedef boost::dynamic_bitset<> antigentring;
typedef std::vector<Antigen> antigenvector;
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
void Pathogen::setNewPathogen(int num_of_loci, int antigen_size, int mhcSize, 
        int species, int timeStamp){
    Species = species;
    NumOfHostsInfected = 0;
    SelectedToReproduct = 0;
    for(int i = 0; i < num_of_loci; ++i){
        PathogenProts.push_back(Antigen());
        PathogenProts.back().setNewAntigen(antigen_size, mhcSize, timeStamp);
    }
}

void Pathogen::setNewPathogenNthSwap(int num_of_loci, anigenstring antigenn, 
        int mhcSize, int species, int timeStamp, int Nth){
    Species = species;
    NumOfHostsInfected = 0;
    SelectedToReproduct = 0;
    Tagging_system* pTagging_system = Tagging_system::getInstance();
    unsigned long int Tag = pTagging_system->getTag();
    for(int i = 0; i < num_of_loci; ++i){
        PathogenProts.push_back(Antigen());
        PathogenProts.back().setAntigenFlipedPositions(antigenn, Tag, Nth, mhcSize, timeStamp);
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
void Pathogen::chromoMutProcess(double mut_probabl, int mhcSize, int timeStamp){
    for(int i = 0; i < PathogenProts.size(); ++i){
//        ChromosomePat[i].mutateGeneWhole(mut_probabl);
//        ChromosomePat[i].mutateGeneWhole(mut_probabl, LowerGeneValue,
//                UpperGeneValue, timeStamp);
        PathogenProts[i].mutateAntigenBitByBit(mut_probabl, mhcSize, timeStamp);
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
void Pathogen::chromoMutProcessWithRestric(double mut_probabl, int mhcSize, 
        int timeStamp, std::set<int>& noMutts){
    for(int i = 0; i < PathogenProts.size(); ++i){
        PathogenProts[i].mutateAntgBitByBitWithRes(mut_probabl, mhcSize, 
                                                   timeStamp, noMutts);
    }
}

/**
 * @brief Core method. Fetches a single antigen from a genome in a bit-string format.
 *
 * @param indx - index number of gene in a chromosome
 * @return - a Boost's dynamic bit string object
 */
antigentring Pathogen::getSingleAntigen(int indx){
    if(indx < PathogenProts.size()){
        return PathogenProts[indx].getBitAntigen();
    }else{
        std::cout << "Error in Host::getSingleGene(): Index out of the range of"\
                " the chromosome size. Fetching the last gene." << std::endl;
        return PathogenProts.back().getBitAntigen();
    }
}

/**
 * @brief Core method. Fetches the pathogene's antigens.
 *
 * @return - vector of genes object (a chromosome).
 */
antigenvector Pathogen::getAllAntigens(){
    return PathogenProts;
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
    for(int i = 0; i < PathogenProts.size(); i++){
        boost::to_string(PathogenProts[i].getBitAntigen(), bitAntigen);
        genomeString += sttr(bitAntigen) + sttr("\tch_pat\t")
                   + std::to_string(PathogenProts[i].timeOfOrigin) + sttr("\t")
                   + std::to_string(PathogenProts[i].AntigenTag);
        for (int j = 0; j < PathogenProts[i].ParentTags.size(); ++j){
            genomeString += sttr("\t") + std::to_string(PathogenProts[i].ParentTags[j]);
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
    for(int i = 0; i < PathogenProts.size(); i++){
        std::cout << PathogenProts[i].getBitAntigen() << std::endl;
    }
}
