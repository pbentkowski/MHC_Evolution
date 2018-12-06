/*
 * File:   Gene.cpp
 * Author: Piotr Bentkowski : bentkowski.piotr@gmail.com
 *
 * Created on 24 November 2018
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
 * 
 */

#include <iostream>
#include "boost/dynamic_bitset.hpp"

#include "Gene.h"

typedef boost::dynamic_bitset<> genestring;

Gene::Gene() = default;

//Gene::Gene(const Gene& orig) {
//}

Gene::~Gene() = default;

/**
 * @brief Core method. Sets a new gene filling it with a random bit-string of
 * a given length.
 *
 * The value of bit-string is randomly selected from a range of values spanning
 * from 0 to 2^length -1 .
 *
 * @param length - number of bits in a gene.
 * @param timeStamp - current time (current number of the model iteration)
 * @randGen - instance of the PRNG
 */
void Gene::setNewGene(unsigned long length, int timeStamp, Random& randGen, Tagging_system& tag) {
    timeOfOrigin = timeStamp;
    TheParentWas = -1;
    BitStringLength = length;
    GenesTag = tag.getTag();
    TheGene = randGen.getRandomFromUniform(0, (unsigned int)std::pow(2, BitStringLength)-1);
}

/**
 * @brief Core method. Sets a new gene filling it with a random bit-string of
 * a given length.
 *
 * The value of the gene is given here by the user.
 *
 *
 * @param length - number of bits in a gene.
 * @param timeStamp - current time (current number of the model iteration)
 * @param fixedGene - the gene value
 * @param fixedTag - tag value
 */
void Gene::setNewFixedGene(unsigned long length, int timeStamp, unsigned long fixedGene,
        unsigned long fixedTag){
    timeOfOrigin = timeStamp;
    TheParentWas = -1;
    BitStringLength = length;
    GenesTag = fixedTag;
    TheGene = fixedGene;
}

/**
 *  @brief Core method. Sets a new gene filling it with a random bit-string of
 * a given length selected from an interval within given values.
 *
 * The value of bit-string is randomly selected from a range of values defined
 * in the parameter list .
 *
 * @param length - number of bits in a gene
 * @param low_lim - lower limit of the range of permitted bit-string values
 * @param up_lim - upper limit of the range of permitted bit-string values
 * @param timeStamp - current time (current number of the model iteration).
 */
void Gene::setNewGene(unsigned long length, unsigned long low_lim, unsigned long up_lim,
        int timeStamp, Random& randGen, Tagging_system& tag) {
    timeOfOrigin = timeStamp;
    TheParentWas = -1;
    GenesTag = tag.getTag();
    BitStringLength = length;
    unsigned long possible_max = std::pow(2, length)-1;
    if(possible_max < up_lim){
        std::cout << "Error in Gene::setNewGene(): Demanded upper limit"\
                  " is out of possible boundaries for a bit string of this"\
                  " length." << std::endl;
        up_lim = possible_max;
    }
    // assigned a random string of ones and zeros of given length
    TheGene = randGen.getRandomFromUniform(low_lim, up_lim);
//    // below is a added line -remove
//    std::cout << "low: " << low_lim << ", high: " << up_lim \
//              << ", gene: "<< TheGene << std::endl;
}

/**
 * @brief Core method. Mutates a gene by overwriting a whole new bit-string.
 *
 * @param mut_prob_whole - probability of a whole-gene mutation.
 * @param timeStamp - current time (current number of the model iteration).
 */
void Gene::mutateGeneWhole(double mut_prob_whole, int timeStamp, Random& randGen, Tagging_system& tag) {
    if(randGen.getUni() < mut_prob_whole){
        TheParentWas = (int) TheGene;
        ParentTags.push_back(GenesTag);
        MutationTime.push_back(timeOfOrigin);
        GenesTag = tag.getTag();
        TheGene = randGen.getRandomFromUniform(0, (unsigned long) std::pow(2, BitStringLength)-1);
        timeOfOrigin = timeStamp;
    }
}

/**
 * @brief Core method. Mutates a gene by overwriting a whole new bit-string by
 * a number given in a defined range.
 *
 * @param mut_prob_whole - probability of a whole-gene mutation.
 * @param low_lim - lower limit of possible gene value. 
 * @param up_lim - upper limit of possible gene value.
 * @param timeStamp - current time (current number of the model iteration).
 */
void Gene::mutateGeneWhole(double mut_prob_whole, unsigned long low_lim,
                           unsigned long up_lim, int timeStamp) {
    RandomNumbs * p_RandomNumbs = RandomNumbs::getInstance();
    if(p_RandomNumbs->NextReal(0.0, 1.0) < mut_prob_whole){
        unsigned long possible_max = (unsigned long) std::pow(2, BitStringLength)-1;
        if(possible_max < up_lim){
            std::cout << "Error in Gene::mutateGeneWhole(): Demanded upper limit"\
                    " is out of possible boundaries for a bit string of this"\
                    " length." << std::endl;
            up_lim = possible_max;
        }
        TheParentWas = (int) TheGene;
        ParentTags.push_back(GenesTag);
        MutationTime.push_back(timeOfOrigin);
        Tagging_system* pTagging_system = Tagging_system::getInstance();
        GenesTag = pTagging_system->getTag();
        TheGene = p_RandomNumbs->NextLongInt(low_lim, up_lim);
        timeOfOrigin = timeStamp;
    }
}

/**
 * @brief Core method. Iterates through a gene sequence and (if selected so)
 * flips the value of a single bit to an opposite one.
 *
 * @param pm_mut_probabl - probability of mutating a single bit.
 * @param timeStamp - current time (current number of the model iteration).
 */
void Gene::mutateGeneBitByBit(double pm_mut_probabl, int timeStamp) {
    unsigned long currentGene = TheGene;
    boost::dynamic_bitset<> bitgene(BitStringLength, TheGene);
    RandomNumbs * p_RandomNumbs = RandomNumbs::getInstance();
    boost::dynamic_bitset<>::size_type bitgeneSize = bitgene.size();
    for(boost::dynamic_bitset<>::size_type i = 0; i < bitgeneSize; ++i) {
        if(p_RandomNumbs->NextReal(0.0, 1.0) < pm_mut_probabl) {
            bitgene[i].flip();
        }
    TheGene = bitgene.to_ulong();
    }
    if(currentGene != TheGene){
        ParentTags.push_back(GenesTag);
        MutationTime.push_back(timeOfOrigin);
        Tagging_system* pTagging_system = Tagging_system::getInstance();
        GenesTag = pTagging_system->getTag();
        timeOfOrigin = timeStamp;
    }
}

/**
 * @brief Core method. Iterates through a gene sequence and (if selected so)
 * flips the value of a single bit to an opposite one. But certain bits are not
 * permitted to mutate - these are specified in a global std::set.
 * 
 * @param mut_probabl - probability of mutating a single bit.
 * @param timeStamp - current time (current number of the model iteration).
 * @param noMutts - a std::set containing indices of residues of the bit-string 
 * that are not allowed to change.
 */
void Gene::mutateBitByBitWithRestric(double pm_mut_probabl, int timeStamp, 
        std::set<unsigned long >& noMutts){
    unsigned long currentGene = TheGene;
    bool exists;
    boost::dynamic_bitset<> bitgene(BitStringLength, TheGene);
    RandomNumbs * p_RandomNumbs = RandomNumbs::getInstance();
    boost::dynamic_bitset<>::size_type bitgeneSize = bitgene.size();
    for(boost::dynamic_bitset<>::size_type i = 0; i < bitgeneSize; ++i) { 
        exists = noMutts.find((int) i) != noMutts.end();
        if(!exists and p_RandomNumbs->NextReal(0.0, 1.0) < pm_mut_probabl) {
            bitgene[i].flip();
        }
    TheGene = bitgene.to_ulong();
    }
    if(currentGene != TheGene){
        ParentTags.push_back(GenesTag);
        MutationTime.push_back(timeOfOrigin);
        Tagging_system* pTagging_system = Tagging_system::getInstance();
        GenesTag = pTagging_system->getTag();
        timeOfOrigin = timeStamp;
    }
}

/**
 * @brief Core method. Returns the gene in a bit-string format, can be used to
 * pass the gene string to an another method.
 *
 * return bitgene - a gene in a bit-string format, a boost::dynamic_bitset object.
 */
genestring Gene::getBitGene(){
    boost::dynamic_bitset<> bitgene(BitStringLength, TheGene);
    return bitgene;
}

/**
 * @brief Core method. Returns the gene it it integer form. The gene is a private
 * data field to prevent it from being changed accidently.
 *
 * @return integer representation of a gene.
 */
unsigned long int Gene::getTheRealGene(){
    return TheGene;
}

/**
 * @brief Auxiliary method useful for debugging. Prints a gene to the screen.
 */
void Gene::printGeneToScreen(){
    boost::dynamic_bitset<> bitgene(BitStringLength, TheGene);
    std::cout << bitgene << std::endl;
}
