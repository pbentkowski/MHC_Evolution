/* 
 * File:   Antigen.cpp
 * Author: Piotr Bentkowski : bentkowski.piotr@gmail.com
 * 
 * Created on November 18, 2015, 18:23
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
#include "boost/dynamic_bitset.hpp"

#include "Antigen.h"

typedef boost::dynamic_bitset<> antigenstring;
typedef std::vector<unsigned long int> longIntVec;

Antigen::Antigen() = default;

Antigen::~Antigen() = default;


/**
 * @brief Core method. Translates a raw bit string representing the antigen into
 * series of epitopes.
 * 
 * Takes a frame of length equal to length of a bit string representing MHC 
 * protein and moves it along the antigen transforming bit strings inside the
 * frame into unsigned long integers useful for fast looking up if a pathogen
 * gets presented by its host, as MHC are also bit strings transformed into 
 * u_long int. Comparing the content of u_long int vectors is faster then 
 * matching bit strings.
 * 
 * @param mhcSize - length of a bit string representing the MHC protein.
 */
void Antigen::calculateEpitopes(unsigned long mhcSize){
    unsigned long vecSize =  TheAntigen.size() - mhcSize;
    boost::dynamic_bitset<> bitEpitope(mhcSize);
    longIntVec tmpEpis(vecSize);
    unsigned long tmpEpisSize = tmpEpis.size();
    for(unsigned long i = 0; i < tmpEpisSize; ++i){
        for(unsigned long j = i; j < (i + mhcSize); ++j){
            bitEpitope[j-i] = TheAntigen[j];
        }
        tmpEpis[i] = bitEpitope.to_ulong();
    }
    Epitopes = tmpEpis;
}


/**
 * @brief Core method. Sets a new antigen filling it with a random bits. 
 * 
 * Generates a bit string of a given length and fills it with random sequence
 * of bits. 
 * 
 * @param lenght - number of bits in the antigen (usually a lot)
 * @param mhcSize - length of a bit string representing the MHC protein.
 * @param timeStamp - current time (current number of the model iteration)
 */
void Antigen::setNewAntigen(unsigned long length, unsigned long mhcSize, int timeStamp, Random& randGen, Tagging_system& tag){
    timeOfOrigin = timeStamp;
    TheParentWas = -1;
    BitStringLength = length;
    AntigenTag = tag.getTag();
    TheAntigen.clear();
    boost::dynamic_bitset<> tmpAntig(length);
    boost::dynamic_bitset<>::size_type tmpAntigSize = tmpAntig.size();
    for(boost::dynamic_bitset<>::size_type i = 0; i < tmpAntigSize; ++i){
        tmpAntig[i] = randGen.getUni() < 0.5 ? true : false;
    }
    TheAntigen = tmpAntig;
    calculateEpitopes(mhcSize);
}

/**
 * @brief Core method. Sets a new antigen filling it with a given bit strings. 
 * 
 * Sets a antigen object filling it with a user-defined bit string.
 * 
 * @param bitgene - a bit string representing an antigen
 * @param Tag - antigen's individual ID tag
 * @param Nth - step at each a bit should be flipped
 * @param mhcSize - length of a bit string representing the MHC protein.
 * @param timeStamp - current time (current number of the model iteration)
 */
void Antigen::setAntigenFlipedPositions(antigenstring bitgene, unsigned long int Tag,
        int Nth, unsigned long mhcSize, int timeStamp){
    TheAntigen = bitgene;
    antigenstring::size_type TheAntigenSize = TheAntigen.size();
    for(antigenstring::size_type i = 0; i < TheAntigenSize; i += Nth){
        TheAntigen[i].flip();
    }
    timeOfOrigin = timeStamp;
    TheParentWas = -1;
    BitStringLength = bitgene.size();
    AntigenTag = Tag;
    calculateEpitopes(mhcSize);
}


/**
 * @brief Core method. Mutates antigen one bit by one bit.
 * 
 * @param pm_mut_probabl - probability of mutating a single bit.
 * @param mhcSize - length of a bit string representing the MHC protein.
 * @param timeStamp - current time (current number of the model iteration).
 */
void Antigen::mutateAntigenBitByBit(double pm_mut_probabl, unsigned long mhcSize, int timeStamp,
                                    Random& randGen, Tagging_system& tag){
    boost::dynamic_bitset<> bitgene;
    bitgene = TheAntigen;
    boost::dynamic_bitset<>::size_type bitgeneSize = bitgene.size();
    for(boost::dynamic_bitset<>::size_type i = 0; i < bitgeneSize; ++i) {
        if(randGen.getUni() < pm_mut_probabl) {
            bitgene[i].flip();
        }
    }
    if(TheAntigen != bitgene){
        ParentTags.push_back(AntigenTag);
        MutationTime.push_back(timeOfOrigin);
        AntigenTag = tag.getTag();
        timeOfOrigin = timeStamp;
        TheAntigen = bitgene;
        calculateEpitopes(mhcSize);
    }
}


/**
 * @brief Core method. Mutates antigen one bit by one bit but leaves predefined
 * positions on the antigen intact to make pathogen species a bit different.
 * 
 * @param pm_mut_probabl - probability of mutating a single bit.
 * @param mhcSize - length of a bit string representing the MHC protein.
 * @param timeStamp - current time (current number of the model iteration).
 * @param noMutts - a STL set of positions that should remain intact during
 * the mutation process, a way to define a species.  
 */
void Antigen::mutateAntgBitByBitWithRes(double pm_mut_probabl, unsigned long mhcSize,
        int timeStamp, std::set<unsigned long>& noMutts, Random& randGen, Tagging_system& tag){
    boost::dynamic_bitset<> bitgene;
    bitgene = TheAntigen;
    unsigned long bitgeneSize = bitgene.size();
    for(unsigned long i = 0; i < bitgeneSize; ++i) {
        if(noMutts.count(i) == 0){
            if(randGen.getUni() < pm_mut_probabl) {
                bitgene[i].flip(); 
           }
        }
    }
    if(TheAntigen != bitgene){
        ParentTags.push_back(AntigenTag);
        MutationTime.push_back(timeOfOrigin);
        AntigenTag = tag.getTag();
        timeOfOrigin = timeStamp;
        TheAntigen = bitgene;
        calculateEpitopes(mhcSize);
    }
}


/**
 * @brief Core method. Returns the antigen so other methods can use it.
 * 
 * @return The antigen in it's native bit format
 */
antigenstring Antigen::getBitAntigen(){
    return TheAntigen;
}

/**
 * @brief Core method. Returns the epitope with the given index.
 * 
 * @param idx - index of the interesitng epitope
 * @return the epitope (as a number)
 */
unsigned long int Antigen::getOneEpitope(unsigned long idx){
    if(idx < Epitopes.size() and idx >= 0){
        return Epitopes[idx];
    } else {
        return 0;
    }
}

/**
 * @brief Core method. Returns the epitopes that can be generated from one
 * antigen. 
 * 
 * It is a vector containing epitopes created by a sliding frame of a size equal
 * to the size of MHC bit string. Epitopes are represented by long unsigned 
 * integers.
 * 
 * @return a vector of long unsigned int numbers.
 */
longIntVec Antigen::getEpitopes(){
    return Epitopes;
}


/**
 * @brief Auxiliary method. Prints antigens to screen. Useful when debugging.
 */
void Antigen::printAntigenToScreen(){
    std::cout << TheAntigen << std::endl;
    unsigned long EpitopesSize = Epitopes.size();
    for(unsigned long i = 0; i < EpitopesSize; ++i){
      std::cout <<  Epitopes[i] << " ";
    }
    std::cout << std::endl;
}