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
#include "Tagging_system.h"
#include "RandomNumbs.h"

typedef boost::dynamic_bitset<> antigentring;
typedef std::vector<unsigned long int> longIntVec;

Antigen::Antigen() {
}

Antigen::~Antigen() {
}


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
void Antigen::calculateEpitopes(int mhcSize){
    int vecSize =  TheAntigen.size() - mhcSize;
    boost::dynamic_bitset<> bitEpitope(mhcSize);
    longIntVec tmpEpis(vecSize);
    for(int i = 0; i < tmpEpis.size(); ++i){
        for(int j = i; j < (i + mhcSize); ++j){
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
void Antigen::setNewAntigen(int length, int mhcSize, int timeStamp){
    timeOfOrigin = timeStamp;
    TheParentWas = -1;
    BitStringLength = length;
    Tagging_system* pTagging_system = Tagging_system::getInstance();
    AntigenTag = pTagging_system->getTag();
    TheAntigen.clear();
    RandomNumbs * p_RandomNumbs = RandomNumbs::getInstance();
    boost::dynamic_bitset<> tmpAntig(length);
    for(boost::dynamic_bitset<>::size_type i = 0; i < tmpAntig.size(); ++i){
        if(p_RandomNumbs->NextReal(0.0, 1.0) < 0.5){
            tmpAntig[i] = true;
        } else {
            tmpAntig[i] = false;
        }
    }
    TheAntigen = tmpAntig;
    calculateEpitopes(mhcSize);
}


/**
 * @brief Core method. Mutates antigen one bit by one bit.
 * 
 * @param pm_mut_probabl - probability of mutating a single bit.
 * @param mhcSize - length of a bit string representing the MHC protein.
 * @param timeStamp - current time (current number of the model iteration).
 */
void Antigen::mutateAntigenBitByBit(double pm_mut_probabl, int mhcSize, int timeStamp){
    boost::dynamic_bitset<> bitgene;
    bitgene = TheAntigen;
    RandomNumbs * p_RandomNumbs = RandomNumbs::getInstance();
    for(boost::dynamic_bitset<>::size_type i = 0; i < bitgene.size(); ++i) {
        if(p_RandomNumbs->NextReal(0.0, 1.0) < pm_mut_probabl) {
            bitgene[i].flip();
        }
    }
    if(TheAntigen != bitgene){
        ParentTags.push_back(AntigenTag);
        MutationTime.push_back(timeOfOrigin);
        Tagging_system* pTagging_system = Tagging_system::getInstance();
        AntigenTag = pTagging_system->getTag();
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
void Antigen::mutateAntgBitByBitWithRes(double pm_mut_probabl, int mhcSize, 
        int timeStamp, std::set<int>& noMutts){
    boost::dynamic_bitset<> bitgene;
    bool exists;
    bitgene = TheAntigen;
    RandomNumbs * p_RandomNumbs = RandomNumbs::getInstance();
    for(boost::dynamic_bitset<>::size_type i = 0; i < bitgene.size(); ++i) {
        exists = noMutts.find(i) != noMutts.end();
        if(exists == false and p_RandomNumbs->NextReal(0.0, 1.0) < pm_mut_probabl) {
            bitgene[i].flip();
        }
    }
    if(TheAntigen != bitgene){
        ParentTags.push_back(AntigenTag);
        MutationTime.push_back(timeOfOrigin);
        Tagging_system* pTagging_system = Tagging_system::getInstance();
        AntigenTag = pTagging_system->getTag();
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
antigentring Antigen::getBitAntigen(){
    return TheAntigen;
}


unsigned long int Antigen::getOneEpitope(int idx){
    if(Epitopes.size()){
        return Epitopes[idx];
    } else {
        return -1;
    }
}

longIntVec Antigen::getEpitopes(){
    return Epitopes;
}


/**
 * @brief Auxiliary method. Prints antigens to screen. Useful when debugging.
 */
void Antigen::printAntigenToScreen(){
    std::cout << TheAntigen << std::endl;
    for(int i = 0; i < Epitopes.size(); ++i){
      std::cout <<  Epitopes[i] << " ";
    }
    std::cout << std::endl;
}