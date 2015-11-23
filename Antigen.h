/* 
 * File:   Antigen.h
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

#ifndef ANTIGEN_H
#define ANTIGEN_H

#include <iostream>
#include <vector>
#include <set>

#include "boost/dynamic_bitset.hpp"
#include "RandomNumbs.h"

typedef boost::dynamic_bitset<> antigentring;
typedef std::vector<unsigned long int> longIntVec;

class Antigen {
public:
    Antigen();
    virtual ~Antigen();
    void calculateEpitopes(int mhcSize);
    void setNewAntigen(int length, int mhcSize, int timeStamp);
    void setNewFixedAntigen(int length, int timeStamp, int fixedGene,
                            unsigned long int fixedTag);
    void mutateAntigenBitByBit(double pm_mut_probabl, int mhcSize, int timeStamp);
    void mutateAntgBitByBitWithRes(double pm_mut_probabl, int mhcSize, 
                                   int timeStamp, std::set<int>& noMutts);
    antigentring getBitAntigen();
    unsigned long int getOneEpitope(int idx);
    longIntVec getEpitopes();
    // === Data harvesting ===
    int timeOfOrigin;
    int TheParentWas;
    std::vector<unsigned long int> ParentTags;
    std::vector<int> MutationTime;
    unsigned long int AntigenTag;
    void printAntigenToScreen();
private:
    antigentring TheAntigen;
    std::vector<unsigned long int> Epitopes;
    int BitStringLength;
};

#endif /* ANTIGEN_H */

