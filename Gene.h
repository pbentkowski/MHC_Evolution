/* 
 * File:   Gene.h
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
 */

#include <iostream>
#include <vector>
#include <set>

#include "boost/dynamic_bitset.hpp"
#include "Random.h"
#include "Tagging_system.h"

#ifndef GENE_H
#define	GENE_H

typedef boost::dynamic_bitset<> genestring;

/**
 * @brief Core class. Stores and handles a single gene object. Has methods to 
 * access and to mutate a gene. This gene is used by the hosts class object to
 * represent a single MHC.
 */
class Gene {
public:
    Gene();
//    Gene(const Gene& orig);
    // === Core stuff ===
    virtual ~Gene();
    void setNewGene(unsigned long length, int timeStamp, Random& randGen, Tagging_system& tag);
    void setNewGene(unsigned long  length, unsigned long low_lim,
                    unsigned long up_lim, int timeStamp, Random& randGen, Tagging_system& tag);
    void setNewFixedGene(unsigned long length, int timeStamp, unsigned long fixedGene,
                         unsigned long int fixedTag);
    void mutateGeneWhole(double mut_prob_whole, int timeStamp, Random& randGen, Tagging_system& tag);
    void mutateGeneWhole(double mut_prob_whole, unsigned long low_lim,
                         unsigned long up_lim, int timeStamp, Random& randGen, Tagging_system& ta);
    void mutateGeneBitByBit(double pm_mut_probabl, int timeStamp, Random& randGen, Tagging_system& tag);
    void mutateBitByBitWithRestric(double pm_mut_probabl, int timeStamp,
                                   std::set<unsigned long>& noMutts, Random& randGen, Tagging_system& tag);
    genestring getBitGene();
    unsigned long int getTheRealGene();
    // === Data harvesting ===
    int timeOfOrigin;
    int TheParentWas;
    std::vector<unsigned long int> ParentTags;
    std::vector<int> MutationTime;
    unsigned long int GenesTag;
    void printGeneToScreen(std::string tagLine);
private:
    unsigned long TheGene;
    unsigned long BitStringLength;
};

#endif	/* GENE_H */

