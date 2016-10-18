/* 
 * File:   Gene.h
 * Author: Piotr Bentkowski : bentkowski.piotr@gmail.com
 *
 * Created on 13 February 2015, 13:26
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
#ifndef GENE_H
#define	GENE_H

#include <iostream>
#include <vector>
#include <set>

#include "boost/dynamic_bitset.hpp"
#include "RandomNumbs.h"

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
    void setNewGene(unsigned long lenght, int timeStamp);
    void setNewGene(unsigned long  lenght, unsigned long low_lim,
                    unsigned long up_lim, int timeStamp);
    void setNewFixedGene(unsigned long lenght, int timeStamp, unsigned long fixedGene,
                         unsigned long int fixedTag);
    void mutateGeneWhole(double mut_prob_whole, int timeStamp);
    void mutateGeneWhole(double mut_prob_whole, unsigned long low_lim,
                         unsigned long up_lim, int timeStamp);
    void mutateGeneBitByBit(double pm_mut_probabl, int timeStamp);
    void mutateBitByBitWithRestric(double pm_mut_probabl, int timeStamp,
                                   std::set<unsigned long>& noMutts);
    genestring getBitGene();
    unsigned long int getTheRealGene();
    // === Data harvesting ===
    int timeOfOrigin;
    int TheParentWas;
    std::vector<unsigned long int> ParentTags;
    std::vector<int> MutationTime;
    unsigned long int GenesTag;
    void printGeneToScreen();
private:
    unsigned long TheGene;
    unsigned long BitStringLength;
};

#endif	/* GENE_H */

