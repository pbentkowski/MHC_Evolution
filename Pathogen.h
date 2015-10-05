/* 
 * File:   Pathogen.h
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
#ifndef PATHOGEN_H
#define	PATHOGEN_H

#include <cstdlib>
#include <vector>
#include <string>
#include "boost/dynamic_bitset.hpp"

#include "Gene.h"

typedef boost::dynamic_bitset<> genestring;
typedef std::vector<Gene> chromovector;

/**
 * @brief Core class. Stores and handles a single pathogen object. Each pathogen
 * has multiple instances of Gene class objects and stores them in a vector called
 * ChromosomePat.
 */
class Pathogen {
public:
    Pathogen();
//    Pathogen(const Pathogen& orig);
    virtual ~Pathogen();
    // === Core methods ===
//    std::vector<unsigned long> HostsInfected; // which host are infected
    unsigned NumOfHostsInfected;  // how many host are infected
    int SelectedToReproduct;
    void setNewPathogen(int num_of_loci, int gene_size, int species, int timeStamp);
    void setNewPathoFixedGene(int num_of_loci, int gene_size, int species, int timeStamp,
                              int fixedGene, unsigned long int fixedTag);
    void setNewPathogen(int num_of_loci, int gene_size, int species,
                        int low_lim, int up_lim, int timeStamp);
    chromovector getChomosome();
    void chromoMutProcess(double mut_probabl, int timeStamp);
    void chromoMutProcessWithRestric(double mut_probabl, int timeStamp,
                                     std::set<int>& noMutts);
    genestring getSingleGene(int indx);
    int getSpeciesTag();
    void clearInfections();
    // === Data harvesting methods ===
    std::string stringGenesFromGenome();
    // === Auxiliary methods ===
    void printGenesFromGenome();
private:
    std::vector<Gene> ChromosomePat;
    int Species;
    int LowerGeneValue;
    int UpperGeneValue;
};

#endif	/* PATHOGEN_H */

