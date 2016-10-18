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

#include "Antigen.h"

typedef boost::dynamic_bitset<> anigenstring;
typedef std::vector<Antigen> antigenvector;

/**
 * @brief Core class. Stores and handles a single pathogen object. Each pathogen
 * can have multiple instances of Antigen class objects and stores them in a 
 * vector called PathogenProts.
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
    void setNewPathogen(int num_of_loci, unsigned long antigen_size, unsigned long mhcSize,
                        int species, int timeStamp);
    void setNewPathogenNthSwap(int num_of_loci, anigenstring antigen, unsigned long mhcSize,
                               int species, int timeStamp, int Nth);
    antigenvector getAllAntigens();
    void chromoMutProcess(double mut_probabl, unsigned long mhcSize, int timeStamp);
    void chromoMutProcessWithRestric(double mut_probabl, unsigned long mhcSize, int timeStamp,
                                     std::set<unsigned long>& noMutts);
    void setNewSpeciesNumber(int new_spp_num);
    anigenstring getSingleAntigen(int indx);
    int getSpeciesTag();
    void clearInfections();
    // === Data harvesting methods ===
    std::string stringGenesFromGenome();
    // === Auxiliary methods ===
    void printGenesFromGenome();
private:
    std::vector<Antigen> PathogenProts;
    int Species;
};

#endif	/* PATHOGEN_H */

