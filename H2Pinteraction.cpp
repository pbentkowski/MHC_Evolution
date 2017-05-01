/* 
 * File:   H2Pinteraction.cpp
 * Author: Piotr Bentkowski : bentkowski.piotr@gmail.com
 * 
 * Created on 20 February 2015, 18:48
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

#include <vector>
#include <iostream>

#include "H2Pinteraction.h"

typedef std::vector<unsigned long int> longIntVec;
typedef std::vector<Gene> chromovector;
typedef std::vector<Antigen> antigenvector;

H2Pinteraction::H2Pinteraction() {
}

//H2Pinteraction::H2Pinteraction(const H2Pinteraction& orig) {
//}

H2Pinteraction::~H2Pinteraction() {
}

/**
 * @brief Core method. Checks if the antigen is presented by a given MHC.
 * 
 * Takes the the vector with epitopes (series of long unsigned integers) and
 * checks if any of them is the same as the MHC gene (also presented as long
 * unsigned integer). See function Antigen::calculateEpitopes(int mhcSize) that
 * generates epitopes from a bit-string shaped antigen.
 * 
 * @param hostgen - a MHC gene
 * @param antigen - vector of epitopes in a form of long unsigned integers
 * @return 'true' if gene is presented, 'false' if it's not
 */
bool H2Pinteraction::presentAntigen(unsigned long int hostgen, longIntVec antigen){
    if(antigen.size() > 0){
        for(unsigned long i = 0; i < antigen.size(); ++i){
            if(antigen[i] == hostgen){ return true; }            
        }
        return false;
    } else {
        std::cout << "Error in H2Pinteraction::presentAntigen(): epitope "\
                  << "vector is empty." << std::endl;
        return false;
    }
}

/**
 * @brief Core method. Checks if a host gets infected with a pathogen. 
 * Heterozygote has an advantage here over homozygote and each species is 
 * allowed to infect a host only ONES.
 * 
 * Iterates through the host genome and trough the pathogen antigens checking 
 * are there any antigens which are being presented. If they are, then strike one
 * for the host and pathogen gets rejected, if there are not, then the host gets
 * infected and a point for the pathogen. If a species is already found in the 
 * host then the procedure is abandoned.
 * 
 * @param host - a Host-class object
 * @param patho - a Pathogen-class object
 */
void H2Pinteraction::doesInfectedHeteroOnePerSpec(Host& host, Pathogen& patho){
    antigenvector tmppatho = patho.getAllAntigens();
    if (host.PathoSpecInfecting.size()){
        // Making sure a pathogen species infects only ones
        for (unsigned long w = 0; w < host.PathoSpecInfecting.size(); ++w){
            if(host.PathoSpecInfecting[w] == patho.getSpeciesTag()) return;
        }
    }
    for(int i = 0; i < host.getChromoOneSize(); ++i){
        for(int j = 0; j < tmppatho.size(); ++j){
            if(presentAntigen(host.getChromosomeOne()[i].getTheRealGene(), tmppatho[j].getEpitopes())){
                // the pathogen gets presented, the host evades infection:
                host.NumOfPathogesPresented = host.NumOfPathogesPresented + 1;
                host.PathogesPresented.push_back(patho.getSpeciesTag());
                return;
            }
        }
    }
    for(int i = 0; i < host.getChromoTwoSize(); ++i){
        for(int j = 0; j < tmppatho.size(); ++j){
            if(presentAntigen(host.getChromosomeTwo()[i].getTheRealGene(), tmppatho[j].getEpitopes())){
                // the pathogen gets presented, the host evades infection:
                host.NumOfPathogesPresented = host.NumOfPathogesPresented + 1;
                host.PathogesPresented.push_back(patho.getSpeciesTag());
                return;
            }
        }


    // The host gets infected:
//    host.PathogesInfecting.push_back(patho.getSpeciesTag());
//    patho.HostsInfected.push_back(host.getHostIndvTag());
    host.NumOfPathogesInfecting = host.NumOfPathogesInfecting + 1;
    patho.NumOfHostsInfected = patho.NumOfHostsInfected + 1;
    host.PathoSpecInfecting.push_back(patho.getSpeciesTag());
    return;
}
