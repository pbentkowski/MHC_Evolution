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
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
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
#include "Host.h"
#include "Pathogen.h"
#include "RandomNumbs.h"

typedef boost::dynamic_bitset<> genestring;
typedef std::vector<Gene>  chromovector;


H2Pinteraction::H2Pinteraction() {
}

//H2Pinteraction::H2Pinteraction(const H2Pinteraction& orig) {
//}

H2Pinteraction::~H2Pinteraction() {
}

/**
 * @brief Compares two genes if they have any N bits similar to each other in
 * coresponding positions.
 * 
 * @param hostgene - first gene
 * @param pathogene - second gene
 * @param simil_mesure - how many bits need to be similar (N)
 * @return 'true' if gene is presented, 'false' if it's not
 */
bool H2Pinteraction::presentGeneAny(genestring hostgene, genestring pathogene,
        int simil_mesure){
    if(hostgene.size() == pathogene.size()){
        int counter = 0;
        for(int i = 0; i < hostgene.size(); ++i){
            if(hostgene[i] == pathogene[i]) counter++;
        }
        if(counter >= simil_mesure){
            return true;
        } else {
            return false;
        }
    } else {
        std::cout << "Error in compareGenes(): genes are of different lengths."\
                  << " Can\'t compare genes of different lengths." << std::endl;
        return false;
    }
}


/**
 * @brief Compares two genes if they have a span of length of N bits similar 
 * to each other.
 * 
 * @param hostgene - first gene
 * @param pathogene - second gene
 * @param simil_mesure - how many bits need to be similar (N)
 * @return 'true' if gene is presented, 'false' if it's not
 */
bool H2Pinteraction::presentGeneRow(genestring hostgene, genestring pathogene,
        int simil_mesure){
    if(hostgene.size() == pathogene.size()){
        int counter = 0;
        for(int i = 0; i < hostgene.size(); ++i){
            if(hostgene[i] == pathogene[i]){
                counter = counter + 1;
            }else{
                counter = 0;
            }
            if(counter >= simil_mesure){return true;}
        }
        if(counter >= simil_mesure){
            return true;
        } else {
            return false;
        }
    } else {
        std::cout << "Error in H2Pinteraction::presentGeneRow(): genes are of"\
                  << " different lengths. Can\'t compare genes of different"\
                  << " lengths." << std::endl;
        return false;
    }
}


/**
 * @brief Core method. Checks if a host gets infected with a pathogen. 
 * Heterozygote has an advantage here over homozygote.
 * 
 * Iterates through the host genome and trough the pathogen genome checking 
 * are there any genes which are being presented. If they are, then strike one
 * for the host and pathogen gets rejected, if there are not, then the host gets
 * infected and a point for the pathogen.
 * 
 * @param host - a Host-class object
 * @param patho - a Pathogen-class object
 * @param simil_mesure - how many similar bits in a row two genes need to have
 */
void H2Pinteraction::doesInfectedHeteroBetter(Host &host, Pathogen &patho, 
        int simil_mesure){
    H2Pinteraction H2P;
    chromovector tmppatho = patho.getChomosome();
    chromovector tmphost = host.getChromosomeOne();
    for(int i = 0; i < tmphost.size(); ++i){
        for(int j = 0; j < tmppatho.size(); ++j){
            if(H2P.presentGeneRow(tmphost[i].getBitGene(),
                    tmppatho[j].getBitGene(), simil_mesure)){
                // the pathogen gets presented, the host evades infection:
//                host.PathogesPresented.push_back(patho.getSpeciesTag());
                host.NumOfPathogesPresented = host.NumOfPathogesPresented + 1;
                return;
            }
        }
    }
    tmphost = host.getChromosomeTwo();
    for(int i = 0; i < tmphost.size(); ++i){
        for(int j = 0; j < tmppatho.size(); ++j){
            if(presentGeneRow(tmphost[i].getBitGene(),
                    tmppatho[j].getBitGene(), simil_mesure)){
                // the pathogen gets presented, the host evades infection:
//                host.PathogesPresented.push_back(patho.getSpeciesTag());
                host.NumOfPathogesPresented = host.NumOfPathogesPresented + 1;
                return;
            }
        }
    }
    // The host gets infected:
//    host.PathogesInfecting.push_back(patho.getSpeciesTag());
//    patho.HostsInfected.push_back(host.getHostIndvTag());
    host.NumOfPathogesInfecting = host.NumOfPathogesInfecting + 1;
    patho.NumOfHostsInfected = patho.NumOfHostsInfected + 1;
    return;
}

/**
 * @brief Core method. Checks if a host gets infected with a pathogen. 
 * Heterozygote has NO advantage here over homozygote.
 * 
 * Iterates through the host genome and trough the pathogen genome checking 
 * are there any genes which are being presented. If they are, then good
 * for the host and pathogen gets rejected, if there are not, then the host gets
 * infected and a point for the pathogen. If the host is a homozygote, then
 * it scores double on the fitness scale then a homozygote.
 * 
 * @param host - a Host-class object
 * @param patho - a Pathogen-class object
 * @param simil_mesure - how many similar bits in a row two genes need to have
 */
void H2Pinteraction::doesInfectedHomoBetter(Host& host, Pathogen& patho,
        int simil_mesure){
    H2Pinteraction H2P;
    chromovector tmppatho = patho.getChomosome();
    chromovector tmphost = host.getChromosomeOne();
    bool flagIfInfects = true;  // pathogen does infect 
    for(int i = 0; i < tmphost.size(); ++i){
        for(int j = 0; j < tmppatho.size(); ++j){
            if(H2P.presentGeneRow(tmphost[i].getBitGene(),
                    tmppatho[j].getBitGene(), simil_mesure)){
                // the pathogen gets presented, the host evades infection:
//                host.PathogesPresented.push_back(patho.getSpeciesTag());
                host.NumOfPathogesPresented = host.NumOfPathogesPresented + 1;
                flagIfInfects = false; // pathogen gets presented 
                goto nextChromosome;
            }
        }
    }
    nextChromosome:
    tmphost = host.getChromosomeTwo();
    for(int i = 0; i < tmphost.size(); ++i){
        for(int j = 0; j < tmppatho.size(); ++j){
            if(presentGeneRow(tmphost[i].getBitGene(),
                    tmppatho[j].getBitGene(), simil_mesure)){
                // the pathogen gets presented, the host evades infection:
//                host.PathogesPresented.push_back(patho.getSpeciesTag());
                host.NumOfPathogesPresented = host.NumOfPathogesPresented + 1;
                flagIfInfects = false; // pathogen gets presented 
                return;
            }
        }
    }
    // The host gets infected:
//    host.PathogesInfecting.push_back(patho.getSpeciesTag());
//    patho.HostsInfected.push_back(host.getHostIndvTag());
    if(flagIfInfects){ 
        host.NumOfPathogesInfecting = host.NumOfPathogesInfecting + 1;
        patho.NumOfHostsInfected = patho.NumOfHostsInfected + 1;
        return;
    }
}

void H2Pinteraction::doesInfectedAllToAll(Host& host, Pathogen& patho,
        int simil_mesure){
     H2Pinteraction H2P;
    chromovector tmppatho = patho.getChomosome();
    chromovector tmphost = host.getChromosomeOne();
    bool flagIfInfects = true;  // pathogen does infect 
    for(int i = 0; i < tmphost.size(); ++i){
        for(int j = 0; j < tmppatho.size(); ++j){
            if(H2P.presentGeneRow(tmphost[i].getBitGene(),
                    tmppatho[j].getBitGene(), simil_mesure)){
                // the pathogen gets presented, the host evades infection:
//                host.PathogesPresented.push_back(patho.getSpeciesTag());
                host.NumOfPathogesPresented = host.NumOfPathogesPresented + 1;
                flagIfInfects = false; // pathogen gets presented 
            }
        }
    }
    tmphost = host.getChromosomeTwo();
    for(int i = 0; i < tmphost.size(); ++i){
        for(int j = 0; j < tmppatho.size(); ++j){
            if(presentGeneRow(tmphost[i].getBitGene(),
                    tmppatho[j].getBitGene(), simil_mesure)){
                // the pathogen gets presented, the host evades infection:
//                host.PathogesPresented.push_back(patho.getSpeciesTag());
                host.NumOfPathogesPresented = host.NumOfPathogesPresented + 1;
                flagIfInfects = false; // pathogen gets presented 
            }
        }
    }
    // The host gets infected:
//    host.PathogesInfecting.push_back(patho.getSpeciesTag());
//    patho.HostsInfected.push_back(host.getHostIndvTag());
    if(flagIfInfects){ 
        host.NumOfPathogesInfecting = host.NumOfPathogesInfecting + 1;
        patho.NumOfHostsInfected = patho.NumOfHostsInfected + 1;
        return;
    }
}