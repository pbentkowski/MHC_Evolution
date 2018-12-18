/* 
 * File:   DataHarvester.h
 * Author: Piotr Bentkowski : bentkowski.piotr@gmail.com
 *
 * Created on 13 March 2015, 13:16
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

#ifndef DATAHARVESTER_H
#define	DATAHARVESTER_H

#include <cstdlib>
#include <vector>
#include <boost/lexical_cast.hpp>
//#include "boost/dynamic_bitset.hpp"

#include "Environment.h"

/**
 * @brief Data harvesting class. A class that has methods to collect data and
 * writes them to files.
 */
class DataHandler {
public:
    DataHandler();
//    DataHarvester(const DataHarvester& orig);
    virtual ~DataHandler();
    bool checkParamsIfWrong(unsigned int numberOfThreads, unsigned long geneLength, unsigned long antigenLength,
        int hostPopSize, int pathoPopSize, int patho_sp, unsigned long hostGeneNumbb,
        int pathoGeneNumb, int patoPerHostGeneration, int numOfHostGenerations,
        double hostMutationProb, double pathoMutationProb, int HeteroHomo,
        double hostDeletion, double hostDuplication, unsigned long maxGene, double alpha,
        double fixedAntigPosit);
    bool checkParamsIfWrong(unsigned int numberOfThreads, unsigned long geneLength, unsigned long antigenLength,
        int hostPopSize, int pathoPopSize, int patho_sp, unsigned long hostGeneNumbb,
        int pathoGeneNumb, int patoPerHostGeneration, int numOfHostGenerations,
        double hostMutationProb, double pathoMutationProb, int HeteroHomo,
        double hostDeletion, double hostDuplication, unsigned long maxGene, int numberOfMates,
        double fixedAntigPosit);
    bool checkParamsIfWrong(unsigned int numberOfThreads, unsigned long geneLength, int hostPopSize,
        int hostGeneNumbb, int numOfHostGenerations, double hostMutationProb,
        int HeteroHomo, double hostDeletion, double hostDuplication, int maxGene,
        unsigned long numberOfMates);
    bool checkParamsIfWrong(unsigned int numberOfThreads, unsigned long geneLength, unsigned long antigenLength,
        int hostPopSize, int pathoPopSize, int patho_sp, unsigned long hostGeneNumbb,
        int pathoGeneNumb, int patoPerHostGeneration, int numOfHostGenerations,
        double hostMutationProb, double pathoMutationProb, int HeteroHomo,
        double hostDeletion, double hostDuplication, unsigned long maxGene, double alpha,
        int numberOfMates); // for alpha and optimal mating scenario
    void inputParamsToFile(unsigned int numberOfThreads, unsigned long geneLength, unsigned long antigenLength,
                           int hostPopSize, int pathoPopSize, int patho_sp, unsigned long hostGeneNumbb,
                           int pathoGeneNumb, int patoPerHostGeneration, int numOfHostGenerations,
                           double hostMutationProb, double pathoMutationProb, int HeteroHomo,
                           double hostDeletion, double hostDuplication, unsigned long maxGene, double alpha,
                           double fixedAntigPosit);
    void inputParamsToFile(unsigned int numberOfThreads, unsigned long geneLength, unsigned long antigenLength,
                           int hostPopSize, int pathoPopSize, int patho_sp, unsigned long hostGeneNumbb,
                           int pathoGeneNumb, int patoPerHostGeneration, int numOfHostGenerations,
                           double hostMutationProb, double pathoMutationProb, int HeteroHomo,
                           double hostDeletion, double hostDuplication, unsigned long maxGene, int numberOfMates,
                           double fixedAntigPosit);
    void inputParamsToFile(unsigned int numberOfThreads, unsigned long geneLength, int hostPopSize,
                           int hostGeneNumbb, int numOfHostGenerations, double hostMutationProb,
                           int HeteroHomo, double hostDeletion, double hostDuplication, unsigned long maxGene,
                           unsigned long numberOfMates);
    void inputParamsToFile(unsigned int numberOfThreads, unsigned long geneLength, unsigned long antigenLength,
                           int hostPopSize, int pathoPopSize, int patho_sp, unsigned long hostGeneNumbb,
                           int pathoGeneNumb, int patoPerHostGeneration, int numOfHostGenerations,
                           double hostMutationProb, double pathoMutationProb, int HeteroHomo,
                           double hostDeletion, double hostDuplication, unsigned long maxGene, double alpha,
                           int numberOfMates); // for alpha and optimal mating scenario
    void setAllFilesAsFirtsTimers();
    void saveNumOfPathoSpeciesToFile(Environment &EnvObj, int tayme);
    void savePathoPopulToFile(Environment &EnvObj, int tayme);
    void saveHostPopulToFile(Environment &EnvObj, int tayme);
    void saveHostGeneticDivers(Environment &EnvObj, int tayme);
    void saveHostGeneNumbers(Environment &EnvObj, int tayme);
    void savePathoNoMuttList(Environment &EnvObj);
    void savePresentedPathos(Environment &EnvObj, int tayme);
    void saveMhcNumbersWhenMating(Environment &EnvObj, int tayme);
private:
    bool ifFirstSpecToFileRun;
    bool ifFirstHostClonesRun;
    bool ifFirstHostGeneDivRun;
    bool ifFirstGeneNumbersTotal;
    bool ifFirstGeneNumbersUnique;
    bool ifNoMuttPathoListUnique;
    bool ifNumberOfPresentedPatho;
    bool ifNumberOfMhcWhenMating;
};

#endif	/* DATAHARVESTER_H */

