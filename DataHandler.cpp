/* 
 * File:   DataHarvester.cpp
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

#include <fstream>
#include <cstdlib>
#include <iostream>
#include <string>
#include <stdio.h>
#include <time.h>
#include <vector>
#include <tuple>

#include "DataHandler.h"
#include"Gene.h"

typedef std::string sttr;

DataHandler::DataHandler() {
}

//DataHarvester::DataHarvester(const DataHarvester& orig) {
//}

DataHandler::~DataHandler() {
}

/**
 * @brief Data harvesting method. Calculates some stats of population genetics:
 * Shannon's index, number of MHC/antigen types, total number of MHC copies. 
 * Runs internally within DataHarvester class.
 * 
 * @return <a href="http://en.cppreference.com/w/cpp/utility/tuple/make_tuple">
 * STD Touple object</a> with 3 numbers: Shannon's index, number of MHC/antigen 
 * types, total number of MHC copies
 */
auto getShannonIndx(std::vector<int> GeneVals){
    int Types = 0;
    double Summ = 0.0;
    double tot_gene_numb = (double) GeneVals.size();
    if(GeneVals.size()){
        int GeneCounter = GeneVals.size();
        bool IfCountedLyst[GeneCounter];
        for (int w = 0; w < GeneCounter; ++w) {
            IfCountedLyst[w] = false;
        }
        double abundance;
        for (int i = 0; i < GeneCounter; ++i) {
            if (IfCountedLyst[i] == false) {
                Types += 1;
                abundance = 1.0;
                for (int j = i + 1; j < GeneCounter; ++j) {
                    if (i != j && IfCountedLyst[j] == false 
                               && GeneVals[i] == GeneVals[j]){
                        abundance += 1.0;
                        IfCountedLyst[j] = true;
                    }
                }
                // the proper calculations
                Summ = Summ - ((abundance / tot_gene_numb) 
                        * std::log(abundance / tot_gene_numb));
            }
        }
    }
    // ShannIndx (double), typesOfGenes (int), numberOfAllGenes (double)
    auto RTRN = std::make_tuple( Summ, Types, tot_gene_numb );
    return RTRN;
}

/**
 * @brief Data harvesting method. Sets status of all data files as "brand new".
 */
void DataHandler::setAllFilesAsFirtsTimers(){
    ifFirstSpecToFileRun = true;
    ifFirstHostClonesRun = true;
    ifFirstHostGeneDivRun = true;
    ifFirstGeneNumbersTotal = true;
    ifFirstGeneNumbersUnique = true;
    ifNoMuttPathoListUnique = true;
}

/** 
 * @brief Data harvesting method. Gets current date/time, format is YYYY-MM-DD.HH:mm:ss
 * 
 * Visit http://en.cppreference.com/w/cpp/chrono/c/strftime for more information 
 * about date/time format.
 * 
 * @return date and time in a string format ready to print out.
 */
const std::string currentDateTime(){
    time_t now = time(0);
    struct tm tstruct;
    char buf[80];
    tstruct = *localtime(&now);
    strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);
    return buf;
}

/**
 * @brief Input params validation method. Does the basic check if the entered 
 * parameters are free of total nonsense.
 * 
 * @param rndSeed
 * @param geneLength
 * @param antigenLength
 * @param hostPopSize
 * @param pathoPopSize
 * @param patho_sp
 * @param hostGeneNumbb
 * @param pathoGeneNumb
 * @param patoPerHostGeneration
 * @param numOfHostGenerations
 * @param hostMutationProb
 * @param pathoMutationProb
 * @param HeteroHomo
 * @param hostDeletion
 * @param hostDuplication
 * @param maxGene
 * @param alpha
 * @param fixedAntigPosit
 * @return 'true' if something is wrong, 'false' if no errors were found.
 */
bool DataHandler::checkParamsIfWrong(int rndSeed, int geneLength, int antigenLength,
        int hostPopSize, int pathoPopSize, int patho_sp, int hostGeneNumbb,
        int pathoGeneNumb, int patoPerHostGeneration, int numOfHostGenerations,
        double hostMutationProb, double pathoMutationProb, int HeteroHomo,
        double hostDeletion, double hostDuplication, int maxGene, double alpha,
        double fixedAntigPosit){
    bool ifError = false;
    if (rndSeed < 0){
        std::cout << "\nError in RNG seed. It has to be a positive integer!." << std::endl;
        ifError = true;
    }
    if (geneLength > 31){
        std::cout << "\nError in number of bits per gene. ";
        std::cout << geneLength << " bits in genes is bit too much. ";
        std::cout << "Try something less radical, e.g. smaller than 31."  << std::endl;
        ifError =  true;
    } 
    if (geneLength > antigenLength){
        std::cout << "\nError in size of an antigen. ";
        std::cout << " Number of bits in antigen cannot be smaller";
        std::cout << " than the number of bits in MHC gene." << std::endl;
        ifError =  true;
    }
    
    if (hostMutationProb < 0.0 or hostMutationProb > 1.0){
        std::cout << "\nError in the hosts' mutation probability. It has to be " <<
                "within the range [0, 1]." << std::endl;
        ifError = true;
    }
    if (pathoMutationProb < 0.0 or pathoMutationProb > 1.0){
        std::cout << "\nError in the pathogens' mutation probability. It has " <<
                "to be within the range [0, 1]." << std::endl;
        ifError = true;
    }
    if (hostDeletion < 0.0 or hostDeletion > 1.0){
        std::cout << "\nError in the hosts' probability of deletion of a gene. " <<
                "It has to be within the range [0, 1]." << std::endl;
        ifError = true;
    }
    if (hostDuplication < 0.0 or hostDuplication > 1.0){
        std::cout << "\nError in the hosts' duplication of a gene probability. " <<
                "It has to be within the range [0, 1]." << std::endl;
        ifError = true;
    }
    if (HeteroHomo != 10 and HeteroHomo != 11){
        std::cout << "\nError in the heterozygote advantage / lack of advantage " <<
                "mode. It has to be 10 for heterozygote advantage or 11 for " <<
                "lack of thereof." << std::endl;
        ifError = true;
    }
    if (maxGene < 1.0 or maxGene < hostGeneNumbb){
        std::cout << "\nError in the hosts' maximal number of genes per " <<
                "chromosome. It has to be at least one, but not less then the " <<
                "number used to initialize the system. "<< std::endl;
        ifError = true;
    }
    if (alpha < 0.0 or alpha > 1.0){
        std::cout << "\nError in the hosts' alpha factor for the host fitness function. " <<
                "It has to be within the range [0, 1]." << std::endl;
        ifError = true;
    }
    if (fixedAntigPosit < 0.0 or fixedAntigPosit > 1.0){
        std::cout << "\nError in the parameter for fraction of antigen bits " <<
                "being fixed. It has to be within the range [0, 1]." << std::endl;
        ifError = true;
    }
    return ifError;
}

/**
 * @brief Data harvesting method. Writes all the input params and some run stats
 * to a file. Must be run only ones per run of the model.
 * 
 * It's better to run it after running DataHarvester::checkParamsIfWrong() which
 * will check if parameters make any sense.
 * 
 * @param rndSeed
 * @param geneLength
 * @param antigenLength
 * @param hostPopSize
 * @param pathoPopSize
 * @param patho_sp
 * @param hostGeneNumbb
 * @param pathoGeneNumb
 * @param patoPerHostGeneration
 * @param numOfHostGenerations
 * @param hostMutationProb
 * @param pathoMutationProb
 * @param HeteroHomo
 * @param hostDeletion
 * @param hostDuplication
 * @param maxGene
 * @param alpha
 * @param fixedAntigPosit
 */
void DataHandler::inputParamsToFile(int rndSeed, int geneLength, int antigenLength,
        int hostPopSize, int pathoPopSize, int patho_sp, int hostGeneNumbb,
        int pathoGeneNumb, int patoPerHostGeneration, int numOfHostGenerations,
        double hostMutationProb, double pathoMutationProb, int HeteroHomo,
        double hostDeletion, double hostDuplication, int maxGene, double alpha,
        double fixedAntigPosit){
    std::ofstream InputParams;
    InputParams.open("InputParameters.csv");

    InputParams << "# Runtime properties:" << std::endl;
    InputParams << "\trun_start_date_and_time = " << currentDateTime() << std::endl;
    InputParams << "# Model's core parameters:" << std::endl;
    InputParams << "\trandom_number_seed = " << rndSeed << std::endl;
    InputParams << "\tnumber_of_bits_per_gene = " << geneLength << std::endl;
    InputParams << "\tnumber_of_bits_per_antigen = " << antigenLength << std::endl;
    InputParams << "\thost_population_size = " << hostPopSize << std::endl;
    InputParams << "\tpathogen_population_size = " << pathoPopSize << std::endl;
    InputParams << "\tnumber_of_pathogen_species = " << patho_sp << std::endl;
    InputParams << "\tnumber_of_genes_per_host_one_chromosome = " << 
        hostGeneNumbb << std::endl;
    InputParams << "\tnumber_of_antigens_per_pathogen = " << pathoGeneNumb << std::endl;
    InputParams << "\tnumber_of_pathogen_generation_per_one_host_generation = " <<
        patoPerHostGeneration << std::endl;
    InputParams << "\tnumber_of_host_generations = " << numOfHostGenerations << std::endl;
    InputParams << "\tmutation_probability_in_host = " << 
            hostMutationProb << std::endl;
    InputParams << "\tmutation_probability_in_pathogen = " << 
            pathoMutationProb << std::endl;
    if (HeteroHomo == 10){
        InputParams << "\theterozygote_advantage = YES" << std::endl;
    }else if (HeteroHomo == 11){
        InputParams << "\theterozygote_advantage = NO" << std::endl;
    } else {
        InputParams << "\theterozygote_advantage = ERROR" << std::endl;
    }
    InputParams << "\thost_gene_deletion_probability = " <<
            hostDeletion << std::endl;
    InputParams << "\thost_gene_duplication_probability = " <<
            hostDuplication << std::endl;
    InputParams << "\thost_maximal_number_of_genes_in_chromosome = " <<
            maxGene << std::endl;
    InputParams << "\tAlpha_factor_for_the_host_fitness_function = " <<
            alpha << std::endl;
    InputParams << "\tFraction_of_antigen_bits_getting_fixed = " <<
            fixedAntigPosit << std::endl;
    InputParams.close();
}
    
/**
 * @brief Data harvesting method. Writes to a file population sizes of all 
 * pathogen species in a given time.
 * 
 * @param EnvObj - the Environment object
 * @param tayme - time stamp 
 */
void DataHandler::saveNumOfPathoSpeciesToFile(Environment& EnvObj, int tayme){
    if(ifFirstSpecToFileRun){
        std::ofstream PathoPopulFile;
        PathoPopulFile.open("PathoPopSizes.csv");
        PathoPopulFile << "#time subsequent_species_pop_size" << std::endl;
        PathoPopulFile.close();
        ifFirstSpecToFileRun = false;
    }
    std::ofstream PathoPopulFile;
    PathoPopulFile.open("PathoPopSizes.csv",
            std::ios::out | std::ios::ate | std::ios::app);
    PathoPopulFile << tayme;
    for(unsigned j = 0; j < EnvObj.getPathoNumOfSpecies(); ++j){
        PathoPopulFile << " " << EnvObj.getPathoSpeciesPopSize(j);      
    }
    PathoPopulFile << std::endl;
    PathoPopulFile.close();
}

/**
 * @brief Data harvesting method. Writes to a file all pathogens with their
 * genomes in a human-readable format.
 * 
 * @param EnvObj - the Environment object
 * @param tayme - time stamp
 */
void DataHandler::savePathoPopulToFile(Environment& EnvObj, int tayme){
    sttr theFilename = sttr("PathoGenomesFile.") + std::to_string(tayme) + sttr(".csv");
    std::ofstream PathogGenomeFile;
    PathogGenomeFile.open(theFilename);
    PathogGenomeFile << "#Genomes_of_all_pathogens_at_time = " << tayme << std::endl;
    PathogGenomeFile << "#bit-gene\tchromosome\ttime_of_origin\tgene_own_tag\tAll_parental_tags"
            << std::endl;
    for(int i = 0; i < EnvObj.getPathoNumOfSpecies(); ++i){
        for(int j = 0; j < EnvObj.getPathoSpeciesPopSize(i); ++j){
            PathogGenomeFile << EnvObj.getPathoGenesToString(i, j);
        }
    }
    PathogGenomeFile.close();
}

/**
 * @brief Data harvesting method. Writes to a file all hosts with their
 * genomes in a human-readable format.
 * 
 * @param EnvObj - the Environment object
 * @param tayme - time stamp
 */
void DataHandler::saveHostPopulToFile(Environment& EnvObj, int tayme){
    sttr theFilename = sttr("HostGenomesFile.") + std::to_string(tayme) + sttr(".csv");
    std::ofstream HostGenomesFile;
    HostGenomesFile.open(theFilename);
    HostGenomesFile << "#Genomes_of_all_host_at_time = " << tayme << std::endl;
    HostGenomesFile << "#bit-gene\tchromosome\ttime_of_origin\tgene_own_tag" <<
            "\tTime_of_parental_mutation\tAll_parental_tags etc."
            << std::endl;
    for(int i = 0; i < EnvObj.getHostsPopSize(); ++i){
        HostGenomesFile << EnvObj.getHostGenesToString(i);
    }
    HostGenomesFile.close();
}

/**
 * @brief Data harvesting method. Calculates and writes to a file some stats 
 * about the hosts population genetic. Will create one file per run.
 * 
 * You call it ones per host generation iteration. In each call it adds a line
 * with some statistics regarding population genetics of the host population.
 * Data regarding individuals are not being stored here. We wrote custom Python
 * scripts to visualize this dataset.
 * 
 * @param EnvObj - the Environment object
 * @param tayme - time stamp
 */
void DataHandler::saveHostGeneticDivers(Environment& EnvObj, int tayme){
    if(ifFirstHostGeneDivRun){
        std::ofstream HostGeneDivFile;
        HostGeneDivFile.open("HostsGeneDivers.csv");
        HostGeneDivFile << "#time pop_size tot_num_of_genes num_of_MHC_types" <<
                " Shannon_indx mean_fitness std_fitness" << std::endl;
        HostGeneDivFile.close();
        ifFirstHostGeneDivRun = false;
    }
    std::vector<int> AllTheGeneVals;
    std::vector<double> Fitness;
    double popSize = (double) EnvObj.getHostsPopSize();
    int homoLociNum = -1; // not applicable at the moment!
//    double tot_gene_numb = 0;
    Fitness.clear();
    for(int i = 0; i < EnvObj.getHostsPopSize(); ++i){
//        tot_gene_numb += EnvObj.getSingleHostGenomeSize(i);
        Fitness.push_back(EnvObj.getHostFitness(i));
        // Harvesting genes from Chromosome One in all hosts
        for(int j = 0; j < EnvObj.getSingleHostChromoOneSize(i); ++j){
            AllTheGeneVals.push_back(EnvObj.getSingleHostRealGeneOne(i, j));
        }
        // Harvesting genes from Chromosome Two in all hosts
        for(int m = 0; m < EnvObj.getSingleHostChromoTwoSize(i); ++m){
            AllTheGeneVals.push_back(EnvObj.getSingleHostRealGeneTwo(i, m));
        }
    }
    // Calculating the Shannon index et al. plus extracting it from the tuple
    auto ShOut = getShannonIndx(AllTheGeneVals);
    double Summ = std::get < 0 >( ShOut );
    int mhcTypes = std::get < 1 >( ShOut );
    double tot_gene_numb = std::get < 2 >( ShOut );
    
    // calculating coefficient of variation of hosts fitness
    double cv_sum = std::accumulate(Fitness.begin(), Fitness.end(), 0.0);
    double cv_mean = cv_sum / Fitness.size();
    double sq_sum = std::inner_product(Fitness.begin(), Fitness.end(), 
                                       Fitness.begin(), 0.0);
    double stdev = std::sqrt(sq_sum / Fitness.size() - cv_mean * cv_mean);
    
    std::ofstream HostGenomesFile;
    HostGenomesFile.open("HostsGeneDivers.csv",
            std::ios::out | std::ios::ate | std::ios::app);
    HostGenomesFile << tayme << " " << popSize << " " << tot_gene_numb << " " <<
            mhcTypes << " " << Summ << " " << cv_mean <<
            " " << stdev << std::endl;
    HostGenomesFile.close();
}

/**
 * @brief  Data harvesting method. Record the total number of genes and types
 * of MHCs in hosts. All the freaking hosts!
 * 
 * Writes to separate files (a) the number of  genes of all the host in each 
 * time step, (b) the number of unique MHC types corresponding to the numbers
 *  in the first file. All together you get two files witch data sync to each
 * other.
 * 
 * @param EnvObj - the Environment class object
 * @param tayme - time stamp
 */
void DataHandler::saveHostGeneNumbers(Environment& EnvObj, int tayme){
    if(ifFirstGeneNumbersTotal){
        std::ofstream HostGeneNumbTotal;
        HostGeneNumbTotal.open("HostGeneNumbTotal_ChrOne.csv");
        HostGeneNumbTotal << "#time total_number_of_genes_in_all_host_cells"
                << std::endl;
        HostGeneNumbTotal.close();
        ifFirstGeneNumbersTotal = false;
    }
    if(ifFirstGeneNumbersUnique){
        std::ofstream HostMHCNumbUnique;
        HostMHCNumbUnique.open("HostMHCsNumbUniq_ChrOne.csv");
        HostMHCNumbUnique << "#time number_of_unique_MHCs_in_all_host_cells" 
                << std::endl;
        HostMHCNumbUnique.close();
        ifFirstGeneNumbersUnique = false;
    }
    std::vector<unsigned int> AllGenomesSize;
    std::vector<int> UniqueMHCs;
    std::vector<int> TheGeneVals;
    for(int i = 0; i < EnvObj.getHostsPopSize(); ++i){
        TheGeneVals.clear();
        // Harvesting genes from Chromosome One in all hosts
        for(int j = 0; j < EnvObj.getSingleHostChromoOneSize(i); ++j){
            TheGeneVals.push_back(EnvObj.getSingleHostRealGeneOne(i, j));
        }
        // Harvesting genes from Chromosome Two in all hosts
//        for(int m = 0; m < EnvObj.getSingleHostChromoTwoSize(i); ++m){
//            TheGeneVals.push_back(EnvObj.getSingleHostRealGeneTwo(i, m));
//        }
        AllGenomesSize.push_back(TheGeneVals.size());
        auto ShOut = getShannonIndx(TheGeneVals);
        int mhcTypes = std::get < 1 >( ShOut );
        UniqueMHCs.push_back(mhcTypes);
    }
    std::ofstream HostGeneNumbTotal;
    HostGeneNumbTotal.open("HostGeneNumbTotal_ChrOne.csv",
                            std::ios::out | std::ios::ate | std::ios::app);
    HostGeneNumbTotal << tayme;
    for(unsigned int ii = 0; ii < AllGenomesSize.size(); ++ii){
        HostGeneNumbTotal << " " << AllGenomesSize[ii];
    }
    HostGeneNumbTotal << std::endl;
    HostGeneNumbTotal.close();
    
    std::ofstream HostMHCNumbUnique;
    HostMHCNumbUnique.open("HostMHCsNumbUniq_ChrOne.csv",
                            std::ios::out | std::ios::ate | std::ios::app);
    HostMHCNumbUnique << tayme;
    for(unsigned int jj = 0; jj < UniqueMHCs.size(); ++jj){
        HostMHCNumbUnique << " " << UniqueMHCs[jj];
    }
    HostMHCNumbUnique << std::endl;
    HostMHCNumbUnique.close();
}

/**
 * @brief Data harvesting method. Record the indices of fixed bits in all
 * pathogen species. Run this just ones. 
 * 
 * @param EnvObj - the Environment class object
 */
void DataHandler::savePathoNoMuttList(Environment& EnvObj){
    if(ifNoMuttPathoListUnique){
        std::ofstream NoMuttPathoList;
        NoMuttPathoList.open("NoMutationInPathoList.csv");
        NoMuttPathoList << "#list_of_bits_excluded_from_mutating_Each_line_is_a_spp"
                << std::endl;
        NoMuttPathoList.close();
        ifNoMuttPathoListUnique = false;
    }
    std::ofstream NoMuttPathoList;
    NoMuttPathoList.open("NoMutationInPathoList.csv",
                            std::ios::out | std::ios::ate | std::ios::app);
    NoMuttPathoList << EnvObj.getFixedBitsInAntigens();
    NoMuttPathoList.close();
}
