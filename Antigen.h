/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Antigen.h
 * Author: piotr
 *
 * Created on November 18, 2015, 6:23 PM
 */

#ifndef ANTIGEN_H
#define ANTIGEN_H

#include <iostream>
#include <vector>
#include <set>

#include "boost/dynamic_bitset.hpp"
#include "RandomNumbs.h"

typedef boost::dynamic_bitset<> genestring;

class Antigen {
public:
    Antigen();
    virtual ~Antigen();
    void setNewAntigen(int lenght, int timeStamp);
    void setNewFixedAntigen(int lenght, int timeStamp, int fixedGene,
                            unsigned long int fixedTag);
    void mutateAntigenBitByBit(double pm_mut_probabl, int timeStamp);
    void mutateAntgBitByBitWithRes(double pm_mut_probabl, int timeStamp,
                                   std::set<int>& noMutts);
    
    // === Data harvesting ===
    int timeOfOrigin;
    int TheParentWas;
    std::vector<unsigned long int> ParentTags;
    std::vector<int> MutationTime;
    unsigned long int GenesTag;
    void printAntigenToScreen();
private:
    genestring Antigen;
    std::vector<int> Epitopes;
};

#endif /* ANTIGEN_H */

