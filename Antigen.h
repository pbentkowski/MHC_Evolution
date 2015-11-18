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
private:

};

#endif /* ANTIGEN_H */

