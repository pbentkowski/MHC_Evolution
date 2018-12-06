
//
// Created by piotr on 04/12/18.
//

#include <vector>
#include <iostream>
#include "omp.h"

//#include "Gene.h"
//#include "Random.h"
#include "Tagging_system.h"



int main(int argc, char** argv) {
    unsigned long int vect[1000];
//    std::vector<unsigned long int> vect[1000];
    Tagging_system tag;
    #pragma omp parallel for
    for (int i = 0; i < 1000; ++i) {
        vect[i] = tag.getTag();
    }
//    for (auto it = vect.begin(); it != vect.end(); ++it){
//        std::cout << *it << " ";
//    }
    for (int j = 0; j < 1000; ++j) {
        std::cout << vect[j] << " ";
    }
    std::cout << std::endl;
}
