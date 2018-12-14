
//
// Created by piotr on 04/12/18.
// compile: g++ Tagging_system.cpp Gene.cpp Random.cpp mainTest.cpp -fopenmp -std=c++14
//     run: ./a.out
//

#include <vector>
#include <iostream>
#include <thread>     // for reading the number of concurrent threads supported
#include "omp.h"

//#include "Gene.h"
//#include "Random.h"
#include "Tagging_system.h"
#include "Gene.h"



int main(int argc, char** argv) {
    std::vector<Gene> chrom;
    std::vector<unsigned int> tagLine;
    Tagging_system tag;
    unsigned int numberOfThreads = 0;
    unsigned int NN = 50000;
    Random* mRandGenArr;
    if(numberOfThreads == 0)
    {
        numberOfThreads = std::thread::hardware_concurrency();
        if(numberOfThreads == 0) // if the value is not well defined or not computable, set at least 1 thread
            numberOfThreads = 1;
    }
    std::cout << "We have " << numberOfThreads << " threads" << std::endl;
    omp_set_num_threads(numberOfThreads);
    mRandGenArr = new Random[numberOfThreads];

    for(unsigned int i = 0; i < numberOfThreads; ++i)
        mRandGenArr[i].reseed(10*i);

    Random* randGen_ptr = mRandGenArr;
    for(unsigned int k = 0; k < NN; ++k) {
        chrom.push_back(Gene());
        tagLine.push_back(999);
    }
//#pragma omp parallel default(none) shared(randGen_ptr, tag)
#pragma omp parallel shared(randGen_ptr, tag)
    for (int i = 0; i < NN; ++i) {
        chrom[i].setNewGene(16, 0, randGen_ptr[omp_get_thread_num()], tag);
        tagLine[i] = omp_get_thread_num();
    }
//    for (auto it = vect.begin(); it != vect.end(); ++it){
//        std::cout << *it << " ";
//    }
    for (unsigned int j = 0; j < NN; ++j) {
        chrom[j].printGeneToScreen(std::to_string(tagLine[j]));
    }
    std::cout << std::endl;
}
