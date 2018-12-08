
//
// Created by piotr on 04/12/18.
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
    Tagging_system tag;
    unsigned int numberOfThreads = 0;
    Random* mRandGenArr;
    if(numberOfThreads == 0)
    {
        numberOfThreads = std::thread::hardware_concurrency() - 1;
        if(numberOfThreads == 0) //if the value is not well defined or not computable, set at least 1 thread
            numberOfThreads = 1;
    }
    std::cout << "We have " << numberOfThreads << "threads" << std::endl;
    omp_set_num_threads(numberOfThreads);
    mRandGenArr = new Random[numberOfThreads];

    for(unsigned int i = 0; i < numberOfThreads; ++i)
        mRandGenArr[i].reseed(10*i);

    for (int k = 0; k < 1000; ++k) {
        chrom.push_back(Gene());
    }
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
