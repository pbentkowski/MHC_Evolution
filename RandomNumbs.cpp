/* 
 * File:   RandomNumbs.cpp
 * Author: Piotr Bentkowski : bentkowski.piotr@gmail.com
 * 
 * Created on 12 February 2015, 17:34
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

#include <iostream>
#include <random>
#include "RandomNumbs.h"


RandomNumbs *RandomNumbs::s_instance = 0;

/**
 * @brief Core method. Seeds the PRNG. Remember to call as the first function
 * after initializing this class.
 */ 
void RandomNumbs::SetSeed(int seed){
    rg.seed( seed );
}

/**
 * @brief Core method. Gets the instance of a singleton class.
 */
 RandomNumbs* RandomNumbs::getInstance() {
    if (!s_instance)
        s_instance = new RandomNumbs;
    return s_instance;
}
 
/**
 * @brief Core method. Generates one integer random number from a given interval
 * with an uniform distribution.
 */
int RandomNumbs::NextInt(int lowerLim, int upperLim) {
    std::uniform_int_distribution<> distrib(lowerLim, upperLim);
    return distrib(rg);
}

/**
 * @brief Core method. Generates one real random number from a given interval
 * with an uniform distribution.
 */
double RandomNumbs::NextReal(double lowerLim, double upperLim) {
    std::uniform_real_distribution<> distrib(lowerLim, upperLim);
    return distrib(rg);
}


