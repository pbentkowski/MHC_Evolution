/* 
 * File:   RandomNumbs.h
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
#ifndef RANDOMNUMBS_H
#define	RANDOMNUMBS_H

#include <iostream>
#include <random>

/**
 * @brief Core method. Provides the Marsenne Twister random number generator and
 * methods of accessing it. 
 * 
 * Generates random numbers. It's a singleton class which runs the Marsenne Twister
 * random number generator implemented in STD C++11 libraries. This way we are
 * sure that all the random numbers come from the very same run of the PRNG
 * and are not correlated with each other. 
 * More about <a href="http://en.cppreference.com/w/cpp/numeric/random">
 * std::random library can be found in this on-line handbook</a>. 
 * More about why randomness is important can be found on 
 * <a href="https://www.random.org/randomness/">www.random.org</a>.
 * 
 * Usage of this PRNG engine:
 * 1. Initialize:    RandomNumbs * p_RandomNumbs = RandomNumbs::getInstance();
 * 2. Seed the PRNG: p_RandomNumbs->SetSeed(seed);
 * 3. get the values: 
 *   a) get an integer: int x = p_RandomNumbs->NextInt(1, 10);
 *   b) get a real num: double x = p_RandomNumbs->NextReal(0.0, 1.5);
 * 4. Your program needs to be compile in C++11 standard, so in g++ compiler use
 *  the flag: -std=c++11
 */
class RandomNumbs {
private:
    static RandomNumbs *s_instance;
public:
    typedef std::mt19937 RandomGeneratorType;
    void SetSeed(int seed);
    static RandomNumbs *getInstance();
    // the proper functions which throw random numbers:
    int NextInt(int lowerLim, int upperLim);
    double NextReal(double lowerLim, double upperLim);
private:
    RandomGeneratorType rg;
};

#endif	/* RANDOMNUMBS_H */

