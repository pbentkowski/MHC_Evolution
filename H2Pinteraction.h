/* 
 * File:   H2Pinteraction.h
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
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *    MA 02110-1301, USA.
 */
#ifndef H2PINTERACTION_H
#define	H2PINTERACTION_H

#include <vector>

#include "Host.h"
#include "Pathogen.h"
#include "RandomNumbs.h"

typedef std::vector<unsigned long int> longIntVec;

/**
 * @brief Core class. Handles interactions between hosts and pathogens.
 */
class H2Pinteraction {
public:
    H2Pinteraction();
//    H2Pinteraction(const H2Pinteraction& orig);
    virtual ~H2Pinteraction();
    bool presentAntigen(unsigned long int hostgen, longIntVec antigen);
    void doesInfectedHeteroOnePerSpecTotalGenome(Host &host, Pathogen &patho);
    void doesInfectedHeteroOnePerSpecUniqe(Host &host, Pathogen &patho);
};

#endif	/* H2PINTERACTION_H */

