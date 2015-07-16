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
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
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

typedef boost::dynamic_bitset<> genestring;

/**
 * @brief Core class. Handles interactions between the hosts and the pathogens.
 */
class H2Pinteraction {
public:
    H2Pinteraction();
//    H2Pinteraction(const H2Pinteraction& orig);
    virtual ~H2Pinteraction();
    bool presentGeneAny(genestring hostgene, genestring pathogene, int simil_mesure);
    bool presentGeneRow(genestring hostgene, genestring pathogene, int simil_mesure);
    void doesInfectedHeteroBetter(Host &host, Pathogen &patho, int simil_mesure);
    void doesInfectedHomoBetter(Host &host, Pathogen &patho, int simil_mesure);
    void doesInfectedAllToAll(Host &host, Pathogen &patho, int simil_mesure);
};

#endif	/* H2PINTERACTION_H */

