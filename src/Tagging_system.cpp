/*      tagging_system.cpp
 *
 *      Copyright 2015 Piotr Bentkowski <bentkowski.piotr@gmail.com>
 *
 *      This program is free software; you can redistribute it and/or modify
 *      it under the terms of the GNU General Public License as published by
 *      the Free Software Foundation; either version 2 of the License, or
 *      (at your option) any later version.
 *
 *      This program is distributed in the hope that it will be useful,
 *      but WITHOUT ANY WARRANTY; without even the implied warranty of
 *      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *      GNU General Public License for more details.
 *
 *      You should have received a copy of the GNU General Public License
 *      along with this program; if not, write to the Free Software
 *      Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *      MA 02110-1301, USA.
 */

#include "Tagging_system.h"

unsigned long int theTag = 0;
/**
 * @brief Data collecting method. Constructor.
 */
Tagging_system::Tagging_system()
{
    theTag = 0;
    omp_init_lock(&lockTagCreation);
}

unsigned long int Tagging_system::getTag()
{
    omp_set_lock(&lockTagCreation);
    theTag++;
    omp_unset_lock(&lockTagCreation);
    return theTag;
}
