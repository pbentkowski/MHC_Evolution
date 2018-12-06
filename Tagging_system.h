/*      tagging_system.h
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
#ifndef _TAGGING_SYSTEM_
#define _TAGGING_SYSTEM_

#include <omp.h>
#include <atomic>

extern unsigned long int theTag;

/**
 * @class Tagging_system
 * 
 * @brief A singleton class which generates unique tags for genes. These tags 
 * are used when keeping a track of evolution of individual alleles.
 * 
 */
class Tagging_system {
private:
    //    unsigned long int theTag;
    omp_lock_t lockTagCreation;

public:
    Tagging_system();
    unsigned long int getTag();
};

#endif  /* TAGGING_SYSTEM */
