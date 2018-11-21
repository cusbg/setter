/*
 Copyright (c) 2013 David Hoksza

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef TMSCORE_HPP_
#define TMSCORE_HPP_

#include "align.hpp"
#include "vector.hpp"
#include "matrix.hpp"
#include "protein.hpp"
#include "export.hpp"



/**
 * TM-score mód
 */
typedef enum
{
    TM_FULL,
    TM_FAST_SSE,
    TM_FAST_ESSE,   // experimental
    TM_FAST_eSSE,   // experimental
    TM_FAST 
}
tmmode_t;



/**
 * iterative TM-score
 * 
 * \param dpiterations  count of DP ieratations                                         (2?)
 * \param dpwidth       half! size of the DP belt                                       (21)
 * \param mode          TM-score mode
 * \param normParam     scaling parameter (TM-score = 1 - nn-TM-score / normParam)      (query->length)
 * \param qprot         query protein
 * \param dprot         database protein
 * \param align         alignment
 * \param stream        vmd script output stream
 */
double itmscore(int dpiterations, int dpwidth, tmmode_t mode, double normParam, cRMSDStruct *qprot, cRMSDStruct *dprot, cRMSDAlign &align, std::ostream *stream = NULL);


#endif /* TMSCORE_HPP_ */
