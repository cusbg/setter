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
#ifndef EXPORT_HPP_
#define EXPORT_HPP_

#include <iostream>
#include "vector.hpp"
#include "matrix.hpp"
#include "protein.hpp"



#define QUERY_PATH              "c:/svnsiret/experiments/RNASimilarity/" 
#define QUERY_ALIGN_COLOR       "0" // blue
#define QUERY_ALIGN_STYLE       "Trace"
#define QUERY_NOTALIGN_COLOR    "0" // blue
#define QUERY_NOTALIGN_STYLE    "Trace"

#define DATA_PATH               "c:/svnsiret/experiments/RNASimilarity/" 
#define DATA_ALIGN_COLOR        "1" // red
#define DATA_ALIGN_STYLE        "Trace"
#define DATA_NOTALIGN_COLOR     "1" // red
#define DATA_NOTALIGN_STYLE     "Trace"

#define ALIGN_LINE_COLOR        "5" // tan
#define BALL_RADIUS             "0.1"
#define LINE_WIDTH              "2"

/*
#define MAXALIGN_LINE_COLOR     "6" // silver
#define MAXBALL_RADIUS          "0.1"
#define MAXLINE_WIDTH           "2"
*/

#define MAXALIGN_LINE_COLOR     "16" // black
#define MAXBALL_RADIUS          "0"
#define MAXLINE_WIDTH           "1"

#define CUTALIGN_LINE_COLOR     "7"
#define CUTBALL_RADIUS          "0.1"
#define CUTLINE_WIDTH           "2"



void generate_vmd(std::ostream &stream, cRMSDStruct *qprot, cRMSDStruct *dprot, cRMSDAlign &align, cRMSDAlign &align_max, cRMSDAlign &cut_max, cRMSDMatrix &u, cRMSD3DCoord &xt, cRMSD3DCoord &yt);


#endif /* VMD_HPP_ */
