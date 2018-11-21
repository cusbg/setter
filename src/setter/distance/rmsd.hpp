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
#ifndef RMSD_HPP_
#define RMSD_HPP_

#include "align.hpp"
#include "vector.hpp"
#include "matrix.hpp"
#include "export.hpp"

double rmsd(cRMSDStruct *qprot, cRMSDStruct *dprot, cRMSDAlign &align, cRMSDMatrix &u, cRMSD3DCoord &xt, cRMSD3DCoord &yt);

/**
 * RMSD
 * 
 * \param qprot     query protein
 * \param dprot     database protein
 * \param align     alignment
 * \param stream    vmd script output stream
 */
inline double
rmsd(cRMSDStruct *qprot, cRMSDStruct *dprot, cRMSDAlign &align, std::ostream *stream = NULL)
{
        cRMSDMatrix u;
        cRMSD3DCoord xt;
        cRMSD3DCoord yt;
        
        double ret = rmsd(qprot, dprot, align, u, xt, yt);
        
        if(stream)
            generate_vmd(*stream, qprot, dprot, align, align, align, u, xt, yt);
        
        return ret;
}


#endif /* RMSD_HPP_ */
