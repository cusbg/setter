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
#ifndef PROTEIN_HPP_
#define PROTEIN_HPP_

#include "vector.hpp"
#include "matrix.hpp"

#include <vector>



/**
 * The cRMSDStruct class represents a protein structure.
 */
class cRMSDStruct
{
public:
	std::string name;           // name of the protein (need for generate_vmd)
    
    unsigned int length;                 // number of resids (max length)
	unsigned int lengthTemp;             // number of resids (lengthTemp <= length)
    cRMSD3DCoord *coord;              // coordinates of resids
    std::vector<std::string> aa;
	std::string sse;                  // type of secondary structure (A,B,T)
    
    int resid_offset;           // number of first resid in PDB (need for generate_vmd)
    
    int family;
    int superfamily;
    int fold;
    int klass;
    
    cRMSDStruct()
    {
		//name = new char[100]; 
		length = 0;
		coord = NULL;
    }
    
    /**
     * Creates a new transformed protein structure.
     * 
     * \param p     template protein
     * \param u     rotation
     * \param t     transposition
     */
    cRMSDStruct(cRMSDStruct *p, const cRMSDMatrix &u, const cRMSD3DCoord &t)
    {
		//name = new char[100]; 

		SetName(p->name);
		//strcpy(name, p->name); 
        length = p->length;
		lengthTemp = p->lengthTemp;
        aa = p->aa;
		sse = p->sse;
        resid_offset = p->superfamily;
        superfamily = p->superfamily;
        fold = p->fold;
        klass = p->klass;
        
        coord = new cRMSD3DCoord[length];
        
        for(unsigned int i = 0; i < length; i++)
            coord[i] = u*p->coord[i] + t;
    }	

	void Reinit()
	{
		if (coord != NULL)
		{
			delete[] coord;
			coord = NULL;
			length = 0;
			name = "";
		}
	}
    
    ~cRMSDStruct()
    {
		if (coord != NULL)
		{
			delete[] coord;			
		}
    }

	void SetName(std::string _name)
	{
		name = _name;
	}
    
    /*
    int
    operator()(cRMSD3DCoord &v)
    {
        int res = -1;
        double best = 0;
        
        for(int i = 0; i < length; i++)
        {
            double d = norm2(coord[i] - v);
            
            if(d < best || res == -1)
            {
                res = i;
                best = d;
            }
        }
        
        return res;
    }
    */
};


#endif /* PROTEIN_HPP_ */
