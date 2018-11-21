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
#ifndef ALIGN_HPP_
#define ALIGN_HPP_

#include "protein.hpp"
#include "vector.hpp"

#define NOWEIGHT

/**
 * Weight of the distance of the resid pair.
 * 
 * \param qprot         query protein
 * \param dprot         database protein
 * \param qres          resid in the query protein
 * \param dres          resid in the database protein
 */
inline double
pair_weight(cRMSDStruct *qprot, cRMSDStruct *dprot, int qres, int dres) {
    return 1;
}



/**
 * Alignment item
 */
typedef struct
{
    int x;      // resid number in query protein (the first has 0)
    int y;      // resid number in database protein (the first has 0)
#ifndef NOWEIGHT
    double w;   // weight
#endif /*NOWEIGHT*/    
}
align_t;



inline bool
operator!=(align_t a, align_t b)
{
    if(a.x != b.x || a.y != b.y)
        return true;
    else
        return false;
}



/**
 * Alignment
 */
class cRMSDAlign
{
    align_t *al;    // data
    int len;        // length
    int *counter;   // data reference counter
    
public:
    
    /**
     * Create an empty alignment with the given length. 
     */
    cRMSDAlign(int length=0) : al(0), counter(0)
    {
        len = length;
        
        if(len)
        {
            al = new align_t[len];
            counter = new int(1);
        }
    }
    
    
    /**
     * Create a shallow copy.
     */
    cRMSDAlign(const cRMSDAlign &align)
    {
        al      = align.al;
        len     = align.len;
        counter = align.counter;
        
        if(counter)
            (*counter)++;

    }
    
    
    /**
     * Create an alignment from a DP path. 
     */
    cRMSDAlign(int plen, int (*path)[2]);
    
    
    /**
     * Destructor.
     */
    ~cRMSDAlign()
    {
        if(counter)
        {
            (*counter)--;
            
            if(! *counter)
            {
                delete counter;
                delete[] al;
            }
        }
    }
    
    
    /**
     * Assignment operator (shallow copy)
     */
    cRMSDAlign &operator=(const cRMSDAlign &align)
    {
        if(counter && al == align.al)
            return *this;
        
        if(counter)
        {
            (*counter)--;
            
            if(! *counter)
            {
                delete counter;
                delete[] al;
            }
        }
        
        al      = align.al;
        len     = align.len;
        counter = align.counter;
        
        if(counter)
            (*counter)++;
        
        return *this;
    }
    
    
    /**
     * Fill the alignment according to a part of other alignment.
     */
    void
    set(cRMSDAlign &base, int from, int to)
    {
        len = to - from;
        
        for(int i = from; i < to; i++)
        {
            al[i-from].x = base.al[i].x;
            al[i-from].y = base.al[i].y;
#ifndef NOWEIGHT
            al[i-from].w = base.al[i].w;
#endif /*NOWEIGHT*/    

        }
    }
    
    
    /**
     * Add to the alignment a part of other alignment.
     */
    void
    add(cRMSDAlign &base, int from, int to)
    {
        int offset = len;
        
        len += to - from;
        
        for(int i = from; i < to; i++)
        {
            al[offset + i-from].x = base.al[i].x;
            al[offset + i-from].y = base.al[i].y;
#ifndef NOWEIGHT
            al[offset + i-from].w = base.al[i].w;
#endif /*NOWEIGHT*/    
        }
    }
    
    
    /**
     * Returns the length of the alignment.
     */
    int
    length()
    {
        return len;
    }
    
    
    /**
     * Set length of the alignment.
     */
    void
    setlength(int len)
    {
        this->len = len;
    }
    
    
    align_t &
    operator[](int i)
    {
        return al[i];
    }
    
    
    bool
    operator==(cRMSDAlign &align)
    {
        return al == align.al;
    }
    
    
    cRMSDAlign
    copy()
    {
        cRMSDAlign align(len);
        align.set(*this, 0, len);
        return align;
    }
    
    
    /**
     * Return an alignment after DP phase.
     */
    cRMSDAlign dp(cRMSDStruct *cprot, cRMSDStruct *dprot, int dpwidth, double d0);
};  


#endif /* ALIGN_HPP_ */
