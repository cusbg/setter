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

#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <algorithm>
#include "align.hpp"
#include "protein.hpp"




cRMSDAlign::cRMSDAlign(int length, int (*path)[2])
{
    len = length;
    al = new align_t[len];
    counter = new int(1);
    
    int alen = 0;
    int prev_x = std::max(path[0][1]-1, 0);
    int prev_y = std::max(path[0][0]-1, 0);
    
    
    for(int i=0; i < length; i++)
    {
        if((path[i][1] == prev_x + 1 && path[i][0] == prev_y + 1))
        {
            (*this)[alen].x = prev_x;
            (*this)[alen].y = prev_y;
#ifndef NOWEIGHT
            (*this)[alen].w = 1;
#endif /*NOWEIGHT*/
            alen++;
        }
        
        prev_x = path[i][1];
        prev_y = path[i][0];
    }
    
    setlength(alen);
}



/**
 * Local DP. Perform DP for a fragment of the quuery protein and a fragment of the quuery protein.
 */
inline void
local_dp(int qprot_min, int qprot_max, int dprot_min, int dprot_max, double *score, double *table, /*bool *gap,*/ double d0, cRMSDStruct *qprot,  cRMSDStruct *dprot)
{
    int sr = dprot->length;
    int tr = dprot->length+1;
    
    for(int i = qprot_min; i <= qprot_max; i++)
    {
        for(int j = dprot_min; j <= dprot_max; j++)
        {
            // score!
            cRMSD3DCoord diff = qprot->coord[i] - dprot->coord[j];
#ifdef NOWEIGHT
            double dis = sqrt(diff * diff);
#else
            double dis = sqrt(diff * diff) * pair_weight(qprot, dprot, i, j);
#endif /*NOWEIGHT*/
            score[i*sr+j] = 1 / (1 + ((dis / d0) * (dis / d0)));    // item score!
            
            
            double va = score[i*sr+j] + table[i*tr+j]; 
            double vt = std::max(table[(i+1)*tr+j] /*+ (gap[(i+1)*tr+j]?EGP:OGP)*/, table[i*tr+j+1] /*+ (gap[i*tr+j+1]?EGP:OGP)*/);
            
            if(va >= vt)
            {
                table[(i+1)*tr+j+1] = va;
                //gap[(i+1)*tr+j+1] = false;
            }
            else
            {
                table[(i+1)*tr+j+1] = vt;
                //gap[(i+1)*tr+j+1] = true;
            }
        }
    }
}



/**
 * DP
 */
cRMSDAlign
cRMSDAlign::dp(cRMSDStruct *qprot, cRMSDStruct *dprot, int dpwidth, double d0)
{
    cRMSDAlign tmpalign(std::min(qprot->length, dprot->length));
    cRMSDAlign align = *this;
    
    double *score = new double[qprot->length * dprot->length];
    int sr = dprot->length;
    
    //bool *gap = new bool[(qprot->length + 1)*(dprot->length + 1)];
    
    double *table = new double[(qprot->length + 1)*(dprot->length + 1)];
    int tr = dprot->length + 1;
    

    for(unsigned int i = 0; i < (qprot->length + 1)*(dprot->length + 1); i++)
    {
        //gap[i] = false;
        table[i] = -1000;
    }
    
    for(unsigned int i = 0; i <= qprot->length; i++)
    {
        //gap[i*tr] = true;
        table[i*tr] = 0;
    }
    
    for(unsigned int j = 1; j <= dprot->length; j++)
    {
        //gap[j] = true;
        table[j] = 0;
    }
    
    
    if(dpwidth > 0)
    {
        int pre_x = align[0].x - 1;
        int pre_y = align[0].y - 1;
    
        
        if((*this)[0].x != 0 || (*this)[0].x != 0)
        {
            int i_R = std::min<int>(align[0].x + dpwidth, qprot->length - 1);
            int dprot_max = std::min<int>(align[0].y + dpwidth, dprot->length - 1);
            
            local_dp(0, i_R, 0, dprot_max, score, table, /*gap,*/ d0, qprot, dprot);
        }
        
        for (int p = 0; p < align.length(); p++)
        {
            if(align[p].x == pre_x + 1 && align[p].y == pre_y + 1)
            {
                int i = align[p].x;
                int dprot_min = std::max<int>(align[p].y - dpwidth, 0);
                int dprot_max = std::min<int>(align[p].y + dpwidth, dprot->length - 1);
                
                
                local_dp(i, i, dprot_min, dprot_max, score, table, /*gap,*/ d0, qprot, dprot);
            }
            else
            {
                int i_L = std::max<int>(pre_x - dpwidth, 0);
                int i_R = std::min<int>(align[p].x + dpwidth, qprot->length - 1);
                
                int dprot_min = std::max<int>(pre_y - dpwidth, 0);
                int dprot_max = std::min<int>(align[p].y + dpwidth, dprot->length - 1);
                
                
                local_dp(i_L, pre_x-1, pre_y, dprot_max, score, table, /*gap,*/ d0, qprot, dprot);
                local_dp(pre_x,align[p].x, dprot_min, dprot_max, score, table, /*gap,*/ d0, qprot, dprot);
                local_dp(align[p].x+1, i_R, dprot_min, align[p].y, score, table, /*gap,*/ d0, qprot, dprot);
            }
            
            pre_x = align[p].x;
            pre_y = align[p].y;
        }
        
        
        
        if(align[align.length()-1].x != qprot->length - 1 || align[align.length()-1].y != dprot->length - 1)
        {
            int i_L = std::max(align[align.length()-1].x - dpwidth, 0);
            int dprot_min = std::max(align[align.length()-1].y - dpwidth, 0);
            
            local_dp(i_L, qprot->length-1, dprot_min, dprot->length-1, score, table, /*gap,*/ d0, qprot, dprot);
        }
    }
    else
    {
        local_dp(0, qprot->length-1, 0, dprot->length-1, score, table, /*gap,*/ d0, qprot, dprot);
    }
    
    
    /*********************************/
    int i = qprot->length - 1;
    int j = dprot->length - 1;
    int len = 0;
    
    while(i != -1 && j != -1)
    {
        if(table[(i+1)*tr + j+1] == table[i*tr+j] + score[i*sr+j])
        {
            tmpalign[len].x = i;
            tmpalign[len].y = j;
            len++;
            
            i--;
            j--;
        }
        else
        {
            if(table[(i+1)*tr + j+1] == table[i*tr+j+1]/* + (gap[i*tr+j+1]?EGP:OGP)*/)
                i--;
            else if(table[(i+1)*tr + j+1] == table[(i+1)*tr+j]/* + (gap[(i+1)*tr+j]?EGP:OGP)*/)
                j--;
            else
                exit(2);
        }
    }
    
    
    bool state = true;
    
    if(this->len != len)
        state = false;
    
    cRMSDAlign newalign(std::min(qprot->length, dprot->length));
    for(int l=0; l < len; l++)
    {
        newalign[l].x = tmpalign[len-1-l].x;
        newalign[l].y = tmpalign[len-1-l].y;
		#ifndef NOWEIGHT
        newalign[l].w = pair_weight(qprot, dprot, newalign[l].x, newalign[l].y);
		#endif /*NOWEIGHT*/
        
        if(l >= this->len || (*this)[l].x != newalign[l].x || (*this)[l].y != newalign[l].y)
            state = false;
    }
    
    newalign.setlength(len);
    
    
    //delete[] gap;
    delete[] table;
    delete[] score;
    
    if(state)
        return *this;
    else
        return newalign;
}
