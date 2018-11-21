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
#ifndef MATRIX_HPP_
#define MATRIX_HPP_

#include "vector.hpp"

/**
 * The cRMSDMatrix class represents a 3-by-3 matrix. 
 */
class cRMSDMatrix
{
protected:
    cRMSD3DCoord data[3];

public:
    cRMSDMatrix()
    {
        for(int i = 0; i < 3; i++)
            for(int j = 0; j < 3; j++)
                data[i][j] = 0.0;
    }
    
    cRMSDMatrix(double value[3][3])
    {
        for(int i = 0; i < 3; i++)
            for(int j = 0; j < 3; j++)
                data[i][j] = value[i][j];
    }
    
    cRMSDMatrix(const cRMSDMatrix &v)
    {
        for(int i = 0; i < 3; i++)
            data[i] = v.data[i];
    }
    
    cRMSD3DCoord &
    operator[](int index)
    {
        return data[index];
    }
};



inline cRMSDMatrix
operator+(cRMSDMatrix a, cRMSDMatrix b)
{
    cRMSDMatrix m;
    
    for(int i = 0; i < 3; i++)
        m[i] = a[i] + b[i];
    
    return m;
}



inline cRMSDMatrix
operator-(cRMSDMatrix a, cRMSDMatrix b)
{
    cRMSDMatrix m;
    
    for(int i = 0; i < 3; i++)
        m[i] = a[i] - b[i];
    
    return m;
}

inline double
operator*(cRMSDMatrix a, cRMSDMatrix b)
{
    double product = 0;
    
    for(int i = 0; i < 3; i++)
        product += a[i]*b[i];
    
    return product;
}



inline cRMSDMatrix
operator*(double a, cRMSDMatrix b)
{
    cRMSDMatrix m;
    
    for(int i = 0; i < 3; i++)
        m[i] = a * b[i];
    
    return m;
}



inline cRMSDMatrix
operator*(cRMSDMatrix a, double b)
{
    cRMSDMatrix m;
    
    for(int i = 0; i < 3; i++)
        m[i] = a[i] * b;
    
    return m;
}



inline cRMSD3DCoord
operator*(cRMSDMatrix a, cRMSD3DCoord b)
{
    cRMSD3DCoord v;
    
    for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
            v[i] += a[i][j] * b[j];
    
    return v;
}



/**
 * Computes the vector direct product of vectors.
 */
inline cRMSDMatrix
operator^(cRMSD3DCoord a, cRMSD3DCoord b)
{
    cRMSDMatrix m;
    
    for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
            m[i][j] = a[i] * b[j];
    
    return m;
}



inline cRMSDMatrix
operator/(cRMSDMatrix a, double b)
{
    cRMSDMatrix m;
    
    for(int i = 0; i < 3; i++)
        m[i] = a[i] / b;
    
    return m;
}


inline double
norm2(cRMSDMatrix a)
{
	return norm2(a[0] + a[1] + a[2])/3;
}

#endif /* MATRIX_HPP_ */
