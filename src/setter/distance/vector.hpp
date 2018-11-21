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
#ifndef VECTOR_HPP_
#define VECTOR_HPP_

#include <iostream>



/**
 * The cRMSD3DCoord class represents a 3-dimensional vector.
 */
class cRMSD3DCoord
{
	friend class RMatrix;
protected:
    double data[3];
    
public:
    cRMSD3DCoord()
    {
        for(int i = 0; i < 3; i++)
            data[i] = 0.0;
    }
    
    cRMSD3DCoord(double value[3])
    {
        for(int i = 0; i < 3; i++)
            data[i] = value[i];
    }
    
    cRMSD3DCoord(double x, double y, double z)
    {
        data[0] = x;
        data[1] = y;
        data[2] = z;
    }
    
    cRMSD3DCoord(const cRMSD3DCoord &v)
    {
        for(int i = 0; i < 3; i++)
            data[i] = v.data[i];
    }
    
    double &
    operator[](int index)
    {
        return data[index];
    }

};



inline cRMSD3DCoord
operator+(cRMSD3DCoord a, cRMSD3DCoord b)
{
    cRMSD3DCoord v;
    
    for(int i = 0; i < 3; i++)
        v[i] = a[i] + b[i];
    
    return v;
}



inline cRMSD3DCoord
operator-(cRMSD3DCoord a, cRMSD3DCoord b)
{
    cRMSD3DCoord v;
    
    for(int i = 0; i < 3; i++)
        v[i] = a[i] - b[i];
    
    return v;
}



inline double
operator*(cRMSD3DCoord a, cRMSD3DCoord b)
{
    double d = 0.0;
    
    for(int i = 0; i < 3; i++)
        d += a[i] * b[i];
    
    return d;
}



inline cRMSD3DCoord
operator*(double a, cRMSD3DCoord b)
{
    cRMSD3DCoord v;
    
    for(int i = 0; i < 3; i++)
        v[i] = a * b[i];
    
    return v;
}



inline cRMSD3DCoord
operator*(cRMSD3DCoord a, double b)
{
    cRMSD3DCoord v;
    
    for(int i = 0; i < 3; i++)
        v[i] = a[i] * b;
    
    return v;
}



inline cRMSD3DCoord
operator/(cRMSD3DCoord a, double b)
{
    cRMSD3DCoord v;
    
    for(int i = 0; i < 3; i++)
        v[i] = a[i] / b;
    
    return v;
}



inline std::ostream &
operator<<(std::ostream &stream, cRMSD3DCoord a)
{
    stream << "{ ";
    
    for(int i = 0; i < 3; i++)
        stream << a[i] << " ";
    
    stream << "}";
    return stream;
}



/**
 * Computes the square of the norm of the vector.
 */
inline double
norm2(cRMSD3DCoord a)
{
    return a*a;
}


#endif /* VECTOR_HPP_ */
