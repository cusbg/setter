/*
 Copyright (c) 2013 David Hoksza, Peter Szepe

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
#pragma once

#include "distance/matrix.hpp"
#include <math.h>

#define PI 3.14159265


class RMatrix : cRMSDMatrix
{
public:
	double GetTrace();
	double GetAngle();
	cRMSD3DCoord GetUnitVector();
	RMatrix GetMatrix(double theta, cRMSD3DCoord unitVect);
	RMatrix Transpose();
	bool AllNull();
	RMatrix GetHalfRotation();

public:
	RMatrix();
	RMatrix(cRMSDMatrix &m);

	double &operator()(int i, int j);
	double GetData(int i, int j);
	RMatrix operator*(const RMatrix &other);
	cRMSD3DCoord operator*(const cRMSD3DCoord &v);
};

class RotationData
{
	RMatrix mRotM;
	double mX1, mX2;
	double mY1, mY2;
	double mZ1, mZ2;

	double mX, mY, mZ;

	RMatrix mRX;
	RMatrix mRY;
	RMatrix mRZ;

	void RotateX(double angle);
	void RotateY(double angle);
	void RotateZ(double angle);
public:
	RotationData();
	RotationData(RMatrix &rotM);

	void Rotate(double x, double y, double z, RMatrix &res);
	void GetRotQ(RMatrix &rot);
	void GetRotDB(RMatrix &rot);

	RMatrix GetMatrix();

	void Test();
};