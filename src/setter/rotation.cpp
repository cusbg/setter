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
#include "rotation.h"
#include <cmath>

using namespace std;

RMatrix::RMatrix() : cRMSDMatrix()
{
	// no-op
}

RMatrix::RMatrix(cRMSDMatrix &m) : cRMSDMatrix(m)
{
	// no=op
}

double &RMatrix::operator()(int i, int j)
{
	int ii, jj;
	if (i < 0) {
		ii = 0;
	}
	else if (i > 2) {
		ii = 2;
	}
	else {
		ii = i;
	}

	if (j < 0) {
		jj = 0;
	}
	else if (j > 2) {
		jj = 2;
	}
	else {
		jj = j;
	}

	return data[ii][jj];
}

double GetData(int i, int j)
{
	return 0;
}

RMatrix RMatrix::operator*(const RMatrix &other)
{
	RMatrix result;
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			double pRes = 0;
			for (int k = 0; k < 3; ++k) {
				pRes += data[i][k] * other.data[k].data[j];
			}
			result.data[i][j] = pRes;
		}
	}
	return result;
}

cRMSD3DCoord RMatrix::operator*(const cRMSD3DCoord &v)
{
	cRMSD3DCoord result;

	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			result[i] += data[i][j] * v.data[j];

	return result;
}

RMatrix RMatrix::Transpose()
{
	RMatrix result;
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			result.data[i][j] = data[j][i];
		}
	}
	return result;
}

/////// fontos

double RMatrix::GetTrace()
{
	return data[0][0] + data[1][1] + data[2][2];
}

double RMatrix::GetAngle()
{
	return acos((GetTrace() - 1) / 2.0);
}

cRMSD3DCoord RMatrix::GetUnitVector()
{
	double **a = new double *[3];
	for (int i = 0; i < 3; ++i) {
		a[i] = new double[3];
	}

	cRMSD3DCoord result;

	// copy matrix
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			if (abs(data[i][j]) < 0.0000001) {
				a[i][j] = 0.0;
			}
			else a[i][j] = data[i][j];
		}
	}

	int row;
	for (row = 0; row < 3; ++row) {
		if (a[row][row] != 0.0) {
			break;
		}
	} // find the row

	for (int i = 0; i < 3; ++i) {
		if (i == row) {
			continue;
		}
		if (a[i][row] == 0.0) {
			continue;
		}

		double mul = (a[row][row] - 1.0) / a[i][row];
		for (int j = 0; j < 3; ++j) {
			a[i][j] *= mul;
			a[i][j] -= a[row][j];
			a[i][j] /= mul;
		}
		a[i][row] = 0;
	} // null the rowth position in each rows

	// create new matrix
	double **b = new double *[2];
	for (int i = 0; i < 2; ++i) {
		b[i] = new double[2];
	}

	// copy the not null parts
	int map[2];
	int ii = 0;
	for (int i = 0; i < 3; ++i) {
		if (i == row) {
			continue;
		}
		map[ii] = i;
		int jj = 0;
		for (int j = 0; j < 3; ++j) {
			if (j == row) {
				continue;
			}
			b[ii][jj] = a[i][j];
			++jj;
		}
		++ii;
	}

	// compute with b
	if (b[0][1] == 0) {
		if (abs(b[0][0] - 1.0) < 0.0000001) {
			result[map[0]] = 1.0;
		}
		else {
			result[map[0]] = 0.0;
			if (abs(b[1][1] - 1.0) < 0.0000001) {
				result[map[1]] = 1.0;
			}
			else {
				result[map[1]] = 0.0;
			}
		}
	}
	else if (b[1][0] == 0) {
		if (abs(b[1][1] - 1.0) < 0.0000001) {
			result[map[1]] = 1.0;
		}
		else {
			result[map[1]] = 0.0;
			if (abs(b[0][0] - 1.0) < 0.0000001) {
				result[map[0]] = 1.0;
			}
			else {
				result[map[0]] = 0.0;
			}
		}
	}
	else {
		result[map[1]] = 1.0;
		result[map[0]] = (1 - b[1][1]) / b[1][0];
	}

	// delete b
	for (int i = 0; i < 2; ++i) {
		delete [] b[i];
	}
	delete [] b;

	// sompute result[row]
	double sum = 0.0;
	for (int i = 0; i < 3; ++i) {
		if (i == row) {
			continue;
		}
		sum += a[row][i] * result[i];
	}

	if (abs(sum) < 0.0000001) {
		if (abs(a[row][row] - 1.0) < 0.0000001) {
			result[row] = 1.0;
		}
		else {
			result[row] = 0.0;
		}
	}
	else {
		result[row] = sum / (1 - a[row][row]);
	}


	// delete a
	for (int i = 0; i < 3; ++i) {
		delete [] a[i];
	}
	delete [] a;

	// normalize
	double norm = sqrt(result[0] * result[0] + result[1] * result[1] + result[2] * result[2]);
	for (int i = 0; i < 3; i++) {
		result[i] /= norm;
	}

	return result;
}

RMatrix RMatrix::GetMatrix(double theta, cRMSD3DCoord u)
{
	RMatrix result;

	double sinT = sin(theta);
	double cosT = cos(theta);
	double icosT = 1 - cosT;


	result[0][0] = cosT + u[0] * u[0] * icosT;
	result[0][1] = u[0] * u[1] * icosT - u[2] * sinT;
	result[0][2] = u[0] * u[2] * icosT + u[1] * sinT;

	result[1][0] = u[1] * u[0] * icosT + u[2] * sinT;
	result[1][1] = cosT + u[1] * u[1] * icosT;
	result[1][2] = u[1] * u[2] * icosT - u[0] * sinT;

	result[2][0] = u[2] * u[0] * icosT - u[1] * sinT;
	result[2][1] = u[2] * u[1] * icosT + u[0] * sinT;
	result[2][2] = cosT + u[2] * u[2] * icosT;

	return result;
}

bool RMatrix::AllNull()
{
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			if (data[i][j] != 0.0) {
				return false;
			}
		}
	}
	return true;
}

RMatrix RMatrix::GetHalfRotation()
{
	if (AllNull()) {
		//cout << "null angle" << endl;
		return *this;
	}

	double angle = GetAngle();
	cRMSD3DCoord u = GetUnitVector();

	return GetMatrix(angle / 2.0, u);
}

// will be deleted
RotationData::RotationData()
{
	// no_op
}

RotationData::RotationData(RMatrix &rotM) :
mRotM(rotM)
{
	/*if (abs(mRotM(3,1)) != 1){
		mY1 = -asin(mRotM(3,1));
		mY2 = PI - mY1;

		mX1 = atan2(mRotM(3, 2) / cos(mY1), mRotM(3, 3) / cos(mY1));
		mX2 = atan2(mRotM(3, 2) / cos(mY2), mRotM(3, 3) / cos(mY2));

		mZ1 = atan2(mRotM(2, 1) / cos(mY1), mRotM(1, 1) / cos(mY1));
		mZ2 = atan2(mRotM(2, 1) / cos(mY2), mRotM(1, 1) / cos(mY2));
		} else {
		mZ1 = mZ2 = 0;

		if (mRotM(3, 1) == -1) {
		mY1 = mY2 = PI / 2;
		mX1 = mX2 = atan2(mRotM(1, 2), mRotM(1, 3));
		} else {
		mY1 = mY2 = -PI / 2;
		mX1 = mX2 = atan2(-mRotM(1, 2), -mRotM(1, 3));
		}
		}*/

	mX = mX1;
	mY = mY1;
	mZ = mZ1;
}

void RotationData::RotateX(double angle)
{
	mRX(0, 0) = 1;
	mRX(1, 1) = cos(angle);		mRX(1, 2) = -sin(angle);
	mRX(2, 1) = sin(angle);		mRX(2, 2) = cos(angle);
}
void RotationData::RotateY(double angle)
{
	mRY(0, 0) = cos(angle);		mRY(0, 2) = sin(angle);
	mRY(1, 1) = 1;
	mRY(2, 0) = -sin(angle);	mRY(2, 2) = cos(angle);
}
void RotationData::RotateZ(double angle)
{
	mRZ(0, 0) = cos(angle);		mRZ(0, 1) = -sin(angle);
	mRZ(1, 0) = sin(angle);		mRZ(1, 1) = cos(angle);
	mRZ(2, 2) = 1;
}

void RotationData::Rotate(double x, double y, double z, RMatrix &rot)
{
	RotateX(x); RotateY(y); RotateZ(z);
	rot = mRX * mRY * mRZ;
}

void RotationData::GetRotQ(RMatrix &rot)
{
	Rotate(-mX / 2, -mY / 2, -mZ / 2, rot);
}
void RotationData::GetRotDB(RMatrix &rot)
{
	Rotate(mX / 2, mY / 2, mZ / 2, rot);
}

void RotationData::Test()
{
	Rotate(0.1, 0.3, 0.0, this->mRotM);
	RotationData(this->mRotM);
}

RMatrix RotationData::GetMatrix()
{
	return mRotM;
}
