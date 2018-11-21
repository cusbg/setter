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
#include <cmath>
#include "align.hpp"
#include "vector.hpp"
#include "matrix.hpp"
#include "rmsd.hpp"

using namespace std;

/**
 * RMSD
 *
 * \param qprot     query protein
 * \param dprot     database protein
 * \param align     alignment
 * \param u         output RMSD rotation of the database protein (after its transposition)
 * \param xt        output RMSD transposition of the query protein
 * \param yt        output RMSD transposition of the database protein
 */
double
rmsd(cRMSDStruct *qprot, cRMSDStruct *dprot, cRMSDAlign &al, cRMSDMatrix &u, cRMSD3DCoord &xt, cRMSD3DCoord &yt)
{
	double sigma;

	cRMSDMatrix r;
	cRMSDMatrix a;
	cRMSDMatrix b;

	cRMSD3DCoord e;
	double d;
	double spur;
	double det;
	double cof;
	double p;
	double h;
	double g;
	double cth;
	double sth;
	double sqrth;

	double rr[6];
	double ss[6];

	double sqrt3 = 1.73205080756888E+00;
	double tol = 1.0E-2;
	double zero = 0.0E+00;
	double one = 1.0E+00;
	double two = 2.0E+00;
	double three = 3.0E+00;

	int ip[9] = { 0, 1, 3, 1, 2, 4, 3, 4, 5 };
	int ip2312[4] = { 1, 2, 0, 1 };

	double rms = 0.0;
	double e0 = zero;



	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			d = (i == j) ? one : zero;

			u[i][j] = d;
			a[i][j] = d;
			r[i][j] = zero;
		}
	}


	xt = cRMSD3DCoord();
	yt = cRMSD3DCoord();

#ifndef NOWEIGHT
	double w2 = 0;
#endif /*NOWEIGHT*/

	for (int m = 0; m < al.length(); m++)
	{
#ifdef NOWEIGHT
		xt = xt + qprot->coord[al[m].x];
		yt = yt + dprot->coord[al[m].y];
#else
		xt = xt + qprot->coord[al[m].x] * al[m].w * al[m].w;
		yt = yt + dprot->coord[al[m].y] * al[m].w * al[m].w;
		w2 += al[m].w * al[m].w;
#endif /*NOWEIGHT*/        
	}

#ifdef NOWEIGHT
	xt = xt / al.length();
	yt = yt / al.length();
#else    
	xt = xt / w2;
	yt = yt / w2;
#endif /*NOWEIGHT*/


	for (int m = 0; m < al.length(); m++)
	{
#ifdef NOWEIGHT
		e0 = e0 + (norm2(qprot->coord[al[m].x] - xt) + norm2(dprot->coord[al[m].y] - yt));
		r = r + ((qprot->coord[al[m].x] - xt) ^ (dprot->coord[al[m].y] - yt));
#else
		e0 = e0 + al[m].w * al[m].w * (norm2(qprot->coord[al[m].x] - xt) + norm2(dprot->coord[al[m].y] - yt));
		r = r   + al[m].w * al[m].w * ((qprot->coord[al[m].x] - xt) ^ (dprot->coord[al[m].y] - yt));
#endif /*NOWEIGHT*/
	}


	det = ((r[0][0] * ((r[1][1] * r[2][2]) - (r[1][2] * r[2][1]))) - (r[0][1] * ((r[1][0] * r[2][2]) - (r[1][2] * r[2][0])))) + (r[0][2] * ((r[1][0] * r[2][1]) - (r[1][1] * r[2][0])));

	sigma = det;

	int m = 0;

	for (int j = 0; j < 3; j++)
		for (int i = 0; i <= j; i++)
			rr[m++] = (r[0][i] * r[0][j]) + (r[1][i] * r[1][j]) + (r[2][i] * r[2][j]);

	spur = ((rr[0] + rr[2]) + rr[5]) / three;

	cof = ((((((rr[2] * rr[5]) - (rr[4] * rr[4])) + (rr[0] * rr[5])) - (rr[3] * rr[3])) + (rr[0] * rr[2])) - (rr[1] * rr[1])) / three;
	det = det * det;

	for (int i = 0; i < 3; i++)
		e[i] = spur;


	if (spur <= zero)
		goto j40;



	d = spur * spur;
	h = d - cof;

	g = (((spur * cof) - det) / two) - (spur * h);

	if (h > zero)
	{
		sqrth = sqrt(h);
		d = ((h * h) * h) - (g * g);

		if (d < zero)
			d = zero;

		d = atan2(sqrt(d), -g) / three;
		cth = sqrth * cos(d);
		sth = (sqrth * sqrt3) * sin(d);
		e[0] = (spur + cth) + cth;
		e[1] = (spur - cth) + sth;
		e[2] = (spur - cth) - sth;


		for (int l = 0; l < 3; l += 2)
		{
			d = e[l];

			ss[0] = ((d - rr[2]) * (d - rr[5])) - (rr[4] * rr[4]);
			ss[1] = ((d - rr[5]) * rr[1]) + (rr[3] * rr[4]);
			ss[2] = ((d - rr[0]) * (d - rr[5])) - (rr[3] * rr[3]);
			ss[3] = ((d - rr[2]) * rr[3]) + (rr[1] * rr[4]);
			ss[4] = ((d - rr[0]) * rr[4]) + (rr[1] * rr[3]);
			ss[5] = ((d - rr[0]) * (d - rr[2])) - (rr[1] * rr[1]);

			int j;

			if (abs(ss[0]) >= abs(ss[2]))
			{
				if (abs(ss[0]) < abs(ss[5]))
					j = 3;
				else
					j = 1;
			}
			else
			{
				if (abs(ss[2]) < abs(ss[5]))
					j = 3;
				else
					j = 2;
			}


			d = zero;

			j = 3 * (j - 1);

			for (int i = 0; i < 3; i++)
			{
				int k = ip[i + j];
				a[i][l] = ss[k];
				d = d + (ss[k] * ss[k]);
			}

			if (d > zero)
				d = one / sqrt(d);

			for (int i = 0; i < 3; i++)
				a[i][l] = a[i][l] * d;
		}



		d = ((a[0][0] * a[0][2]) + (a[1][0] * a[1][2])) + (a[2][0] * a[2][2]);


		int m;
		int m1;

		if ((e[0] - e[1]) > (e[1] - e[2]))
		{
			m1 = 2;
			m = 0;
		}
		else
		{
			m1 = 0;
			m = 2;
		}

		p = zero;


		for (int i = 0; i < 3; i++)
		{
			a[i][m1] = a[i][m1] - (d * a[i][m]);
			p = p + (a[i][m1] * a[i][m1]);
		}


		if (p > tol)
		{
			p = one / sqrt(p);

			for (int i = 0; i < 3; i++)
				a[i][m1] = a[i][m1] * p;
		}
		else
		{
			p = one;

			int j = 0;

			for (int i = 0; i < 3; i++)
			{
				if (p < abs(a[i][m]))
					continue; //?

				p = abs(a[i][m]);
				j = i;
			}


			int k = ip2312[j];
			int l = ip2312[j + 1];
			p = sqrt((a[k][m] * a[k][m]) + (a[l][m] * a[l][m]));

			if (p <= tol)
				goto j40;

			a[j][m1] = zero;
			a[k][m1] = -(a[l][m] / p);
			a[l][m1] = a[k][m] / p;
		}


		a[0][1] = (a[1][2] * a[2][0]) - (a[1][0] * a[2][2]);
		a[1][1] = (a[2][2] * a[0][0]) - (a[2][0] * a[0][2]);
		a[2][1] = (a[0][2] * a[1][0]) - (a[0][0] * a[1][2]);
	}



	for (int l = 0; l < 2; l++)
	{
		d = zero;

		for (int i = 0; i < 3; i++)
		{
			b[i][l] = ((r[i][0] * a[0][l]) + (r[i][1] * a[1][l])) + (r[i][2] * a[2][l]);
			d = d + (b[i][l] * b[i][l]);
		}

		if (d > zero)
			d = one / sqrt(d);

		for (int i = 0; i < 3; i++)
			b[i][l] *= d;
	}



	d = ((b[0][0] * b[0][1]) + (b[1][0] * b[1][1])) + (b[2][0] * b[2][1]);
	p = zero;

	for (int i = 0; i < 3; i++)
	{
		b[i][1] = b[i][1] - (d * b[i][0]);
		p = p + (b[i][1] * b[i][1]);
	}


	if (p > tol)
	{
		p = one / sqrt(p);

		for (int i = 0; i < 3; i++)
			b[i][1] = b[i][1] * p;
	}
	else
	{
		p = one;

		int j = 0;

		for (int i = 0; i < 3; i++)
		{
			if (p >= abs(b[i][0]))
			{
				p = abs(b[i][0]);
				j = i;
			}
		}

		int k = ip2312[j];
		int l = ip2312[j + 1];
		p = sqrt((b[k][0] * b[k][0]) + (b[l][0] * b[l][0]));

		if (p <= tol)
			goto j40;

		b[j][1] = zero;
		b[k][1] = -(b[l][0] / p);
		b[l][1] = b[k][0] / p;
	}

	b[0][2] = (b[1][0] * b[2][1]) - (b[1][1] * b[2][0]);
	b[1][2] = (b[2][0] * b[0][1]) - (b[2][1] * b[0][0]);
	b[2][2] = (b[0][0] * b[1][1]) - (b[0][1] * b[1][0]);



	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			u[i][j] = ((b[i][0] * a[j][0]) + (b[i][1] * a[j][1])) + (b[i][2] * a[j][2]);


j40:


	for (int i = 0; i < 3; i++)
	{
		if (e[i] < zero)
			e[i] = zero;
		else
			e[i] = sqrt(e[i]);
	}

	d = e[2];

	if (sigma < 0.0)
		d = -d;


	d = (d + e[1]) + e[0];
	rms = (e0 - d) - d;

	if (rms < 0.0)
		rms = 0.0;


	return rms;
}
