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
#pragma once

#include "common.h"

void ConvertHairpinToRNA(
	cRNAStructure &rna, 
	cHairpin &hp, 
	cRMSDStructExt &s, 
	int alignActual[3][2], 
	int alignActualConverted[3]);

pair<int, double> IdentifyL2NearestNeighbor(
	cRMSDStructExt &s, 
	cRMSD3DCoord &c, 
	int ixFrom, 
	int ixTo);

pair<int, double> IdentifyL2NearestNeighbor(
	cHairpin &hp, 
	cRNAStructure &rna, 
	cRMSD3DCoord &c);

void SortArray3(int arr[]);

/*
 * cntNN ... number of NN in a given range
*/
double Normalize(
	double score, 
	int lenQ, 
	int lenDB, 
	int cntHpQ, 
	int cntHpDB, 
	int cntNN);

sResMatch MatchHairpinsBasedOnTriplets(
	cRMSDStructExt &sQ, 
	cRNAStructure &rnaQ, 
	cHairpin &hpQ, 
	cRMSDStructExt &sDB, 
	cRNAStructure &rnaDB, 
	cHairpin &hpDB, 
	sParams params, 
	cRMSDAlign &align, 
	int alignActualQ[3][2], 
	int alignActualDB[3][2], 
	double scoreLimit);

sResMatch MatchHairpins(
	cRNAStructure &rnaQ, 
	cHairpin &hpQ, 
	cRNAStructure &rnaDB, 
	cHairpin &hpDB, 
	sParams params, 
	double scoreLimit = ResMatch::MAX);

typedef struct IxsScore
{
	IxsScore(double xscore, int xi, int xj);

	double score;
	int i, j;
} sIxsScore;

bool compareIxsScore(IxsScore i1, IxsScore i2);

sResMatch Match(
	cRNAStructure rnaQ, 
	cRNAStructure rnaDB, 
	sParams params, 
	bool extraSuperposition = false,
	bool nnListNeeded = false, 
	double bestScore = ResMatch::MAX);

sResMatch Match(
	cRNAStructure rnaQ, 
	cRNAStructure rnaDB, 
	sParams params, 
	sResMatch ***outRess, 
	vector<GSSUPair> &outOptPair, 
	bool nnListNeeded = false, 
	double bestScore = ResMatch::MAX);

double Match(
	cRNAStructure rnaQ, 
	string chainQ, 
	cRNAStructure rnaDB, 
	string chainDB, 
	sParams params, 
	double scoreLimit = ResMatch::MAX);