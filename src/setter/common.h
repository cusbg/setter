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

#include <vector>
#include <map>
#include <string>
#include <sstream>

#include <boost/date_time.hpp>

#include "cRNAStructure.h"
#include "distance/rmsd.hpp"

using namespace std;

struct GSSUPair
{
	//! compairs two GSSU pairs according the q
	static bool compareGSSUPair(GSSUPair i1, GSSUPair i2);
	
	GSSUPair();
	//! Construct GSSU pair.
	GSSUPair(int q, int db);

	int q;
	int db;
};

typedef int hairpinIdx;

typedef map<hairpinIdx, cRNAStructure> RNAMap;

typedef vector<cRNAStructure> RNAVect;

//! One type of RNA could behave multiple ways according the distance between rnaQ and rnaDB.
enum InputType
{
	Q,		/*!< query RNA as query RNA */
	DB,		/*!< db RNA as db RNA */
	NO		/*!< query RNA as db RNA */
};

class ScoreRes;
class cRMSDStructExt;

//! Structure copied from SETTER
typedef struct ResMatch
{
	// the maximum value of double
	static const double MAX;

	// default constructor
	ResMatch();

	double												score;
	cRMSDMatrix											rot; //rotation matrix
	cRMSD3DCoord										trans; //translation vector	
	vector< pair < pair < string, string >, double > >	nn[2]; // set of NN, for each residue in Q, its NN in DB and distance
	vector<pair<int, int> >								alignedGSSUs;
	vector<double>										alignedScores;
	vector<double>										rotationPenalties;
	tPair												superpositionTriplet[3]; //triplet used for the overall superposition

	boost::posix_time::ptime							timeStart, timeEnd; //start and end time of the matching process

} sResMatch;

//! Structure copied from SETTER
typedef struct Params
{
	int		neckMaxShift; //8 - number of residues under the neck which will be tried as neck residues (preventing
	double	identicalLetterModificator; //0.1 - if residues having the same nt type are aligned their distance is multiplied by given factor
	double	identicalPairTypeModificator; //1 - if residues being in the same SSE type are aligned their distance is multiplied by given factor
	double	NNDist; //; 6 //angstroms - the distance threshold for which the aligned residues are being counted in
	int		topKForMultiGSSU;	//3 - all-to-all GSSU distances are computed and topK best (with minimum distance) pairs are taken into account for full structure alignemnt
	int		noHeadShiftStepRatioDiviser; //10 - if a structure does not have secondary structure or its GSSU does not contain head, only GSSU_length/noHeadShiftStepRatioDiviser are tried when aligning
	int		earlyTermination; //1 - (lambda) influences how often early termination is applied. The lower value, the more often.
	double  rotationPenalty; //0.9
} sParams;

sParams GetDefaultParams();
sParams GetGlobalParams();
void SetGlobalParams(sParams params);

//! Global params structure to be shared by all parts of the application.
extern sParams gParams;



//! Class copied from SETTER
class ScoreRes
{
	int mIdxQ;
	int mIdxDB;
	double mScore;
	cRMSDMatrix mRot;
	cRMSD3DCoord mTrans;

public:
	static bool Compare(ScoreRes o1, ScoreRes o2);

	ScoreRes(int idxQ, int idxDB, double score, cRMSDMatrix rot, cRMSD3DCoord trans);
	pair<int, int> GetPair();
	double GetScore();
	cRMSDMatrix &GetRot();
	cRMSD3DCoord &GetTrans();
};

//! Class copied from SETTER
class cRMSDStructExt : public cRMSDStruct
{
public:
	vector<tPairResiude> residueInfo;

	cRMSDStructExt();
	cRMSDStructExt(cRMSDStruct *p, const cRMSDMatrix &u, const cRMSD3DCoord &t);
};

//! Calculates the average of coordinates
/*!
	\param atom1	Input atom.
	\param atom2	Input atom.
	\param aver		Output coordinates calculated from (atom1.coords + atom2.coords)/2.


*/
void CalcAver(tPdbAtom &atom1, tPdbAtom &atom2, t3DCoords &aver);

//! Calculates the translation from new coordinates and previous coordinates
/*!
	\param atom			The atom before translation.
	\param newCoord		The new coordinates of the atom.
	\param trans		The translation vectore which translates the atom to new coordinates.
*/
void CalcTrans(tPdbAtom &atom, t3DCoords &newCoord, t3DCoords &trans);

//! Calculates the score (minimum, aritmetical mean, geometrical mean)
double CalcFinalScore(double score1, double scre2);

//! Calculates score, matches rna structures.
sResMatch CalcScore(
	cRNAStructure &rnaQ, 
	cRNAStructure &rnaDB, 
	sParams &params);

//! returns the PDB ID of the RNA structure.
string GetPdbId(cRNAStructure &rna);

//! Returns true if the file exists with path given as parameter.
bool FileExists(string file);

//! Translate coordinates given a rotation matrix and translation vector.
void TranslateCoord(t3DCoords &coord, const cRMSDMatrix &u, const cRMSD3DCoord &t);

//! Returns current time.
string GetTime();

//! Returns current date.
string GetDate();

//! Removes non-digit characters from the input string (e.g. 32A -> 32)
/*!
	\param number		String containing the number to be cleansed.
*/
string CleanseNumber(string number);

vector<string> GetRGBColors(int n);
