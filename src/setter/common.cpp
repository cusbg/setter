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
#include "common.h"
#include "algo.h"
#include "parameters.h"
#include <limits>
#include <math.h>

using namespace std;

sParams gParams;

GSSUPair::GSSUPair() : q(-1), db(-1)
{
	// no-op
}

GSSUPair::GSSUPair(int q, int db) : q(q), db(db)
{
	// no-op
}

bool GSSUPair::compareGSSUPair(GSSUPair i1, GSSUPair i2)
{
	return (i1.q < i2.q);
}

double const ResMatch::MAX = numeric_limits<double>::max();

ResMatch::ResMatch() {
	score = numeric_limits<double>::max();
}

bool ScoreRes::Compare(ScoreRes o1, ScoreRes o2)
{
	// TODO something with limits
	if (o1.mScore == std::numeric_limits<double>::min() && o2.mScore == std::numeric_limits<double>::min()) {
		return true;
	} else if (o1.mScore == std::numeric_limits<double>::max() && o2.mScore == std::numeric_limits<double>::max()) {
		return true;
	} else {
		return o1.mScore <= o2.mScore;
	}
}

ScoreRes::ScoreRes(int idxQ, int idxDB, double score, cRMSDMatrix rot, cRMSD3DCoord trans) :
	mIdxQ(idxQ), mIdxDB(idxDB), mScore(score), mRot(rot), mTrans(trans)
{
	// no-op
}

cRMSDMatrix &ScoreRes::GetRot()
{
	return mRot;
}

cRMSD3DCoord &ScoreRes::GetTrans()
{
	return mTrans;
}

pair<int, int> ScoreRes::GetPair()
{
	return make_pair(mIdxQ, mIdxDB);
}

double ScoreRes::GetScore()
{
	return mScore;
}

cRMSDStructExt::cRMSDStructExt() : cRMSDStruct() 
{
	// no-op
}
cRMSDStructExt::cRMSDStructExt(
	cRMSDStruct *p, 
	const cRMSDMatrix &u, 
	const cRMSD3DCoord &t) : cRMSDStruct(p, u, t) 
{
	// no-op
}

void CalcAver(tPdbAtom &atom1, tPdbAtom &atom2, t3DCoords &aver)
{
	aver.x = (atom1.coords.x + atom2.coords.x) / 2;
	aver.y = (atom1.coords.y + atom2.coords.y) / 2;
	aver.z = (atom1.coords.z + atom2.coords.z) / 2;
}

void CalcTrans(tPdbAtom &atom, t3DCoords &newCoord, t3DCoords &trans)
{
	trans.x = newCoord.x - atom.coords.x;
	trans.y = newCoord.y - atom.coords.y;
	trans.z = newCoord.z - atom.coords.z;
}

/// calculates the final score from two scores
double CalcFinalScore(double score1, double score2)
{
	//cout << "score_compare:" << score1 << ", " << score2 << endl;
	double result;

	Parameters &param = Parameters::GetInstance();

	if (param.GetMean() == Parameters::Mean::MIN) {
		result =  min(score1, score2);
	}

	if (param.GetMean() == Parameters::Mean::ARITHM) {
		result = (score1 + score2) / 2.0;
	}

	if (param.GetMean() == Parameters::Mean::GEOM) {
		result = sqrt(score1 * score2);
	}
	
	return result;
}

sResMatch CalcScore(
	cRNAStructure &rnaQ, 
	cRNAStructure &rnaDB, 
	sParams &params) 
{
	sResMatch res = Match(rnaQ, rnaDB, params);
	double score2 = Match(rnaDB, rnaQ, params).score;

	//cout << "score_compare:" << score1 << ", " << score2 << endl;

	res.score = CalcFinalScore(res.score, score2);

	return res;
}

string GetPdbId(cRNAStructure &rna)
{
	string id = rna.getId();

	int pos = id.find_last_of("/");
	/*if (pos >= 0) id = id.substr(pos + 1);

	pos = id.find_last_of(".pdb");
	if (pos == id.size()-1){
		id = id.substr(0, pos - 3);
	}

	return id;*/
	if (pos == -1) {
		return id;
	}

	return id.substr(pos + 1, 4);
}

bool FileExists(string file) {
	FILE *fp = fopen(file.c_str(),"r");
	if( fp ) {
		fclose(fp);
		return true;
	} else {
		return false;
	}
}

string GetDate()
{
	char cptime[50]; 
	time_t now = time(NULL);
	strftime(cptime, 50, "%Y-%m-%d", localtime(&now)); //uses short month name
	string strTime = cptime;
	return strTime;
}

string GetTime()
{
	char cptime[50]; 
	time_t now = time(NULL);
	strftime(cptime, 50, "%I:%M:%S", localtime(&now));
	string strTime = cptime;
	return strTime;
}

void TranslateCoord(t3DCoords &coord, const cRMSDMatrix &u, const cRMSD3DCoord &t)
{
	cRMSD3DCoord coordRMSD = u * cRMSD3DCoord(coord.x, coord.y, coord.z) + t;
	coord.x = coordRMSD[0];
	coord.y = coordRMSD[1];
	coord.z = coordRMSD[2];
}

string CleanseNumber(string number)
{
	string cleanString;
	for(string::iterator it = number.begin(); it != number.end(); ++it)
		if (isdigit(*it)) cleanString.push_back(*it);
	return cleanString;
}

sParams GetDefaultParams()
{
	sParams params;

	params.neckMaxShift = 8;
	params.identicalLetterModificator = 0.1;
	params.identicalPairTypeModificator = 1;
	params.NNDist = 6;
	params.topKForMultiGSSU = 3;
	params.noHeadShiftStepRatioDiviser = 10;
	params.earlyTermination = 1;
	params.rotationPenalty = 0.9;

	return params;
}

void SetGlobalParams(sParams params)
{
	gParams = params;
}

sParams GetGlobalParams()
{
	return gParams;
}

// This is a subfunction of HSLtoRGB
static void HSLtoRGB_Subfunction(unsigned int& c, const float& temp1, const float& temp2, const float& temp3)
{
	if((temp3 * 6) < 1)
		c = (unsigned int)((temp2 + (temp1 - temp2)*6*temp3)*100);
	else
		if((temp3 * 2) < 1)
			c = (unsigned int)(temp1*100);
		else
			if((temp3 * 3) < 2)
				c = (unsigned int)((temp2 + (temp1 - temp2)*(.66666 - temp3)*6)*100);
			else
				c = (unsigned int)(temp2*100);
	return;
}

// This function converts the "color" object to the equivalent RGB values of
// the hue, saturation, and luminance passed as h, s, and l respectively
string HSLtoRGB(const unsigned int& h, const unsigned int& s, const unsigned int& l)
{
	unsigned int r = 0;
	unsigned int g = 0;
	unsigned int b = 0;

	float L = ((float)l)/100;
	float S = ((float)s)/100;
	float H = ((float)h)/360;

	if(s == 0)
	{
		r = l;
		g = l;
		b = l;
	}
	else
	{
		float temp1 = 0;
		if(L < .50)
		{
			temp1 = L*(1 + S);
		}
		else
		{
			temp1 = L + S - (L*S);
		}

		float temp2 = 2*L - temp1;

		float temp3 = 0;
		for(int i = 0 ; i < 3 ; i++)
		{
			switch(i)
			{
			case 0: // red
				{
					temp3 = H + .33333f;
					if(temp3 > 1)
						temp3 -= 1;
					HSLtoRGB_Subfunction(r,temp1,temp2,temp3);
					break;
				}
			case 1: // green
				{
					temp3 = H;
					HSLtoRGB_Subfunction(g,temp1,temp2,temp3);
					break;
				}
			case 2: // blue
				{
					temp3 = H - .33333f;
					if(temp3 < 0)
						temp3 += 1;
					HSLtoRGB_Subfunction(b,temp1,temp2,temp3);
					break;
				}
			default:
				{

				}
			}
		}
	}
	r = (unsigned int)((((float)r)/100)*255);
	g = (unsigned int)((((float)g)/100)*255);
	b = (unsigned int)((((float)b)/100)*255);
	ostringstream ostr;
	ostr << "[" << r << "," << g << "," << b << "]";

	return ostr.str();
}


vector<string> GetRGBColors(int n)
{
	vector<string> colors;
	// assumes hue [0, 360), saturation [0, 100), lightness [0, 100)


	for(int i = 0; i < 360; i += 360 / n) {
		int hue = i;
		int saturation = 90 + (static_cast <float> (rand()) / static_cast <float> (RAND_MAX))*10;
		int lightness = 50 +  (static_cast <float> (rand()) / static_cast <float> (RAND_MAX))*10;
		
		colors.push_back(HSLtoRGB(hue, saturation, lightness));
	}

    return colors;
}

