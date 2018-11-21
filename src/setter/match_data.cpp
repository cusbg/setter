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
#include <sstream>
#include "match_data.h"
#include "algo.h"

sParams MatchData::mDefParams = GetGlobalParams();

MatchData::MatchData() :
mRnaQ(NULL),
mRnaDB(NULL),
mCntHpQ(0),
mCntHpDB(0)
{
}

MatchData::MatchData(
	cRNAStructure *rnaQ,
	cRNAStructure *rnaDB,
	bool init) :
	mRnaQ(rnaQ),
	mRnaDB(rnaDB),
	mInit(init)
{
	mCntHpQ = rnaQ->GetHairpinsCount();
	mCntHpDB = rnaDB->GetHairpinsCount();

	stringstream sstream;
	sstream << 0;
	sRnaQIdx = sstream.str();

	sstream.str(string());
	sstream.clear();
	sstream << 1;
	sRnaDBIdx = sstream.str();

	if (mInit) {
		Initalize();
	}
}

MatchData::MatchData(
	vector<cRNAStructure *> &rnaVect,
	int rnaQIdx,
	int rnaDBIdx,
	bool init) :
	mRnaQ(rnaVect[rnaQIdx]),
	mRnaDB(rnaVect[rnaDBIdx]),
	mInit(init)
{
	mCntHpQ = mRnaQ->GetHairpinsCount();
	mCntHpDB = mRnaDB->GetHairpinsCount();

	stringstream sstream;
	sstream << "" << rnaQIdx;
	sRnaQIdx = sstream.str();

	sstream.str(string());
	sstream.clear();
	sstream << rnaDBIdx;
	sRnaDBIdx = sstream.str();

	//cout << rnaQIdx << ", " << rnaDBIdx << ": " << endl;
	if (mInit) {
		Initalize();
	}
}

MatchData::MatchData(
	cRNAStructure *rnaQ,
	cRNAStructure *rnaDB,
	string rnaQIdx,
	string rnaDBIdx,
	bool init) :
	mRnaQ(rnaQ),
	mRnaDB(rnaDB),
	sRnaQIdx(rnaQIdx),
	sRnaDBIdx(rnaDBIdx),
	mInit(init)
{
	mCntHpQ = rnaQ->GetHairpinsCount();
	mCntHpDB = rnaDB->GetHairpinsCount();

	//cout << rnaQIdx << ", " << rnaDBIdx << ": " << endl;
	if (init) {
		Initalize();
	}
}

MatchData::MatchData(MatchData &data)
{
	mInit = data.mInit;
	mOptPairs = data.mOptPairs;
	mResult = data.mResult;
	mRnaQ = data.mRnaQ;
	mRnaDB = data.mRnaDB;
	mCntHpQ = data.mCntHpQ;
	mCntHpDB = data.mCntHpDB;
	sRnaQIdx = data.sRnaQIdx;
	sRnaDBIdx = data.sRnaDBIdx;

	if (mInit) {
		mRessMap = data.mRessMap;
	}
}

void MatchData::operator=(MatchData &data)
{
	mInit = data.mInit;
	mOptPairs = data.mOptPairs;
	mResult = data.mResult;
	mRnaQ = data.mRnaQ;
	mRnaDB = data.mRnaDB;
	mCntHpQ = data.mCntHpQ;
	mCntHpDB = data.mCntHpDB;
	sRnaQIdx = data.sRnaQIdx;
	sRnaDBIdx = data.sRnaDBIdx;

	if (mInit) {
		mRessMap = data.mRessMap;
	}
}

sResMatch MatchData::operator()(int i, int j)
{
	if (mRessMap.find(make_pair(i, j)) != mRessMap.end()) {
		return mRessMap[make_pair(i, j)];
	}
	else {
		return sResMatch();
	}
}

MatchData::~MatchData()
{
	// no-op
}

void MatchData::Initalize()
{

	vector<GSSUPair> optP;

	sResMatch **ress = new sResMatch *[mCntHpQ];
	for (int i = 0; i < mCntHpQ; ++i) {
		ress[i] = new sResMatch[mCntHpDB];
	}

	mResult = Match(*mRnaQ, *mRnaDB, GetGlobalParams(), &ress, optP);
	sort(optP.begin(), optP.end(), GSSUPair::compareGSSUPair);

	int last = -1;
	for (int i = 0; i < optP.size(); ++i) {
		if (optP[i].q != last) {
			mOptPairs.push_back(optP[i]);
			last = optP[i].q;
		}
	}

	pair<int, int> actual;
	for (int i = 0; i < mOptPairs.size(); ++i) {
		actual = make_pair(mOptPairs[i].q, mOptPairs[i].db);
		mRessMap[actual] = ress[mOptPairs[i].q][mOptPairs[i].db];

		actual = make_pair(mOptPairs[i].db, mOptPairs[i].q);
		mRessMap[actual] = ress[mOptPairs[i].q][mOptPairs[i].db];
	}

	if (mCntHpQ == 0) {
		return;
	}
	for (int i = 0; i < mCntHpQ; ++i) {
		delete[] ress[i];
	}
	delete[] ress;

	mInit = true;
}

sResMatch &MatchData::GetResult()
{
	if (!mInit) {
		Initalize();
	}
	return mResult;
}