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
#include "parameters.h"

#if TBB

#include "tbb/parallel_for.h"
#include "tbb/parallel_reduce.h"
#include "tbb/blocked_range.h"
#include "tbb/blocked_range2d.h"
#include "tbb/mutex.h"

#include "algo.h"
#include "input.h"
#include "make_aver_rna.h"

//! Intel TBB class to read and parsethe input RNA structures.
class ParallelInput
{
public:
	ParallelInput(
		vector<cRNAStructure *> &rnaVect,
		MatchData ****data);

	void operator()(const tbb::blocked_range2d<size_t> &r) const;

private:
	vector<cRNAStructure *> &mRnaVect;
	MatchData ****mData;
};

//! Intel TBB class to add new RNA to input RNAs, computes the distance matrix.
class ParallelAddRNA
{
public:
	ParallelAddRNA(
		unsigned int size,
		cRNAStructure &rna,
		string &rnaIdx,
		vector<cRNAStructure *> &rnaVect,
		vector<string> &mMap,
		MatchData ****data,
		bool init);

	void operator()(const tbb::blocked_range<size_t> &r) const;

private:
	unsigned int mSize;
	cRNAStructure &mRna;
	string &mRnaIdx;
	vector<cRNAStructure *> &mRnaVect;
	vector<string> &mMap;
	MatchData ****mData;
	bool mInit;
};

//! Intel TBB class to hairpin match, it is not recommented to use.
class ParallelHairpinMatch
{
public:
	ParallelHairpinMatch(
		cRNAStructure &rnaQ, 
		cRNAStructure &rnaDB,
		cHairpin &hpQ,
		cHairpin &hpDB,
		int ixNeckQ,
		int ixNeckDB,
		sParams &params,
		double scoreLimit,
		sResMatch &outRes,
		tbb::mutex &mResMutex);

	void operator()(const tbb::blocked_range<size_t> &r) const;
	void operator()(const tbb::blocked_range2d<size_t> &r) const;
	void operator()(int iMax, int jMax) const;

private:
	cRNAStructure &mRnaQ; 
	cRNAStructure &mRnaDB;
	cHairpin &mHpQ;
	cHairpin &mHpDB;
	sParams &mParams;
	int mIxNeckQ;
	int mIxNeckDB;
	double mScoreLimit;
	sResMatch &mOutRes;
	tbb::mutex &mResMutex;
};

//! Intel TBB class to match RNA structures for testing dataset with average RNA structures.
class ParallelMainBlock
{
public:
	ParallelMainBlock(
		cRNAStructure &aver,
		vector<int> &positions,
		vector<cRNAStructure*> &dataset,
		vector<string> &names,
		vector<sResMatch> &results);

	void operator()(const tbb::blocked_range<size_t> &r) const;
	void PrintOutput();
	inline vector<sResMatch> const &results() const { return mResults; }


private:
	cRNAStructure &mAver;
	vector<cRNAStructure*> &mDataset;
	vector<string> &mNames;
	vector<int> &mPositions;
	vector<sResMatch> &mResults;
	Parameters &param;
};

#endif