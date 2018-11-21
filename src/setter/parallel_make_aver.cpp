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
#include "make_aver_rna.h"

ParallelFullQVect::ParallelFullQVect(
	MatchData &data,
	RNAVect &fullQVect) :
	mData(data),
	mFullQVect(fullQVect)
{
	// no op
}

void ParallelFullQVect::operator()(const tbb::blocked_range<size_t> &r) const
{
	GSSUPair pair;
	for (size_t i = r.begin(); i != r.end(); ++i) {
		cRNAStructure &rnaQ = mFullQVect[i];
		cHairpin hpQ;
		FormAvg::MakeHairpin(mData.mRnaQ->mHairpins[i], hpQ);
		FormAvg::MakeRNAFromHairpin(*mData.mRnaQ, hpQ, i, rnaQ);
	}
}

void ParallelFullQVect::serial(int size) const
{
	GSSUPair pair;
	for (size_t i = 0; i < size; ++i) {
		cRNAStructure &rnaQ = mFullQVect[i];
		cHairpin hpQ;
		FormAvg::MakeHairpin(mData.mRnaQ->mHairpins[i], hpQ);
		FormAvg::MakeRNAFromHairpin(*mData.mRnaQ, hpQ, i, rnaQ);
	}
}

ParallelMerge::ParallelMerge(
	FormAvg &aver,
	RNAVect &rnaQVect,
	RNAVect &rnaDBVect,
	RNAVect &mergedRnaVect,
	vector<string> &output) :
	mAver(aver),
	mRnaQVect(rnaQVect),
	mRnaDBVect(rnaDBVect),
	mMergedRnaVect(mergedRnaVect),
	mOutput(output),
	param(Parameters::GetInstance())
{
	// no-op
}

void ParallelMerge::operator()(const tbb::blocked_range<size_t> &r) const
{
	GSSUPair pair;
	for (size_t i = r.begin(); i != r.end(); ++i) {
		pair = mAver.mData.mOptPairs[i];
		// create backup
		{
			cRNAStructure &rnaQ = mRnaQVect[i];
			cHairpin hpQ;
			FormAvg::MakeHairpin(mAver.mData.mRnaQ->mHairpins[pair.q], hpQ);
			FormAvg::MakeRNAFromHairpin(*mAver.mData.mRnaQ, hpQ, pair.q, rnaQ);

			cRNAStructure &rnaDB = mRnaDBVect[i];
			cHairpin hpDB;
			FormAvg::MakeHairpin(mAver.mData.mRnaDB->mHairpins[pair.db], hpDB);
			FormAvg::MakeRNAFromHairpin(*mAver.mData.mRnaDB, hpDB, pair.db, rnaDB);
		}
		// rotates newRnaQ
		{
			cRNAStructure &newRnaQ = mAver.mRnaQVect[i];
			cHairpin newHpQ;
			if (mAver.mData(pair.q, pair.db).score == sResMatch::MAX) {
				newRnaQ = mRnaQVect[i];
			}
			else {
				FormAvg::MakeHairpin(mAver.mData.mRnaQ->mHairpins[pair.q], newHpQ);
				FormAvg::MakeRNAFromHairpin(*mAver.mData.mRnaQ, newHpQ, pair.q, newRnaQ);
				mAver.Rotate(newRnaQ, pair.q, pair.db, InputType::Q);
			}
		}
		// rotates newRnaDB
		{
			cRNAStructure &newRnaDB = mAver.mRnaDBVect[i];
			cHairpin newHpDB;
			if (mAver.mData(pair.q, pair.db).score == sResMatch::MAX) {
				newRnaDB = mRnaQVect[i];
			}
			else {
				FormAvg::MakeHairpin(mAver.mData.mRnaDB->mHairpins[pair.db], newHpDB);
				FormAvg::MakeRNAFromHairpin(*mAver.mData.mRnaDB, newHpDB, pair.db, newRnaDB);
				mAver.Rotate(newRnaDB, pair.q, pair.db, InputType::DB);
			}
		}
		// merges pairs
		{
			cRNAStructure &rna = mMergedRnaVect[i];

			stringstream ss;
			double distance = mAver.mData(pair.q, pair.db).score;
			if (distance > mAver.mMean * param.GetDiverMean()) {
				rna = mRnaQVect[i];
			}
			else {
				mAver.MergePairToRNA(i, rna);
				sParams params = GetGlobalParams();
				double match11 = Match(rna, mRnaQVect[i], params).score;
				double match12 = Match(mRnaQVect[i], rna, params).score;

				double match21 = Match(rna, mRnaQVect[i], params).score;
				double match22 = Match(mRnaQVect[i], rna, params).score;

				double match1 = CalcFinalScore(match11, match12);
				double match2 = CalcFinalScore(match21, match22);

				if (match1 + match2 > ResMatch::MAX - 10) {
					rna = mRnaQVect[i];
				}
				else if (match1 + match2 > distance * param.GetDiverParent()) {
					rna = mRnaQVect[i];
				}
			}

			mOutput[i] = ss.str();
		}
	}
}

void ParallelMerge::serial(int size) const
{
	GSSUPair pair;
	for (size_t i = 0; i < size; ++i) {
		pair = mAver.mData.mOptPairs[i];
		// create backup
		{
			cRNAStructure &rnaQ = mRnaQVect[i];
			cHairpin hpQ;
			FormAvg::MakeHairpin(mAver.mData.mRnaQ->mHairpins[pair.q], hpQ);
			FormAvg::MakeRNAFromHairpin(*mAver.mData.mRnaQ, hpQ, pair.q, rnaQ);

			cRNAStructure &rnaDB = mRnaDBVect[i];
			cHairpin hpDB;
			FormAvg::MakeHairpin(mAver.mData.mRnaDB->mHairpins[pair.db], hpDB);
			FormAvg::MakeRNAFromHairpin(*mAver.mData.mRnaDB, hpDB, pair.db, rnaDB);
		}
		// rotates newRnaQ
		{
			cRNAStructure &newRnaQ = mAver.mRnaQVect[i];
			cHairpin newHpQ;
			if (mAver.mData(pair.q, pair.db).score == sResMatch::MAX) {
				newRnaQ = mRnaQVect[i];
			}
			else {
				FormAvg::MakeHairpin(mAver.mData.mRnaQ->mHairpins[pair.q], newHpQ);
				FormAvg::MakeRNAFromHairpin(*mAver.mData.mRnaQ, newHpQ, pair.q, newRnaQ);
				mAver.Rotate(newRnaQ, pair.q, pair.db, InputType::Q);
			}
		}
		// rotates newRnaDB
		{
			cRNAStructure &newRnaDB = mAver.mRnaDBVect[i];
			cHairpin newHpDB;
			if (mAver.mData(pair.q, pair.db).score == sResMatch::MAX) {
				newRnaDB = mRnaQVect[i];
			}
			else {
				FormAvg::MakeHairpin(mAver.mData.mRnaDB->mHairpins[pair.db], newHpDB);
				FormAvg::MakeRNAFromHairpin(*mAver.mData.mRnaDB, newHpDB, pair.db, newRnaDB);
				mAver.Rotate(newRnaDB, pair.q, pair.db, InputType::DB);
			}
		}
		// merges pairs
		{
			cRNAStructure &rna = mMergedRnaVect[i];

			double distance = mAver.mData(pair.q, pair.db).score;
			if (distance > mAver.mMean * param.GetDiverMean()) {
				rna = mRnaQVect[i];
			}
			else {
				mAver.MergePairToRNA(i, rna);
				double match11 = Match(rna, mRnaQVect[i], GetGlobalParams()).score;
				double match12 = Match(mRnaQVect[i], rna, GetGlobalParams()).score;

				double match21 = Match(rna, mRnaQVect[i], GetGlobalParams()).score;
				double match22 = Match(mRnaQVect[i], rna, GetGlobalParams()).score;

				double match1 = CalcFinalScore(match11, match12);
				double match2 = CalcFinalScore(match21, match22);

				if (match1 + match2 > distance * param.GetDiverParent()) {
					rna = mRnaQVect[i];
				}
			}
		}
	}
}

void ParallelMerge::Print()
{
	for (int i = 0; i < mOutput.size(); ++i) {
		cout << mOutput[i] << endl;
	}
}

ParallelNewRnaDis::ParallelNewRnaDis(
	cRNAStructure &newRna,
	cRNAStructure &rnaQ,
	cRNAStructure &rnaDB,
	vector<double> &dist) :
	mNewRna(newRna),
	mRnaQ(rnaQ),
	mRnaDB(rnaDB),
	mDist(dist)
{
	// no-op
}

void ParallelNewRnaDis::operator()(const tbb::blocked_range<size_t> &r) const
{
	for (size_t i = r.begin(); i != r.end(); ++i) {
		switch (i) {
		case 0: mDist[0] = Match(mNewRna, mRnaQ, GetGlobalParams()).score; break;
		case 1: mDist[1] = Match(mRnaQ, mNewRna, GetGlobalParams()).score; break;
		case 2: mDist[2] = Match(mNewRna, mRnaDB, GetGlobalParams()).score; break;
		case 3: mDist[3] = Match(mRnaDB, mNewRna, GetGlobalParams()).score; break;
		}
	}
}

void ParallelNewRnaDis::serial(int size) const
{
	for (size_t i = 0; i < size; ++i) {
		switch (i) {
		case 0: mDist[0] = Match(mNewRna, mRnaQ, GetGlobalParams()).score; break;
		case 1: mDist[1] = Match(mRnaQ, mNewRna, GetGlobalParams()).score; break;
		case 2: mDist[2] = Match(mNewRna, mRnaDB, GetGlobalParams()).score; break;
		case 3: mDist[3] = Match(mRnaDB, mNewRna, GetGlobalParams()).score; break;
		}
	}
}