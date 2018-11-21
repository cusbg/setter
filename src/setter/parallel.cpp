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
#include "parameters.h"
#include "parallel.h"

#if TBB

ParallelInput::ParallelInput(
	vector<cRNAStructure *> &rnaVect,
	MatchData ****data) :
	mRnaVect(rnaVect),
	mData(data)
{
	// no-op;
}

void ParallelInput::operator()(const tbb::blocked_range2d<size_t> &r) const
{
	for (size_t i = r.rows().begin(); i != r.rows().end(); ++i){
		for (size_t j = r.cols().begin(); j != r.cols().end(); ++j) {
			(*mData)[i][j] = new MatchData(mRnaVect, i, j);
		}
	}
}

ParallelAddRNA::ParallelAddRNA(
	unsigned int size,
	cRNAStructure &rna,
	string &rnaIdx,
	vector<cRNAStructure *> &rnaVect,
	vector<string> &map,
	MatchData ****data,
	bool init) :
	mSize(size),
	mRna(rna),
	mRnaIdx(rnaIdx),
	mRnaVect(rnaVect),
	mMap(map),
	mData(data),
	mInit(init)
{
	// no-op
}

void ParallelAddRNA::operator()(const tbb::blocked_range<size_t> &r) const
{
	for (size_t i = r.begin(); i != r.end(); ++i) {
		(*mData)[i][mSize] = new MatchData(mRnaVect[i], &mRna, mMap[i], mRnaIdx, mInit);
		(*mData)[mSize][i] = new MatchData(&mRna, mRnaVect[i], mRnaIdx, mMap[i], mInit);
	}
}


ParallelHairpinMatch::ParallelHairpinMatch(
	cRNAStructure &rnaQ,
	cRNAStructure &rnaDB,
	cHairpin &hpQ,
	cHairpin &hpDB,
	int ixNeckQ,
	int ixNeckDB,
	sParams &params,
	double scoreLimit,
	sResMatch &outRes,
	tbb::mutex &resMutex) :
	mRnaQ(rnaQ),
	mRnaDB(rnaDB),
	mHpQ(hpQ),
	mHpDB(hpDB),
	mIxNeckQ(ixNeckQ),
	mIxNeckDB(ixNeckDB),
	mParams(params),
	mScoreLimit(scoreLimit),
	mOutRes(outRes),
	mResMutex(resMutex)
{
	// no_op
}

void ParallelHairpinMatch::operator()(const tbb::blocked_range<size_t> &r) const
{
	int jMax = mHpDB.head.residues.size();

	cRMSDStructExt sQ;
	cRMSDStructExt sDB;
	sQ.length = sDB.length = 3;
	sQ.coord = new cRMSD3DCoord[sQ.length];
	sDB.coord = new cRMSD3DCoord[sDB.length];

	const unsigned int alignLen = 3;
	cRMSDAlign align[2] = { cRMSDAlign(alignLen), cRMSDAlign(alignLen) };
	for (unsigned int i = 0; i < 3; i++)
	{
		align[0][i].x = align[0][i].y = i;
		align[1][i].x = i;
		align[1][i].y = 2 - i;
	}

	//sResMatch resOpt;
	//cout << r.begin() << ";" << r.end() << endl;

	tPdbAtom atom;
	for (size_t i = r.begin(); i != r.end(); ++i){
		for (int j = 0; j < jMax; ++j) {
			for (int kQ = mIxNeckQ; kQ < mIxNeckQ + mParams.neckMaxShift; kQ++) {
				if (kQ < 0) continue;
				if (kQ >(int)mHpQ.stem.pairs.size() - 1
					||
					mHpQ.stem.pairs[kQ].r1.ix_residue < 0
					||
					mHpQ.stem.pairs[kQ].r2.ix_residue < 0) {
					continue;
				}

				for (int kDB = mIxNeckDB; kDB < mIxNeckDB + mParams.neckMaxShift; kDB++) {
					//e.g. in case of 1zih and 1zif we need to "move the neck" to obtain correct alignment
					if (kDB < 0) continue;
					if (kDB >(int)mHpDB.stem.pairs.size() - 1
						|| mHpDB.stem.pairs[kDB].r1.ix_residue < 0
						|| mHpDB.stem.pairs[kDB].r2.ix_residue < 0) {
						continue;
					}

					tPairResiude prQ[3] = { mHpQ.stem.pairs[kQ].r1, mHpQ.head.residues[i].r1, mHpQ.stem.pairs[kQ].r2 };
					tPairResiude prDB[3] = { mHpDB.stem.pairs[kDB].r1, mHpDB.head.residues[j].r1, mHpDB.stem.pairs[kDB].r2 };

					atom = mRnaQ.GetAtom(prQ[0]); sQ.coord[0] = cRMSD3DCoord(atom.coords.x, atom.coords.y, atom.coords.z);
					atom = mRnaQ.GetAtom(prQ[1]); sQ.coord[1] = cRMSD3DCoord(atom.coords.x, atom.coords.y, atom.coords.z);
					atom = mRnaQ.GetAtom(prQ[2]); sQ.coord[2] = cRMSD3DCoord(atom.coords.x, atom.coords.y, atom.coords.z);

					atom = mRnaDB.GetAtom(prDB[0]); sDB.coord[0] = cRMSD3DCoord(atom.coords.x, atom.coords.y, atom.coords.z);
					atom = mRnaDB.GetAtom(prDB[1]); sDB.coord[1] = cRMSD3DCoord(atom.coords.x, atom.coords.y, atom.coords.z);
					atom = mRnaDB.GetAtom(prDB[2]); sDB.coord[2] = cRMSD3DCoord(atom.coords.x, atom.coords.y, atom.coords.z);

					//alignActual arrays contain positions of aligned atoms in the original PDB (RNA) structure
					//thus one can access aligned atoms' properties (such as coordinates)
					int alignActualQ[3][2] = {
						{ prQ[0].ix_chain, prQ[0].ix_residue },
						{ prQ[1].ix_chain, prQ[1].ix_residue },
						{ prQ[2].ix_chain, prQ[2].ix_residue }
					};
					int alignActualDB[3][2] = {
						{ prDB[0].ix_chain, prDB[0].ix_residue },
						{ prDB[1].ix_chain, prDB[1].ix_residue },
						{ prDB[2].ix_chain, prDB[2].ix_residue }
					};

					//according to one of the bioinformatics referees, it does not make sense to try align 5'->3' agains 3'->5'
					// in such case, l should be upper limited by 1, else by 2
					//for (int l = 0; l < 2; l++)
					for (int l = 0; l < 1; l++) {

						sResMatch res = MatchHairpinsBasedOnTriplets(sQ, mRnaQ, mHpQ, sDB, mRnaDB, mHpDB, mParams, align[l], alignActualQ, alignActualDB, min(mScoreLimit, mOutRes.score));

						mResMutex.lock();
						if (res.score < mOutRes.score) {
							//resOpt = res;
							mOutRes = res;
						}
						mResMutex.unlock();
					}
				}
			}
		}
	}
}

void ParallelHairpinMatch::operator()(const tbb::blocked_range2d<size_t> &r) const
{
	cRMSDStructExt sQ;
	cRMSDStructExt sDB;
	sQ.length = sDB.length = 3;
	sQ.coord = new cRMSD3DCoord[sQ.length];
	sDB.coord = new cRMSD3DCoord[sDB.length];

	const unsigned int alignLen = 3;
	cRMSDAlign align[2] = { cRMSDAlign(alignLen), cRMSDAlign(alignLen) };
	for (unsigned int i = 0; i < 3; i++)
	{
		align[0][i].x = align[0][i].y = i;
		align[1][i].x = i;
		align[1][i].y = 2 - i;
	}

	sResMatch resOpt;

	tPdbAtom atom;
	for (size_t i = r.rows().begin(); i != r.rows().end(); ++i){
		for (size_t j = r.cols().begin(); j != r.cols().end(); ++j) {
			for (int kQ = mIxNeckQ; kQ < mIxNeckQ + mParams.neckMaxShift; kQ++) {
				if (kQ < 0) continue;
				if (kQ >(int)mHpQ.stem.pairs.size() - 1
					||
					mHpQ.stem.pairs[kQ].r1.ix_residue < 0
					||
					mHpQ.stem.pairs[kQ].r2.ix_residue < 0) {
					continue;
				}

				for (int kDB = mIxNeckDB; kDB < mIxNeckDB + mParams.neckMaxShift; kDB++) {
					//e.g. in case of 1zih and 1zif we need to "move the neck" to obtain correct alignment
					if (kDB < 0) continue;
					if (kDB >(int)mHpDB.stem.pairs.size() - 1
						|| mHpDB.stem.pairs[kDB].r1.ix_residue < 0
						|| mHpDB.stem.pairs[kDB].r2.ix_residue < 0) {
						continue;
					}

					tPairResiude prQ[3] = { mHpQ.stem.pairs[kQ].r1, mHpQ.head.residues[i].r1, mHpQ.stem.pairs[kQ].r2 };
					tPairResiude prDB[3] = { mHpDB.stem.pairs[kDB].r1, mHpDB.head.residues[j].r1, mHpDB.stem.pairs[kDB].r2 };

					atom = mRnaQ.GetAtom(prQ[0]); sQ.coord[0] = cRMSD3DCoord(atom.coords.x, atom.coords.y, atom.coords.z);
					atom = mRnaQ.GetAtom(prQ[1]); sQ.coord[1] = cRMSD3DCoord(atom.coords.x, atom.coords.y, atom.coords.z);
					atom = mRnaQ.GetAtom(prQ[2]); sQ.coord[2] = cRMSD3DCoord(atom.coords.x, atom.coords.y, atom.coords.z);

					atom = mRnaDB.GetAtom(prDB[0]); sDB.coord[0] = cRMSD3DCoord(atom.coords.x, atom.coords.y, atom.coords.z);
					atom = mRnaDB.GetAtom(prDB[1]); sDB.coord[1] = cRMSD3DCoord(atom.coords.x, atom.coords.y, atom.coords.z);
					atom = mRnaDB.GetAtom(prDB[2]); sDB.coord[2] = cRMSD3DCoord(atom.coords.x, atom.coords.y, atom.coords.z);

					//alignActual arrays contain positions of aligned atoms in the original PDB (RNA) structure
					//thus one can access aligned atoms' properties (such as coordinates)
					int alignActualQ[3][2] = {
						{ prQ[0].ix_chain, prQ[0].ix_residue },
						{ prQ[1].ix_chain, prQ[1].ix_residue },
						{ prQ[2].ix_chain, prQ[2].ix_residue }
					};
					int alignActualDB[3][2] = {
						{ prDB[0].ix_chain, prDB[0].ix_residue },
						{ prDB[1].ix_chain, prDB[1].ix_residue },
						{ prDB[2].ix_chain, prDB[2].ix_residue }
					};

					//according to one of the Bioinformatics journal referees, it does not make sense to try align 5'->3' agains 3'->5'
					// in such case, l should be upper limited by 1, else by 2
					//for (int l = 0; l < 2; l++)
					for (int l = 0; l < 1; l++) {
						mResMutex.lock();
						resOpt = mOutRes;
						mResMutex.unlock();

						sResMatch res = MatchHairpinsBasedOnTriplets(sQ, mRnaQ, mHpQ, sDB, mRnaDB, mHpDB, mParams, align[l], alignActualQ, alignActualDB, min(mScoreLimit, resOpt.score));

						mResMutex.lock();
						if (res.score < mOutRes.score) {
							mOutRes = res;
						}
						mResMutex.unlock();
					}
				}
			}
		}
	}
}

void ParallelHairpinMatch::operator()(int iMax, int jMax) const
{
	cRMSDStructExt sQ;
	cRMSDStructExt sDB;
	sQ.length = sDB.length = 3;
	sQ.coord = new cRMSD3DCoord[sQ.length];
	sDB.coord = new cRMSD3DCoord[sDB.length];

	const unsigned int alignLen = 3;
	cRMSDAlign align[2] = { cRMSDAlign(alignLen), cRMSDAlign(alignLen) };
	for (unsigned int i = 0; i < 3; i++)
	{
		align[0][i].x = align[0][i].y = i;
		align[1][i].x = i;
		align[1][i].y = 2 - i;
	}

	sResMatch resOpt;

	tPdbAtom atom;
	for (int i = 0; i < iMax; ++i){
		for (int j = 0; j < jMax; ++j) {
			for (int kQ = mIxNeckQ; kQ < mIxNeckQ + mParams.neckMaxShift; kQ++) {
				if (kQ < 0) continue;
				if (kQ >(int)mHpQ.stem.pairs.size() - 1
					||
					mHpQ.stem.pairs[kQ].r1.ix_residue < 0
					||
					mHpQ.stem.pairs[kQ].r2.ix_residue < 0) {
					continue;
				}

				for (int kDB = mIxNeckDB; kDB < mIxNeckDB + mParams.neckMaxShift; kDB++) {
					//e.g. in case of 1zih and 1zif we need to "move the neck" to obtain correct alignment
					if (kDB < 0) continue;
					if (kDB >(int)mHpDB.stem.pairs.size() - 1
						|| mHpDB.stem.pairs[kDB].r1.ix_residue < 0
						|| mHpDB.stem.pairs[kDB].r2.ix_residue < 0) {
						continue;
					}

					tPairResiude prQ[3] = { mHpQ.stem.pairs[kQ].r1, mHpQ.head.residues[i].r1, mHpQ.stem.pairs[kQ].r2 };
					tPairResiude prDB[3] = { mHpDB.stem.pairs[kDB].r1, mHpDB.head.residues[j].r1, mHpDB.stem.pairs[kDB].r2 };

					atom = mRnaQ.GetAtom(prQ[0]); sQ.coord[0] = cRMSD3DCoord(atom.coords.x, atom.coords.y, atom.coords.z);
					atom = mRnaQ.GetAtom(prQ[1]); sQ.coord[1] = cRMSD3DCoord(atom.coords.x, atom.coords.y, atom.coords.z);
					atom = mRnaQ.GetAtom(prQ[2]); sQ.coord[2] = cRMSD3DCoord(atom.coords.x, atom.coords.y, atom.coords.z);

					atom = mRnaDB.GetAtom(prDB[0]); sDB.coord[0] = cRMSD3DCoord(atom.coords.x, atom.coords.y, atom.coords.z);
					atom = mRnaDB.GetAtom(prDB[1]); sDB.coord[1] = cRMSD3DCoord(atom.coords.x, atom.coords.y, atom.coords.z);
					atom = mRnaDB.GetAtom(prDB[2]); sDB.coord[2] = cRMSD3DCoord(atom.coords.x, atom.coords.y, atom.coords.z);

					//alignActual arrays contain positions of aligned atoms in the original PDB (RNA) structure
					//thus one can access aligned atoms' properties (such as coordinates)
					int alignActualQ[3][2] = {
						{ prQ[0].ix_chain, prQ[0].ix_residue },
						{ prQ[1].ix_chain, prQ[1].ix_residue },
						{ prQ[2].ix_chain, prQ[2].ix_residue }
					};
					int alignActualDB[3][2] = {
						{ prDB[0].ix_chain, prDB[0].ix_residue },
						{ prDB[1].ix_chain, prDB[1].ix_residue },
						{ prDB[2].ix_chain, prDB[2].ix_residue }
					};

					//according to one of the bioinformatics referees, it does not make sense to try align 5'->3' agains 3'->5'
					// in such case, l should be upper limited by 1, else by 2
					//for (int l = 0; l < 2; l++)
					for (int l = 0; l < 1; l++) {
						sResMatch res = MatchHairpinsBasedOnTriplets(sQ, mRnaQ, mHpQ, sDB, mRnaDB, mHpDB, mParams, align[l], alignActualQ, alignActualDB, min(mScoreLimit, resOpt.score));
						if (res.score < resOpt.score) {
							resOpt = res;
						}
					}
				}
			}
		}
	}
}


ParallelMainBlock::ParallelMainBlock(
	cRNAStructure &aver,
	vector<int> &positions,
	vector<cRNAStructure*> &dataset,
	vector<string> &names,
	vector<sResMatch> &results) :
	mAver(aver),
	mPositions(positions),
	mDataset(dataset),
	mNames(names),
	mResults(results),
	param(Parameters::GetInstance())
{
	// no-oper
}

void ParallelMainBlock::operator()(const tbb::blocked_range<size_t> &r) const
{
	for (size_t i = r.begin(); i != r.end(); ++i) {
		int pos = mPositions[i];
		cRNAStructure &rna = *mDataset[pos];

		sParams params = GetGlobalParams();
		sResMatch res = CalcScore(mAver, rna, params);

		mResults[i] = res;

	}
}

void ParallelMainBlock::PrintOutput()
{
	for (int i = 0; i < mResults.size(); ++i) {
		cout << mNames[mPositions[i]] << " - aver : " << mResults[i].score << "; " << endl;
	}
}

#endif