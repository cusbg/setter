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
#include <math.h>
#include "make_aver_rna.h"
#include "rotation.h"
#include "parameters.h"
#include "algo.h"

#include <time.h>

FormAvg::FormAvg(MatchData &data) :
mData(data)
{
	Calc();
}

FormAvg::FormAvg(double mean, MatchData &data, int &cnt, double &result) :
mMean(mean),
mData(data)
{
	if (mData.GetResult().score < 0.0000000001) {
		cout << "no need to compute" << endl;
		mResultRna = *mData.mRnaQ;
	}
	else {
		Calc();
	}

	vector<double> dist(4);

	//int start = clock();
	ParallelNewRnaDis parallelNewRnaDis(mResultRna, *mData.mRnaQ, *mData.mRnaDB, dist);

	tbb::parallel_for(
		tbb::blocked_range<size_t>(0, 4), parallelNewRnaDis);

	//mResultRna.RenameWithSeqOrder();

	//int end = clock();

	//cout << "distance calck: " << end - start << endl;

	double match1 = CalcFinalScore(dist[0], dist[1]);
	double match2 = CalcFinalScore(dist[2], dist[3]);

	if (match1 == sResMatch::MAX || match2 == sResMatch::MAX) {
		mResultRna = *mData.mRnaQ;
	}

	++cnt;
	result += data.GetResult().score;
}

void FormAvg::Calc()
{
	int start, end;

	GSSUPair pair;

	// contains the original hairpins as RNAs
	RNAVect rnaQVect(mData.mOptPairs.size());
	RNAVect rnaDBVect(mData.mOptPairs.size());
	RNAVect fullQVect(mData.mCntHpQ);

	mRnaQVect.resize(mData.mOptPairs.size());
	mRnaDBVect.resize(mData.mOptPairs.size());

	RNAVect mergedRnaVect(mData.mOptPairs.size());
	vector<string> output(mData.mOptPairs.size());
	ParallelFullQVect parallelFullQVect(mData, fullQVect);
	ParallelMerge parallelMerge(*this, rnaQVect, rnaDBVect, mergedRnaVect, output);

#if PARALLEL_INPUT
	start = clock();

	tbb::parallel_for(
		tbb::blocked_range<size_t>(0, mData.mCntHpQ, 10), parallelFullQVect);

	tbb::parallel_for(
		tbb::blocked_range<size_t>(0, mData.mOptPairs.size(), 4), parallelMerge);

	end = clock();
	//cout << "parallel section ends " << end - start << endl;
#else
	parallelFullQVect.serial(mData.mCntHpQ);
	parallelMerge.serial(mData.mOptPairs.size());
#endif
	int pos = 0;
	int residueIdx;
	if (mData.mOptPairs.size() > pos) {
		residueIdx = mData.mOptPairs[pos].q;
	}
	else {
		residueIdx = -1;
	}

	Parameters &param = Parameters::GetInstance();

	start = clock();

	for (unsigned i = 0; i < mData.mCntHpQ; ++i) {
		cRNAStructure rna;
		if (i == residueIdx) {
			rna = mergedRnaVect[pos];
			++pos;
			if (mData.mOptPairs.size() > pos) {
				residueIdx = mData.mOptPairs[pos].q;
			}
			else {
				residueIdx = -1;
			}
		}
		else {
			rna = fullQVect[i];
		}

		if (i == 0) {
			mResultRna = rna;
		}
		else {
			AppendRNA(mResultRna, rna);
		}
	}

	end = clock();
}

tPairResiude FormAvg::GetFirstResidue(cHairpin &hp, int pos)
{
	int count = 0;

	tPairResiude first;
	// search in stem.r2 side
	for (int idx = hp.stem.pairs.size() - 1; idx >= 0; --idx) {
		first = hp.stem.pairs[idx].r1;

		if (first.ix_residue != -1) {
			if (count == pos) {
				return first;
			}
			else {
				++count;
			}
		}
	}

	// search in head
	for (int idx = 0; idx < hp.head.residues.size(); ++idx) {
		first = hp.head.residues[idx].r1;

		if (first.ix_residue != -1) {
			if (count == pos) {
				return first;
			}
			else {
				++count;
			}
		}
	}

	// search in stem.r1 side
	for (int idx = 0; idx < hp.stem.pairs.size(); ++idx) {
		first = hp.stem.pairs[idx].r2;

		if (first.ix_residue != -1) {
			if (count == pos) {
				return first;
			}
			else {
				++count;
			}
		}
	}

	// no result;
	return tPairResiude();
}


tPairResiude FormAvg::GetLastResidue(cHairpin &hp, int pos)
{
	int count = 0;

	tPairResiude last;
	// search in stem.r2 side
	for (int idx = hp.stem.pairs.size() - 1; idx >= 0; --idx) {
		last = hp.stem.pairs[idx].r2;

		if (last.ix_residue != -1) {
			if (count == pos) {
				return last;
			}
			else {
				++count;
			}
		}
	}

	// search in head
	for (int idx = hp.head.residues.size() - 1; idx >= 0; --idx) {
		last = hp.head.residues[idx].r1;

		if (last.ix_residue != -1) {
			if (count == pos) {
				return last;
			}
			else {
				++count;
			}
		}
	}

	// search in stem.r1 side
	for (int idx = 0; idx < hp.stem.pairs.size(); ++idx) {
		last = hp.stem.pairs[idx].r1;

		if (last.ix_residue != -1) {
			if (count == pos) {
				return last;
			}
			else {
				++count;
			}
		}
	}

	// no result;
	return tPairResiude();
}


cRNAStructure FormAvg::getResult()
{
	return mResultRna;
}

void FormAvg::MakeHairpin(
	const cHairpin &hp,
	cHairpin &outHp)
{
	vector<tPair>::const_iterator iter;
	tPair residuePair;
	for (iter = hp.stem.pairs.begin(); iter != hp.stem.pairs.end(); ++iter) {
		residuePair = *iter;
		//if (residuePair.r1.ix_residue >= 0 && residuePair.r2.ix_residue >= 0) {
		// paired residues
		outHp.stem.pairs.push_back(residuePair);
		//}
	}

	for (iter = hp.head.residues.begin(); iter != hp.head.residues.end(); iter++) {
		residuePair = *iter;
		outHp.head.residues.push_back(residuePair);
	}
}

void FormAvg::MakeRNAFromHairpin(
	cRNAStructure &rna,
	cHairpin &hp,
	hairpinIdx idx,
	cRNAStructure &outRna)
{
	stringstream sstream;
	int ixS = 0;
	string sxS;

	vector<string> chain;
	tPdbModel &m = rna.mModels[0];

	tPairResiude	*tp;
	tPdbModel		newModel;
	tPdbModel		newModelNt;
	tPdbChain		newCh;
	tPdbChain		newChNt;

	tPair p;
	if (hp.stem.pairs.size() > 0) {
		p = hp.stem.pairs.back();
	}
	else if (hp.head.residues.size() > 0) {
		p.r1 = hp.head.residues.front().r1;
		p.r2 = hp.head.residues.back().r1;
	}
	else {
		return;
	}

	string residue;
	tPdbAtom atom;
	tPdbAtom atomIt;
	unsigned j = 0;
	for (int mode = 0; mode < 3; ++mode){
		unsigned size;
		switch (mode) {
		case 0: size = hp.stem.pairs.size();	break;
		case 1: size = hp.head.residues.size(); break;
		case 2: size = hp.stem.pairs.size();	break;
		}
		for (unsigned i = 0; i < size; ++i) {
			switch (mode) {
			case 0: tp = &hp.stem.pairs[size - 1 - i].r1;	break;
			case 1: tp = &hp.head.residues[i].r1;			break;
			case 2: tp = &hp.stem.pairs[i].r2;				break;
			}

			if (tp->ix_residue == -1) {
				tp->ix_chain = 0;
				continue;
			}
			sstream.str(string());
			sstream.clear();
			sstream << ixS + 1;
			sxS = sstream.str();

			try {
				residue = "   ";
				residue.append(tp->residue_label_short);
			}
			catch (std::exception &e) {
				// no-op jet
			}
			chain.push_back(residue);

			atom = rna.GetAtom(*tp);
			tPdbChain &ch = m[tp->ix_chain];

			bool found = false;
			while (!found){
				for (; j < ch.size(); ++j) {
					if (ch[j].residue_num == atom.residue_num) {
						found = true;
						break;
					}
				}
				if (!found) j = 0;
			}
			for (; j < ch.size(); ++j) {
				if (ch[j].residue_num == atom.residue_num) {
					atomIt = ch[j];
					atomIt.residue_num = sxS;
					newCh.push_back(atomIt);
					//cout << atomIt.residue_num << ";";
				}
				else {
					//break;
				}
			}
			atom.residue_num = sxS;
			newChNt.push_back(atom);
			sstream.str(string());
			sstream.clear();
			sstream << ixS;
			tp->residue_position = sstream.str();
			tp->ix_residue = ixS++;
			tp->ix_chain = 0;
		}
	}

	outRna.mPrimarySeqChains.insert(outRna.mPrimarySeqChains.end(), chain);

	newModel.push_back(newCh);
	newModelNt.push_back(newChNt);

	outRna.mModels.push_back(newModel);
	outRna.mModelsNt.push_back(newModelNt);

	string pdbID = rna.mPdbID;
	sstream.str(string());
	sstream.clear();
	sstream << idx;
	outRna.mPdbID = pdbID.append(sstream.str());
	outRna.mHairpins.push_back(hp);
}

void FormAvg::Rotate(
	cRNAStructure &rna,
	hairpinIdx hpQIdx,
	hairpinIdx hpDBIdx,
	InputType inputType)
{
	Parameters &param = Parameters::GetInstance();
	cRMSDMatrix data = mData(hpQIdx, hpDBIdx).rot;
	RMatrix matrix(data);
	RMatrix halfRotation = matrix.GetHalfRotation();


	switch (inputType) {
	case InputType::Q:
		if (param.GetHalfRotation()) {
			halfRotation = halfRotation.Transpose();
			Rotate(rna, halfRotation);
		}
		break;
	case InputType::DB:
	{
		if (param.GetHalfRotation()) {
			Rotate(rna, halfRotation);
		}
		else {
			Rotate(rna, matrix);
		}
		cRMSD3DCoord coords = mData(hpQIdx, hpDBIdx).trans;
		Translate(rna, coords);
		break;
	}
	case InputType::NO:
		// do nothing
		break;
	}
}

void FormAvg::Translate(cRNAStructure &rna, cRMSD3DCoord &tVector)
{
	if (rna.mModels.size() == 0) {
		return;
	}

	cRMSD3DCoord coord, coordTrans;

	tPdbChain &ch = rna.mModels[0][0];
	for (int i = 0; i < ch.size(); ++i) {
		coord[0] = ch[i].coords.x; coord[1] = ch[i].coords.y;	coord[2] = ch[i].coords.z;
		coordTrans = coord + tVector;
		ch[i].coords.x = coordTrans[0]; ch[i].coords.y = coordTrans[1]; ch[i].coords.z = coordTrans[2];
	}

	tPdbChain &chNt = rna.mModelsNt[0][0];
	for (int i = 0; i < chNt.size(); ++i) {
		coord[0] = chNt[i].coords.x; coord[1] = chNt[i].coords.y;	coord[2] = chNt[i].coords.z;
		coordTrans = coord + tVector;
		chNt[i].coords.x = coordTrans[0]; chNt[i].coords.y = coordTrans[1]; chNt[i].coords.z = coordTrans[2];
	}
}

void FormAvg::Rotate(cRNAStructure &rna, RMatrix &rMatrix)
{
	if (rna.mModels.size() == 0) {
		return;
	}

	cRMSD3DCoord coord, coordRot;

	tPdbChain &ch = rna.mModels[0][0];
	for (int i = 0; i < ch.size(); ++i) {
		coord[0] = ch[i].coords.x; coord[1] = ch[i].coords.y;	coord[2] = ch[i].coords.z;
		coordRot = rMatrix * coord;
		ch[i].coords.x = coordRot[0]; ch[i].coords.y = coordRot[1]; ch[i].coords.z = coordRot[2];
	}

	tPdbChain &chNt = rna.mModelsNt[0][0];
	for (int i = 0; i < chNt.size(); ++i) {
		coord[0] = chNt[i].coords.x; coord[1] = chNt[i].coords.y;	coord[2] = chNt[i].coords.z;
		coordRot = rMatrix * coord;
		chNt[i].coords.x = coordRot[0]; chNt[i].coords.y = coordRot[1]; chNt[i].coords.z = coordRot[2];
	}
}
