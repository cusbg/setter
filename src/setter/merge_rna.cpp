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
#include "math.h"
#include "algo.h"

#define STEM_R1 0
#define HEAD	1
#define STEM_R2 2

int roundDouble(double number)
{
	return number < 0.0 ? ceil(number - 0.5) : floor(number + 0.5);
}

void FormAvg::AppendRNA(
	cRNAStructure &appendedRna,
	cRNAStructure &rna)
{
	stringstream sstream;

	appendedRna.mPrimarySeqChains[0].insert(appendedRna.mPrimarySeqChains[0].end(), rna.mPrimarySeqChains[0].begin(), rna.mPrimarySeqChains[0].end());

	tPdbChain &aCh = appendedRna.mModels[0][0];
	tPdbChain &ch = rna.mModels[0][0];
	int res_num = atoi(aCh.back().residue_num.c_str());
	for (int i = 0; i < ch.size(); ++i) {
		aCh.push_back(ch[i]);
		sstream.str(string());
		sstream.clear();
		sstream << res_num + atoi(ch[i].residue_num.c_str());
		aCh.back().residue_num = sstream.str();
	}

	tPdbChain &aChNt = appendedRna.mModelsNt[0][0];
	tPdbChain &chNt = rna.mModelsNt[0][0];
	res_num = atoi(aChNt.back().residue_num.c_str());
	for (int i = 0; i < chNt.size(); ++i) {
		aChNt.push_back(chNt[i]);
		sstream.str(string());
		sstream.clear();
		sstream << res_num + atoi(chNt[i].residue_num.c_str());
		aChNt.back().residue_num = sstream.str();
	}

	tPairResiude residue = GetLastResidue(appendedRna.mHairpins.back(), 0);
	int ix_residue = residue.ix_residue + 1;
	int res_pos = atoi(residue.residue_position.c_str());

	appendedRna.mHairpins.push_back(rna.mHairpins[0]);
	cHairpin &hp = appendedRna.mHairpins.back();

	for (int i = 0; i < hp.stem.pairs.size(); ++i) {
		if (hp.stem.pairs[i].r1.ix_residue != -1) {
			hp.stem.pairs[i].r1.ix_residue += ix_residue;
		}
		if (hp.stem.pairs[i].r2.ix_residue != -1) {
			hp.stem.pairs[i].r2.ix_residue += ix_residue;
		}
		hp.stem.pairs[i].r1.ix_chain = 0;
		hp.stem.pairs[i].r2.ix_chain = 0;
		sstream.str(string());
		sstream.clear();
		sstream << res_pos + atoi(hp.stem.pairs[i].r1.residue_position.c_str());
		hp.stem.pairs[i].r1.residue_position = sstream.str();
		sstream.str(string());
		sstream.clear();
		sstream << res_pos + atoi(hp.stem.pairs[i].r2.residue_position.c_str());
		hp.stem.pairs[i].r2.residue_position = sstream.str();
	}
	for (int i = 0; i < hp.head.residues.size(); ++i) {
		if (hp.head.residues[i].r1.ix_residue != -1) {
			hp.head.residues[i].r1.ix_residue += ix_residue;
		}
		hp.head.residues[i].r1.ix_chain = 0;
		sstream.str(string());
		sstream.clear();
		sstream << res_pos + atoi(hp.head.residues[i].r1.residue_position.c_str());
		hp.head.residues[i].r1.residue_position = sstream.str();
	}

}

bool FormAvg::SwapRNA(
	cRNAStructure &rna1,
	cRNAStructure &rna2,
	int &length1,
	int &length2,
	int mode)
{
	bool swapped = false;
	if (mode == STEM_R1 || mode == STEM_R2) {
		length1 = rna1.mHairpins[0].stem.pairs.size();
		length2 = rna2.mHairpins[0].stem.pairs.size();
	}
	else {
		length1 = rna1.mHairpins[0].head.residues.size();
		length2 = rna2.mHairpins[0].head.residues.size();
	}

	if (length2 > length1) {
		cRNAStructure rSwap = rna1;		int lSwap = length1;
		rna1 = rna2;					length1 = length2;
		rna2 = rSwap;					length2 = lSwap;
		swapped = true;
	}

	return swapped;
}

void CalcAver(tPdbAtom &atom1, tPdbAtom &atom2, t3DCoords &aver, int mode, int pos)
{
	if (mode == STEM_R1 || mode == STEM_R2) {
		aver.x = (pos * atom1.coords.x + atom2.coords.x) / (pos + 1);
		aver.y = (pos * atom1.coords.y + atom2.coords.y) / (pos + 1);
		aver.z = (pos * atom1.coords.z + atom2.coords.z) / (pos + 1);
	}
	else {
		aver.x = (atom1.coords.x + atom2.coords.x) / 2;
		aver.y = (atom1.coords.y + atom2.coords.y) / 2;
		aver.z = (atom1.coords.z + atom2.coords.z) / 2;
	}
}

void FormAvg::MergePairToRNA(
	unsigned idx,
	cRNAStructure &outRna)
{
	if (Parameters::GetInstance().GetMergeAlgorithm() == Parameters::Merge::ALTERNATING) {
		AlternatingMergePairToRNA(idx, outRna);
	}
	else {
		PatternMergePairToRNA(idx, outRna);
	}
}

void FormAvg::PatternMergePairToRNA(
	unsigned idx,
	cRNAStructure &outRna)
{

	cout << "PM";
	stringstream sstream;
	int ixS = 0;
	string sxS;

	vector<string> chain;

	tPairResiude tp, *tp1, *tp2, *tp1other = NULL, *tp2other = NULL;
	tPdbModel m;//, *m1, *m2;
	int length1, length2;

	cHairpin	hp;
	tPdbModel	newModel;
	tPdbModel	newModelNt;
	tPdbChain	newCh;
	tPdbChain	newChNt;

	cRNAStructure &rna1 = mRnaQVect[idx];
	cRNAStructure &rna2 = mRnaDBVect[idx];

	if (rna1.GetHairpinsCount() == 0 && rna1.GetHairpinsCount() == 0) {
		return;
	}
	if (rna1.GetHairpinsCount() == 0) {
		outRna = rna2;
		return;
	}
	if (rna2.GetHairpinsCount() == 0) {
		outRna = rna1;
		return;
	}

	tPdbChain ch;
	string residue;
	tPdbAtom atom, atom1, atom2;
	tPdbAtom atomIt;

	t3DCoords aver;
	t3DCoords trans;
	int transDir;
	int transMult;

	trans.x = trans.y = trans.z = 0.0;

	for (int mode = 0; mode < 3; ++mode) {
		bool swapped = SwapRNA(rna1, rna2, length1, length2, mode);

		int k = 0;
		int ii = -1;
		for (int i = 0; i < length1; ++i) {
			++ii;

			int pos;

			if (length2 > 0) {
				double step = (double)(length2 - 1) / (double)(length1 - 1);
				int j = roundDouble(i * step);

				switch (mode) {
				case STEM_R1:
					pos = length1 - 1 - i;
					tp1 = &rna1.mHairpins[0].stem.pairs[length1 - 1 - i].r1;
					tp2 = &rna2.mHairpins[0].stem.pairs[length2 - 1 - j].r1;

					tp1other = &rna1.mHairpins[0].stem.pairs[length1 - 1 - i].r2;
					tp2other = &rna2.mHairpins[0].stem.pairs[length2 - 1 - j].r2;
					break;
				case HEAD:
					pos = i;
					tp1 = &rna1.mHairpins[0].head.residues[i].r1;
					tp2 = &rna2.mHairpins[0].head.residues[j].r1;
					break;
				case STEM_R2:
					pos = i;
					tp1 = &rna1.mHairpins[0].stem.pairs[i].r2;
					tp2 = &rna2.mHairpins[0].stem.pairs[j].r2;

					tp1other = &rna1.mHairpins[0].stem.pairs[i].r1;
					tp2other = &rna2.mHairpins[0].stem.pairs[j].r1;
					break;
				}

				if (tp1->ix_residue == -1 && tp2->ix_residue == -1) {
					switch (mode) {
					case STEM_R1:
						hp.stem.pairs.insert(hp.stem.pairs.begin(), tPair());
						break;
					case HEAD:
						hp.head.residues.push_back(tPair());
						break;
					case STEM_R2:
						if (hp.stem.pairs[ii].r1.ix_residue == -1) {
							hp.stem.pairs.erase(hp.stem.pairs.begin() + ii);
							--ii;
						}
						break;
					}
					continue;
				}

				if (tp1->ix_residue != -1 && tp2->ix_residue != -1) {

					atom1 = rna1.GetAtom(*tp1);
					atom2 = rna2.GetAtom(*tp2);

					CalcAver(atom1, atom2, aver);

					if ((tp1other->ix_residue != -1 && tp2other->ix_residue != -1) ||
						(tp1other->ix_residue == -1 && tp2other->ix_residue == -1)) {
						if (!swapped) {
							residue = rna1.getPrimarySeq()[tp1->ix_residue];
							tp = *tp1;
							m = rna1.mModels[0];
							atom = atom1;
							// translation from atom1
							transDir = 1;
						}
						else {
							residue = rna2.getPrimarySeq()[tp2->ix_residue];
							tp = *tp2;
							m = rna2.mModels[0];
							atom = atom2;
							// translation from atom2 
							transDir = 2;
						}
					}
					else if (tp1other->ix_residue != -1) {
						residue = rna1.getPrimarySeq()[tp1->ix_residue];
						tp = *tp1;
						m = rna1.mModels[0];
						atom = atom1;
						// translation from atom1
						transDir = 1;
					}
					else if (tp2other->ix_residue != -1) {
						residue = rna2.getPrimarySeq()[tp2->ix_residue];
						tp = *tp2;
						m = rna2.mModels[0];
						atom = atom2;
						// translation from atom2 
						transDir = 2;
					}

					CalcTrans(atom, aver, trans);
					transMult = 1;

				}
				else if (tp1->ix_residue != -1) {
					if (tp2other->ix_residue != -1 && pos % 2 == 1) {
						if (mode == STEM_R1) {
							hp.stem.pairs.insert(hp.stem.pairs.begin(), tPair());
						}
						continue;
					}
					residue = rna1.getPrimarySeq()[tp1->ix_residue];
					tp = *tp1;
					m = rna1.mModels[0];
					atom = rna1.GetAtom(*tp1);

					if (transDir == 1) {
						transMult = 1;
					}
					else {
						transMult = -1;
					}
				}
				else if (tp2->ix_residue != -1) {
					if (tp1other->ix_residue != -1 && pos % 2 == 0) {
						if (mode == STEM_R1) {
							hp.stem.pairs.insert(hp.stem.pairs.begin(), tPair());
						}
						continue;
					}
					residue = rna2.getPrimarySeq()[tp2->ix_residue];
					tp = *tp2;
					m = rna2.mModels[0];
					atom = rna2.GetAtom(*tp2);

					if (transDir == 2) {
						transMult = 1;
					}
					else {
						transMult = -1;
					}
				}
			}
			else {
				switch (mode) {
				case STEM_R1:
					tp = rna1.mHairpins[0].stem.pairs[length1 - 1 - i].r1;
					break;
				case HEAD:
					tp = rna1.mHairpins[0].head.residues[i].r1;
					break;
				case STEM_R2:
					tp = rna1.mHairpins[0].stem.pairs[i].r2;
					break;
				}

				if (tp.ix_residue == -1) {
					switch (mode) {
					case STEM_R1:
						hp.stem.pairs.insert(hp.stem.pairs.begin(), tPair());
						break;
					case HEAD:
						hp.head.residues.push_back(tPair());
						break;
					case STEM_R2:
						if (hp.stem.pairs[ii].r1.ix_residue == -1) {
							hp.stem.pairs.erase(hp.stem.pairs.begin() + ii);
							--ii;
						}
						break;
					}
					continue;
				}
				atom = rna1.GetAtom(tp);

				residue = rna1.getPrimarySeq()[tp.ix_residue];
				m = rna1.mModels[0];
			}

			sstream.str(string());
			sstream.clear();
			sstream << ixS + 1;
			sxS = sstream.str();

			chain.push_back(residue);
			ch = m[tp.ix_chain];

			bool found = false;
			while (!found){
				for (; k < ch.size(); ++k) {
					if (ch[k].residue_num == atom.residue_num) {
						found = true;
						break;
					}
				}
				if (!found) k = 0;
			}
			for (; k < ch.size(); ++k) {
				if (ch[k].residue_num == atom.residue_num) {
					atomIt = ch[k];
					atomIt.residue_num = sxS;
					if (length2 > 0) {
						atomIt.coords.x += transMult * trans.x;
						atomIt.coords.y += transMult * trans.y;
						atomIt.coords.z += transMult * trans.z;
					}
					newCh.push_back(atomIt);
				}
				else {
					//break;
				}
			}

			atom.residue_num = sxS;
			if (length2 > 0) {
				atom.coords.x += transMult * trans.x;
				atom.coords.y += transMult * trans.y;
				atom.coords.z += transMult * trans.z;
			}
			newChNt.push_back(atom);

			sstream.str(string());
			sstream.clear();
			sstream << ixS;
			tp.residue_position = sstream.str();
			tp.ix_residue = ixS++;
			tp.ix_chain = 0;

			switch (mode) {
			case STEM_R1:
				hp.stem.pairs.insert(hp.stem.pairs.begin(), tPair());
				hp.stem.pairs[0].r1 = tp;
				break;
			case HEAD:
				hp.head.residues.push_back(tPair());
				hp.head.residues.back().r1 = tp;
				break;
			case STEM_R2:
				hp.stem.pairs[ii].r2 = tp;
				break;
			}
		}
	}

	outRna.mPrimarySeqChains.insert(outRna.mPrimarySeqChains.end(), chain);

	newModel.push_back(newCh);
	newModelNt.push_back(newChNt);

	outRna.mModels.push_back(newModel);
	outRna.mModelsNt.push_back(newModelNt);

	string pdbID = rna1.mPdbID + rna2.mPdbID;
	sstream.str(string());
	sstream.clear();
	sstream << mData.mOptPairs[idx].q;
	pdbID.append(sstream.str());
	pdbID.append("_");
	sstream.str(string());
	sstream.clear();
	sstream << mData.mOptPairs[idx].db;
	pdbID.append(sstream.str());
	outRna.mPdbID = pdbID;
	outRna.mHairpins.push_back(hp);
}

void FormAvg::AlternatingMergePairToRNA(
	unsigned idx,
	cRNAStructure &outRna)
{
	stringstream sstream;
	int ixS = 0;
	string sxS;

	vector<string> chain;

	tPairResiude tp, *tp1, *tp2, *tp1other = NULL, *tp2other = NULL;
	tPdbModel m;//, *m1, *m2;
	int length1, length2;

	cHairpin	hp;
	tPdbModel	newModel;
	tPdbModel	newModelNt;
	tPdbChain	newCh;
	tPdbChain	newChNt;

	cRNAStructure &rna1 = mRnaQVect[idx];
	cRNAStructure &rna2 = mRnaDBVect[idx];

	if (rna1.GetHairpinsCount() == 0 && rna1.GetHairpinsCount() == 0) {
		return;
	}
	if (rna1.GetHairpinsCount() == 0) {
		outRna = rna2;
		return;
	}
	if (rna2.GetHairpinsCount() == 0) {
		outRna = rna1;
		return;
	}

	tPdbChain ch;
	string residue;
	tPdbAtom atom, atom1, atom2;
	tPdbAtom atomIt;

	t3DCoords aver;
	t3DCoords trans;
	int transDir = 1;
	int transMult = 1;

	trans.x = trans.y = trans.z = 0.0;

	for (int mode = 0; mode < 3; ++mode) {
		SwapRNA(rna1, rna2, length1, length2, mode);

		int k = 0;
		int ii = -1;
		for (int i = 0; i < length1; ++i) {
			++ii;

			int pos;

			if (length2 > 0) {
				double step = (double)( max(1, length2 - 1)) / (double)(max (1, length1 - 1));
				int j = roundDouble(i * step);

				switch (mode) {
				case STEM_R1:
					pos = length1 - 1 - i;
					tp1 = &rna1.mHairpins[0].stem.pairs[length1 - 1 - i].r1;
					tp2 = &rna2.mHairpins[0].stem.pairs[length2 - 1 - j].r1;

					tp1other = &rna1.mHairpins[0].stem.pairs[length1 - 1 - i].r2;
					tp2other = &rna2.mHairpins[0].stem.pairs[length2 - 1 - j].r2;
					break;
				case HEAD:
					pos = i;
					tp1 = &rna1.mHairpins[0].head.residues[i].r1;
					tp2 = &rna2.mHairpins[0].head.residues[j].r1;
					break;
				case STEM_R2:
					pos = i;
					tp1 = &rna1.mHairpins[0].stem.pairs[i].r2;
					tp2 = &rna2.mHairpins[0].stem.pairs[j].r2;

					tp1other = &rna1.mHairpins[0].stem.pairs[i].r1;
					tp2other = &rna2.mHairpins[0].stem.pairs[j].r1;
					break;
				}

				if (tp1->ix_residue == -1 && tp2->ix_residue == -1) {
					switch (mode) {
					case STEM_R1:
						hp.stem.pairs.insert(hp.stem.pairs.begin(), tPair());
						break;
					case HEAD:
						hp.head.residues.push_back(tPair());
						break;
					case STEM_R2:
						if (hp.stem.pairs[ii].r1.ix_residue == -1) {
							hp.stem.pairs.erase(hp.stem.pairs.begin() + ii);
							--ii;
						}
						break;
					}
					continue;
				}

				if (tp1->ix_residue != -1 && tp2->ix_residue != -1) {

					atom1 = rna1.GetAtom(*tp1);
					atom2 = rna2.GetAtom(*tp2);

					CalcAver(atom1, atom2, aver);

					if ((tp1other != NULL &&  tp1other->ix_residue != -1 && tp2other != NULL  && tp2other->ix_residue != -1) ||
						((tp1other == NULL || tp1other->ix_residue == -1) && (tp2other != NULL || tp2other->ix_residue == -1))) {
						if (pos % 2) {
							residue = rna1.getPrimarySeq()[tp1->ix_residue];
							tp = *tp1;
							m = rna1.mModels[0];
							atom = atom1;
							// translation from atom1
							transDir = 1;
						}
						else {
							residue = rna2.getPrimarySeq()[tp2->ix_residue];
							tp = *tp2;
							m = rna2.mModels[0];
							atom = atom2;
							// translation from atom2 
							transDir = 2;
						}
					}
					else if (tp1other->ix_residue != -1) {
						residue = rna1.getPrimarySeq()[tp1->ix_residue];
						tp = *tp1;
						m = rna1.mModels[0];
						atom = atom1;
						// translation from atom1
						transDir = 1;
					}
					else if (tp2other->ix_residue != -1) {
						residue = rna2.getPrimarySeq()[tp2->ix_residue];
						tp = *tp2;
						m = rna2.mModels[0];
						atom = atom2;
						// translation from atom2 
						transDir = 2;
					}

					CalcTrans(atom, aver, trans);
					transMult = 1;

				}
				else if (tp1->ix_residue != -1) {
					if (tp2other->ix_residue != -1 && pos % 2 == 1) {
						if (mode == STEM_R1) {
							hp.stem.pairs.insert(hp.stem.pairs.begin(), tPair());
						}
						continue;
					}
					residue = rna1.getPrimarySeq()[tp1->ix_residue];
					tp = *tp1;
					m = rna1.mModels[0];
					atom = rna1.GetAtom(*tp1);

					if (transDir == 1) {
						transMult = 1;
					}
					else {
						transMult = -1;
					}
				}
				else if (tp2->ix_residue != -1) {
					if (tp1other->ix_residue != -1 && pos % 2 == 0) {
						if (mode == STEM_R1) {
							hp.stem.pairs.insert(hp.stem.pairs.begin(), tPair());
						}
						continue;
					}
					residue = rna2.getPrimarySeq()[tp2->ix_residue];
					tp = *tp2;
					m = rna2.mModels[0];
					atom = rna2.GetAtom(*tp2);

					if (transDir == 2) {
						transMult = 1;
					}
					else {
						transMult = -1;
					}
				}
			}
			else {
				switch (mode) {
				case STEM_R1:
					tp = rna1.mHairpins[0].stem.pairs[length1 - 1 - i].r1;
					break;
				case HEAD:
					tp = rna1.mHairpins[0].head.residues[i].r1;
					break;
				case STEM_R2:
					tp = rna1.mHairpins[0].stem.pairs[i].r2;
					break;
				}

				if (tp.ix_residue == -1) {
					switch (mode) {
					case STEM_R1:
						hp.stem.pairs.insert(hp.stem.pairs.begin(), tPair());
						break;
					case HEAD:
						hp.head.residues.push_back(tPair());
						break;
					case STEM_R2:
						if (hp.stem.pairs[ii].r1.ix_residue == -1) {
							hp.stem.pairs.erase(hp.stem.pairs.begin() + ii);
							--ii;
						}
						break;
					}
					continue;
				}
				atom = rna1.GetAtom(tp);

				residue = rna1.getPrimarySeq()[tp.ix_residue];
				m = rna1.mModels[0];
			}

			sstream.str(string());
			sstream.clear();
			sstream << ixS + 1;
			sxS = sstream.str();

			chain.push_back(residue);
			ch = m[tp.ix_chain];

			/*
			 * After the representants of the residues are merged (each residue is represented
			 * by one atom), it is time to merge the rest of the atoms.
			 */
			bool found = false;
			while (!found){
				//identify where in the chain the atoms of a given residue (atom.residue_num) start
				for (; k < ch.size(); ++k) {
					if (ch[k].residue_num == atom.residue_num) {
						found = true;
						break;
					}
				}
				if (!found) k = 0;
			}
			for (; k < ch.size(); ++k) {
				//do while the residue number does not change
				if (ch[k].residue_num == atom.residue_num) {
					atomIt = ch[k];
					atomIt.residue_num = sxS;
					if (length2 > 0) {
						atomIt.coords.x += transMult * trans.x;
						atomIt.coords.y += transMult * trans.y;
						atomIt.coords.z += transMult * trans.z;
					}
					newCh.push_back(atomIt);
				}
				else {
					//break;
				}
			}

			atom.residue_num = sxS;
			if (length2 > 0) {
				atom.coords.x += transMult * trans.x;
				atom.coords.y += transMult * trans.y;
				atom.coords.z += transMult * trans.z;
			}
			newChNt.push_back(atom);

			sstream.str(string());
			sstream.clear();
			sstream << ixS;
			tp.residue_position = sstream.str();
			tp.ix_residue = ixS++;
			tp.ix_chain = 0;

			switch (mode) {
			case STEM_R1:
				hp.stem.pairs.insert(hp.stem.pairs.begin(), tPair());
				hp.stem.pairs[0].r1 = tp;
				break;
			case HEAD:
				hp.head.residues.push_back(tPair());
				hp.head.residues.back().r1 = tp;
				break;
			case STEM_R2:
				hp.stem.pairs[ii].r2 = tp;
				break;
			}
		}
	}

	outRna.mPrimarySeqChains.insert(outRna.mPrimarySeqChains.end(), chain);

	newModel.push_back(newCh);
	newModelNt.push_back(newChNt);

	outRna.mModels.push_back(newModel);
	outRna.mModelsNt.push_back(newModelNt);

	string pdbID = rna1.mPdbID + rna2.mPdbID;
	sstream.str(string());
	sstream.clear();
	sstream << mData.mOptPairs[idx].q;
	pdbID.append(sstream.str());
	pdbID.append("_");
	sstream.str(string());
	sstream.clear();
	sstream << mData.mOptPairs[idx].db;
	pdbID.append(sstream.str());
	outRna.mPdbID = pdbID;
	outRna.mHairpins.push_back(hp);
}
