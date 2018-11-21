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
#include "algo.h"
#include "parallel.h"
#include "parameters.h"

#include <boost/algorithm/string.hpp>

void ConvertHairpinToRNA(
	cRNAStructure &rna, 
	cHairpin &hp, 
	cRMSDStructExt &s, 
	int alignActual[3][2], 
	int alignActualConverted[3])
{
	int ixS = 0;
	tPairResiude tp;

	for (int i = hp.stem.pairs.size() - 1; i >= 0; i--)	
	{
		if (hp.stem.pairs[i].r1.ix_residue >= 0)
		{
			tp = hp.stem.pairs[i].r1;
			for (int j = 0; j < 3; j++)
			{
				if (alignActual[j][0] == tp.ix_chain && alignActual[j][1] == tp.ix_residue)
				{
					alignActualConverted[j] = ixS;
				}
			}
			// new if
			if (tp.ix_chain > 0) {
				continue;
			}
			tPdbAtom at = rna.GetAtom(tp);
			s.coord[ixS++] = cRMSD3DCoord(at.coords.x, at.coords.y, at.coords.z);
			s.aa.push_back(boost::algorithm::trim_copy(at.res_name));
			s.sse += hp.stem.pairs[i].pair_type;
			s.residueInfo.push_back(tp);

		}
	}
	for (unsigned int i = 0; i < hp.head.residues.size(); i++)
	{
		tp = hp.head.residues[i].r1;
		for (int j = 0; j < 3; j++)
		{
			if (alignActual[j][0] == tp.ix_chain && alignActual[j][1] == tp.ix_residue)
			{
				alignActualConverted[j] = ixS;
			}
		}
		// new if
		if (tp.ix_chain > 0) {
			continue;
		}
		tPdbAtom at = rna.GetAtom(tp);
		s.coord[ixS++] = cRMSD3DCoord(at.coords.x, at.coords.y, at.coords.z);
		s.aa.push_back(boost::algorithm::trim_copy(at.res_name));
		s.sse += "ooooo";
		s.residueInfo.push_back(tp);
	}
	for (int i = 0; i < (int)hp.stem.pairs.size(); i++)
	{
		if (hp.stem.pairs[i].r2.ix_residue >= 0)
		{
			tp = hp.stem.pairs[i].r2;
			for (int j = 0; j < 3; j++)
			{
				if (alignActual[j][0] == tp.ix_chain && alignActual[j][1] == tp.ix_residue)
				{
					alignActualConverted[j] = ixS;
				}
			}
			// new if
			if (tp.ix_chain > 0) {
				continue;
			}
			tPdbAtom at = rna.GetAtom(tp);
			s.coord[ixS++] = cRMSD3DCoord(at.coords.x, at.coords.y, at.coords.z);
			s.aa.push_back(boost::algorithm::trim_copy(at.res_name));
			s.sse += hp.stem.pairs[i].pair_type;
			s.residueInfo.push_back(tp);
		}
	}

	s.lengthTemp = ixS;
	s.SetName(rna.getId());
}

pair<int, double> IdentifyL2NearestNeighbor(
	cRMSDStructExt &s, 
	cRMSD3DCoord &c, 
	int ixFrom, 
	int ixTo)
{
	double min_dist = ResMatch::MAX;
	int min_index = -1;
	for (int i = ixFrom; i < ixTo; i++)
	{
		double dist = sqrt(pow(s.coord[i][0] - c[0], 2) + pow(s.coord[i][1] - c[1], 2) + pow(s.coord[i][2] - c[2], 2));
		if (dist < min_dist)
		{
			min_dist = dist;
			min_index = i;
		}
	}

	return make_pair(min_index, min_dist);
}

pair<int, double> IdentifyL2NearestNeighbor(
	cHairpin &hp, 
	cRNAStructure &rna, 
	cRMSD3DCoord &c)
{
	double min_dist = ResMatch::MAX;
	int min_index = -1;

	int ixS = 0;
	double dist;
	tPdbAtom at;

	for (int i = 0; i < (int)hp.stem.pairs.size(); i++)
	{
		if (hp.stem.pairs[i].r1.ix_residue >= 0)
		{
			ixS++;
			at = rna.GetAtom(hp.stem.pairs[i].r1);			
			dist = sqrt(pow(at.coords.x - c[0], 2) + pow(at.coords.y - c[1], 2) + pow(at.coords.z - c[2], 2));
			if (dist < min_dist)
			{
				min_dist = dist;
				min_index = ixS;
			}
		}
	}
	for (unsigned int i = 0; i < hp.head.residues.size(); i++)
	{
		ixS++;
		at = rna.GetAtom(hp.stem.pairs[i].r1);			
		dist = sqrt(pow(at.coords.x - c[0], 2) + pow(at.coords.y - c[1], 2) + pow(at.coords.z - c[2], 2));
		if (dist < min_dist)
		{
			min_dist = dist;
			min_index = ixS;
		}
	}
	for (int i = hp.stem.pairs.size() - 1; i >= 0; i--)
	{
		if (hp.stem.pairs[i].r2.ix_residue >= 0)
		{
			ixS++;
			at = rna.GetAtom(hp.stem.pairs[i].r1);			
			dist = sqrt(pow(at.coords.x - c[0], 2) + pow(at.coords.y - c[1], 2) + pow(at.coords.z - c[2], 2));
			if (dist < min_dist)
			{
				min_dist = dist;
				min_index = ixS;
			}
		}
	}

	return make_pair(min_index, min_dist);
}

void SortArray3(int arr[])
{
	int temp[3];
	for (int i = 0; i < 3; temp[i] = arr[i++]);
	arr[0] = min(temp[0], min(temp[1], temp[2]));
	arr[2] = max(temp[0], max(temp[1], temp[2]));
	for (int i = 0; i < 3; i++)
	{
		if (temp[i] != arr[0] && temp[i] != arr[2])
		{
			arr[1] = temp[i];
		}
	}
}

double Normalize(
	double score, 
	int lenQ, int 
	lenDB, 
	int cntHpQ, 
	int cntHpDB, 
	int cntNN)
{
	return score / min(lenQ,lenDB) * (1 + abs((int)lenQ - (int)lenDB) / (double)min(lenQ,lenDB)) / max(cntNN, 1);
}

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
	double scoreLimit)
{
	sResMatch result;

	cRMSDMatrix mRotDB, mRotQ;
	cRMSD3DCoord transQ, transDB;

	for(int i = 0; i < 3; i++)
	{
		mRotQ[i][i] = 1;
	}

	double scoreRMSD = rmsd(&sQ, &sDB, align, mRotDB, transQ, transDB);
	
	// The 3 passed to Normalize function is the number of nearest neighbor pairs being in distance at most NNDist. Since we
	// approximate using triplets the number is 3.
	if (Normalize(scoreRMSD, sQ.length, sDB.length, rnaQ.GetHairpinsCount(), rnaDB.GetHairpinsCount(), 3)/params.earlyTermination >= scoreLimit)
	{
		result.score = ResMatch::MAX;
		
		return result;
	}

	//since the hairpin will be conveted into cRNAStructure, the indexes of aligned residues need to be converted too 
	//because, e.g., the hairpin can be composed of two chains which would corrupt the process of aligning corresponding sections (see below)
	//moreover, the order of items in stem does not correspond to the order of the residues (0-th item in the stem is the neck)
	int alignActualConvertedQ[3], alignActualConvertedDB[3]; 

	cRMSDStructExt sHpQ, sHpDB;
	sHpQ.length = rnaQ.GetLengthNt();
	sHpDB.length = rnaDB.GetLengthNt();
	sHpQ.coord = new cRMSD3DCoord[sHpQ.length];
	sHpDB.coord = new cRMSD3DCoord[sHpDB.length];

	//aligned arrays need to be recomputed since the PDB atoms change their indexes (only a part of rnaX is converted to sHpX)
	ConvertHairpinToRNA(rnaQ, hpQ, sHpQ, alignActualQ, alignActualConvertedQ);
	ConvertHairpinToRNA(rnaDB, hpDB, sHpDB, alignActualDB, alignActualConvertedDB);	
	//store the pair triplet responsible for the superposition
	for (int i = 0; i < 3; i++)
	{
		tPair p;
		result.superpositionTriplet[i].r1.ix_chain = alignActualQ[i][0]; 
		result.superpositionTriplet[i].r1.ix_residue = alignActualQ[i][1];
		result.superpositionTriplet[i].r2.ix_chain = alignActualDB[i][0]; 
		result.superpositionTriplet[i].r2.ix_residue = alignActualDB[i][1];		
	}
	result.rot = mRotDB;
	result.trans = transQ-mRotDB*transDB;
	cRMSDStructExt sDBRot(&sHpDB, result.rot, result.trans);
	
	double score = 0;

	///*	
	int cntNNUpToDist = 0;
	for (int i = 0; i < sHpQ.lengthTemp; i++)
	{
		pair<int, double> nn = IdentifyL2NearestNeighbor(sDBRot, sHpQ.coord[i], 0, sDBRot.lengthTemp);
		if (nn.second < params.NNDist)
		{
			cntNNUpToDist++;
		}
		//score += nn.second;
		if (nn.first >= 0)
		{
			if (sHpQ.aa[i] == sDBRot.aa[nn.first]) //the names are taken from PDB format where they take 4 positions
			{
				nn.second *= params.identicalLetterModificator;
			}
			if (sHpQ.sse.length() >= (i+1)*5 && sDBRot.sse.length() >= (nn.first+1)*5 && sHpQ.sse.substr(i*5,5) == sDBRot.sse.substr(nn.first*5,5)) //the names are taken from PDB format where they take 4 positions
			{
				nn.second *= params.identicalPairTypeModificator;
			}
			score += nn.second;
		}
	}

	for (int i = 0; i < sDBRot.lengthTemp; i++)
	{
		pair<int, double> nn = IdentifyL2NearestNeighbor(sHpQ, sDBRot.coord[i], 0, sHpQ.lengthTemp);
		if (nn.second < params.NNDist)
		{
			cntNNUpToDist++;
		}
		//score += nn.second;
		if (nn.first >= 0)
		{
			if (sDBRot.aa[i] == sHpQ.aa[nn.first]) //the names are taken from PDB format where they take 4 positions
			{
				nn.second *= params.identicalLetterModificator;
			}
			if (sDBRot.sse.length() >= (i+1)*5 && sHpQ.sse.length() >= (nn.first+1)*5 && sDBRot.sse.substr(i*5,5) == sHpQ.sse.substr(nn.first*5,5)) //the names are taken from PDB format where they take 4 positions
			{
				nn.second *= params.identicalPairTypeModificator;
			}
			score += nn.second;
		}
	}
	cntNNUpToDist /= 2; //many NN's distances were accounted twice so this is only rough modification so that values of this symmetric measure correspond to the assymetric one
	score /= 2;

	result.score = Normalize(score, sHpQ.lengthTemp, sHpDB.lengthTemp, 1, 1, cntNNUpToDist);
	
	return result;
}

sResMatch MatchHairpins(
	cRNAStructure &rnaQ, 
	cHairpin &hpQ, 
	cRNAStructure &rnaDB, 
	cHairpin &hpDB, 
	sParams params, 
	double scoreLimit)
{
	tPair neckQ, neckDB;
	int ixNeckQ = -1, ixNeckDB = -1;
	neckQ.r1.ix_residue = neckQ.r2.ix_residue = neckDB.r1.ix_residue = neckDB.r2.ix_residue = -1;

	//we don't use neckQ = hpQ.stem.pairs[0] since there is a possibility that the first pair is not complete
	//therefore we use the for cycle
	for (int i = 0; i < (int)hpQ.stem.pairs.size(); i++)
	{
		if (hpQ.stem.pairs[i].r1.ix_residue >= 0 && hpQ.stem.pairs[i].r2.ix_residue >= 0)
		{
			neckQ = hpQ.stem.pairs[i];
			ixNeckQ = i;
			break;
		}
	}
	for (int i = 0; i < (int)hpDB.stem.pairs.size(); i++)
	{
		if (hpDB.stem.pairs[i].r1.ix_residue >= 0 && hpDB.stem.pairs[i].r2.ix_residue >= 0)
		{
			neckDB = hpDB.stem.pairs[i];
			ixNeckDB = i;
			break;
		}
	}

	sResMatch resOpt;
	resOpt.score = ResMatch::MAX;
	const unsigned int alignLen = 3;	
	cRMSDAlign align[2] = {cRMSDAlign(alignLen), cRMSDAlign(alignLen)};
	cRMSDAlign alignBest(alignLen);
	cRMSDStructExt sQ, sDB;

	sQ.length = sDB.length = 3;
	sQ.coord = new cRMSD3DCoord[sQ.length];
	sDB.coord = new cRMSD3DCoord[sDB.length];

	for (unsigned int i = 0; i < 3; i++)
	{
		align[0][i].x = align[0][i].y = i;
		align[1][i].x = i;
		align[1][i].y = 2 - i;
	}

	tPdbAtom atom;

	/**/
#if PARALLEL
	tbb::mutex resMutex;
	ParallelHairpinMatch phm(rnaQ, rnaDB, hpQ, hpDB, ixNeckQ, ixNeckDB, params, scoreLimit, resOpt, resMutex);

	tbb::parallel_for(
		tbb::blocked_range<size_t>(0, hpQ.head.residues.size()), phm);
#else	
	// serial
	for (unsigned int i = 0; i < hpQ.head.residues.size(); i++)
	{
		for (unsigned int j = 0; j < hpDB.head.residues.size(); j++)
		{	
			for (int kQ = ixNeckQ; kQ < ixNeckQ+params.neckMaxShift; kQ++)
			{ 
				if (kQ < 0) continue;

				if (kQ > (int)hpQ.stem.pairs.size() - 1 || hpQ.stem.pairs[kQ].r1.ix_residue < 0 || hpQ.stem.pairs[kQ].r2.ix_residue < 0) continue;
				
				for (int kDB = ixNeckDB; kDB < ixNeckDB+params.neckMaxShift; kDB++)
				{//e.g. in case of 1zih and 1zif we need to "move the neck" to obtain correct alignment
					
					if (kDB < 0) continue;

					if (kDB > (int)hpDB.stem.pairs.size() - 1 || hpDB.stem.pairs[kDB].r1.ix_residue < 0 || hpDB.stem.pairs[kDB].r2.ix_residue < 0) continue;

					
					tPairResiude prQ[3] = {hpQ.stem.pairs[kQ].r1, hpQ.head.residues[i].r1, hpQ.stem.pairs[kQ].r2};
					tPairResiude prDB[3] = {hpDB.stem.pairs[kDB].r1, hpDB.head.residues[j].r1, hpDB.stem.pairs[kDB].r2};

					atom = rnaQ.GetAtom(prQ[0]); sQ.coord[0] = cRMSD3DCoord(atom.coords.x, atom.coords.y, atom.coords.z);
					atom = rnaQ.GetAtom(prQ[1]); sQ.coord[1] = cRMSD3DCoord(atom.coords.x, atom.coords.y, atom.coords.z);
					atom = rnaQ.GetAtom(prQ[2]); sQ.coord[2] = cRMSD3DCoord(atom.coords.x, atom.coords.y, atom.coords.z);

					atom = rnaDB.GetAtom(prDB[0]); sDB.coord[0] = cRMSD3DCoord(atom.coords.x, atom.coords.y, atom.coords.z);
					atom = rnaDB.GetAtom(prDB[1]); sDB.coord[1] = cRMSD3DCoord(atom.coords.x, atom.coords.y, atom.coords.z);
					atom = rnaDB.GetAtom(prDB[2]); sDB.coord[2] = cRMSD3DCoord(atom.coords.x, atom.coords.y, atom.coords.z);		

					//alignActual arrays contain positions of aligned atoms in the original PDB (RNA) structure
					//thus one can access aligned atoms' properties (such as coordinates)
					int alignActualQ[3][2] = {{prQ[0].ix_chain, prQ[0].ix_residue},{prQ[1].ix_chain, prQ[1].ix_residue}, {prQ[2].ix_chain, prQ[2].ix_residue}};
					int alignActualDB[3][2] = {{prDB[0].ix_chain, prDB[0].ix_residue},{prDB[1].ix_chain, prDB[1].ix_residue}, {prDB[2].ix_chain, prDB[2].ix_residue}};

					//according to one of the bioinformatics referees, it does not make sense to try align 5'->3' agains 3'->5'
					// in such case, l should be upper limited by 1, else by 2 ... HOWEVER, this does not hold when one of those structures is average structure!!!
					for (int l = 0; l < 2; l++)
					//for (int l = 0; l < 1; l++)
					{
						sResMatch res = MatchHairpinsBasedOnTriplets(sQ, rnaQ, hpQ, sDB, rnaDB, hpDB, params, align[l], alignActualQ, alignActualDB, min(scoreLimit, resOpt.score));
						if (res.score < resOpt.score)
						{
							resOpt = res;
							/*
							sQ.name = rnaQ.getId();
							sDB.name = rnaDB.getId();
							string nameQ = sQ.name.substr(sQ.name.find_last_of("/")+1);
							nameQ = nameQ.substr(0, nameQ.find_first_of("."));
							string nameDB = sDB.name.substr(sDB.name.find_last_of("/")+1);
							nameDB = nameDB.substr(0, nameDB.find_first_of("."));
							string name = nameQ + string("_") + nameDB + string(".vmd");
							ofstream ofs(name.c_str());
							rmsd(&sQ, &sDB, align[l], &ofs);
							*/
						}						
					}
				}
			}
		}
	}
#endif	

	int ixStepQ = hpQ.stem.pairs.size() / params.noHeadShiftStepRatioDiviser;
	int ixStepDB = hpDB.stem.pairs.size() / params.noHeadShiftStepRatioDiviser;

	if (ixStepQ == 0) ixStepQ = 1;
	if (ixStepDB == 0) ixStepDB = 1;

	if (resOpt.score == ResMatch::MAX && hpQ.GetLength() > 3 && hpDB.GetLength() > 3 )
	{ //if one of the HPs does not have the head and number of residues in HPs is higher than 3 (otherwise the HPs are optimally superposable)
		if ( (neckQ.r1.ix_residue < 0 || neckQ.r2.ix_residue < 0) || (neckDB.r1.ix_residue < 0 || neckDB.r2.ix_residue < 0) )
		{ //the stem does not have any pair
			for (int i1 = 0; i1 < (int)hpQ.stem.pairs.size(); i1+=ixStepQ)
			{
				if (hpQ.stem.pairs[i1].r1.ix_residue < 0) continue; 
				for (int j1 = i1+1; j1 < (int)hpQ.stem.pairs.size(); j1+=ixStepQ)
				{
					if (hpQ.stem.pairs[j1].r1.ix_residue < 0) continue; 
					for (int k1 = j1+1; k1 < (int)hpQ.stem.pairs.size(); k1+=ixStepQ)
					{
						if (hpQ.stem.pairs[k1].r1.ix_residue < 0) continue; 
						for (int i2 = 0; i2 < (int)hpDB.stem.pairs.size(); i2+=ixStepDB)
						{
							if (hpDB.stem.pairs[i2].r1.ix_residue < 0) continue; 
							for (int j2 = i2+1; j2 < (int)hpDB.stem.pairs.size(); j2+=ixStepDB)
							{
								if (hpDB.stem.pairs[j2].r1.ix_residue < 0) continue; 
								for (int k2 = j2+1; k2 < (int)hpDB.stem.pairs.size(); k2+=ixStepDB)
								{									
									if (hpDB.stem.pairs[k2].r1.ix_residue < 0) continue; 

									tPairResiude prQ[3] = {hpQ.stem.pairs[i1].r1, hpQ.stem.pairs[j1].r1, hpQ.stem.pairs[k1].r1};
									tPairResiude prDB[3] = {hpDB.stem.pairs[i2].r1, hpDB.stem.pairs[j2].r1, hpDB.stem.pairs[k2].r1};

									atom = rnaQ.GetAtom(prQ[0]); sQ.coord[0] = cRMSD3DCoord(atom.coords.x, atom.coords.y, atom.coords.z);
									atom = rnaQ.GetAtom(prQ[1]); sQ.coord[1] = cRMSD3DCoord(atom.coords.x, atom.coords.y, atom.coords.z);
									atom = rnaQ.GetAtom(prQ[2]); sQ.coord[2] = cRMSD3DCoord(atom.coords.x, atom.coords.y, atom.coords.z);

									atom = rnaDB.GetAtom(prDB[0]); sDB.coord[0] = cRMSD3DCoord(atom.coords.x, atom.coords.y, atom.coords.z);
									atom = rnaDB.GetAtom(prDB[1]); sDB.coord[1] = cRMSD3DCoord(atom.coords.x, atom.coords.y, atom.coords.z);
									atom = rnaDB.GetAtom(prDB[2]); sDB.coord[2] = cRMSD3DCoord(atom.coords.x, atom.coords.y, atom.coords.z);		

									int alignActualQ[3][2] = {{prQ[0].ix_chain, prQ[0].ix_residue},{prQ[1].ix_chain, prQ[1].ix_residue}, {prQ[2].ix_chain, prQ[2].ix_residue}};
									int alignActualDB[3][2] = {{prDB[0].ix_chain, prDB[0].ix_residue},{prDB[1].ix_chain, prDB[1].ix_residue}, {prDB[2].ix_chain, prDB[2].ix_residue}};

									for (int k = 0; k < 1; k++)
									//for (int k = 0; k < 2; k++)
									{ //in this situation is no reason for using reverse alignment 
										//scoreBest = min(scoreBest, MatchHairpinsBasedOnTriplets(sQ, rnaQ, hpQ, sDB, rnaDB, hpDB, align[k], alignActualQ, alignActualDB));
										sResMatch res = MatchHairpinsBasedOnTriplets(sQ, rnaQ, hpQ, sDB, rnaDB, hpDB, params, align[k], alignActualQ, alignActualDB, min(scoreLimit, resOpt.score));
										if (res.score < resOpt.score)
										{
											resOpt = res;
											/*
											sQ.name = rnaQ.getId();
											sDB.name = rnaDB.getId();
											string nameQ = sQ.name.substr(sQ.name.find_last_of("/")+1);
											nameQ = nameQ.substr(0, nameQ.find_first_of("."));
											string nameDB = sDB.name.substr(sDB.name.find_last_of("/")+1);
											nameDB = nameDB.substr(0, nameDB.find_first_of("."));
											string name = nameQ + string("_") + nameDB + string(".vmd");
											ofstream ofs(name.c_str());
											rmsd(&sQ, &sDB, align[k], &ofs);
											*/
										}						
									}
								}
							}
						}
					}
				}
			}
		}
		else if ( 
			(neckQ.r1.ix_residue >= 0 && neckQ.r2.ix_residue >= 0 && hpQ.head.residues.size() == 0) 
			|| 
			(neckDB.r1.ix_residue >= 0 && neckDB.r2.ix_residue >= 0 && hpDB.head.residues.size() == 0) 
			)
		{ //one of the HPs does not have head (but both have stems) => only stems will be used for superposition
			for (int i1 = 0; i1 < (int)hpQ.stem.pairs.size(); i1+=ixStepQ)
			{
				if (i1 == ixNeckQ)
				{
					continue;
				}
				for (int i2 = 0; i2 < (int)hpDB.stem.pairs.size(); i2+=ixStepDB)
				{
					if (i2 == ixNeckDB)
					{
						continue;
					}

					tPairResiude prQ[3] = {neckQ.r1, hpQ.stem.pairs[i1].r1, neckQ.r2};
					if (prQ[1].ix_residue < 0) prQ[1] = hpQ.stem.pairs[i1].r2;
					tPairResiude prDB[3] = {neckDB.r1, hpDB.stem.pairs[i2].r1, neckDB.r2};
					if (prDB[1].ix_residue < 0) prDB[1] = hpDB.stem.pairs[i2].r2;
					
					if (prQ[1].ix_residue < 0 || prDB[1].ix_residue < 0) {
						continue;
					}

					atom = rnaQ.GetAtom(prQ[0]); sQ.coord[0] = cRMSD3DCoord(atom.coords.x, atom.coords.y, atom.coords.z);
					atom = rnaQ.GetAtom(prQ[1]); sQ.coord[1] = cRMSD3DCoord(atom.coords.x, atom.coords.y, atom.coords.z);
					atom = rnaQ.GetAtom(prQ[2]); sQ.coord[2] = cRMSD3DCoord(atom.coords.x, atom.coords.y, atom.coords.z);

					atom = rnaDB.GetAtom(prDB[0]); sDB.coord[0] = cRMSD3DCoord(atom.coords.x, atom.coords.y, atom.coords.z);
					atom = rnaDB.GetAtom(prDB[1]); sDB.coord[1] = cRMSD3DCoord(atom.coords.x, atom.coords.y, atom.coords.z);
					atom = rnaDB.GetAtom(prDB[2]); sDB.coord[2] = cRMSD3DCoord(atom.coords.x, atom.coords.y, atom.coords.z);		

					int alignActualQ[3][2] = {{prQ[0].ix_chain, prQ[0].ix_residue},{prQ[1].ix_chain, prQ[1].ix_residue}, {prQ[2].ix_chain, prQ[2].ix_residue}};
					int alignActualDB[3][2] = {{prDB[0].ix_chain, prDB[0].ix_residue},{prDB[1].ix_chain, prDB[1].ix_residue}, {prDB[2].ix_chain, prDB[2].ix_residue}};

					//according to one of the bioinformatics referees, it does not make sense to try align 5'->3' agains 3'->5'
					// in such case, k should be upper limited by 1, else by 2
					//for (int k = 0; k < 2; k++)
					for (int k = 0; k < 1; k++)
					{
						//scoreBest = min(scoreBest, MatchHairpinsBasedOnTriplets(sQ, rnaQ, hpQ, sDB, rnaDB, hpDB, align[k], alignActualQ, alignActualDB));
						sResMatch res = MatchHairpinsBasedOnTriplets(sQ, rnaQ, hpQ, sDB, rnaDB, hpDB, params, align[k], alignActualQ, alignActualDB, min(scoreLimit, resOpt.score));
						if (res.score < resOpt.score)
						{
							resOpt = res;
						}						
					}
				}
			}			
		}
		
	}

	return resOpt;
}

IxsScore::IxsScore(double xscore, int xi, int xj): score(xscore), i(xi), j(xj){};

bool compareIxsScore(IxsScore i1, IxsScore i2)
{
	return (i1.score < i2.score);
}

// resolved
sResMatch Match(
	cRNAStructure rnaQ, 
	cRNAStructure rnaDB, 
	sParams params,
	bool extraSuperposition,
	bool nnListNeeded, 
	double bestScore)
{
	boost::posix_time::ptime tStart = boost::posix_time::microsec_clock::local_time();

	int cntHpQ = (int)rnaQ.GetHairpinsCount();
	int cntHpDB = (int)rnaDB.GetHairpinsCount();
	sResMatch **ress = new sResMatch*[cntHpQ];
	vector<sIxsScore> ressOrderedIxs;
	
	for (int i = 0; i < cntHpQ; i++) ress[i] = new sResMatch[cntHpDB];
	sResMatch resOpt;

	

	resOpt.score = ResMatch::MAX;
		
	for (int i = 0; i < cntHpQ; i++)
	{
		for (int j = 0; j < cntHpDB; j++)
		{	
			cHairpin hpQ = rnaQ.GetHairpin(i);
			cHairpin hpDB = rnaDB.GetHairpin(j);
			sResMatch res;			
			ress[i][j] = res = MatchHairpins(rnaQ, hpQ, rnaDB, hpDB, params, bestScore);
			ressOrderedIxs.push_back(sIxsScore(res.score, i, j));			
			if (ressOrderedIxs.size() == params.topKForMultiGSSU || (ressOrderedIxs.size() > params.topKForMultiGSSU && res.score < ressOrderedIxs[params.topKForMultiGSSU-1].score))
			{
				std::sort(ressOrderedIxs.begin(), ressOrderedIxs.end(), compareIxsScore);
				bestScore = ressOrderedIxs[params.topKForMultiGSSU-1].score;
			}
			//if (res.score < resOpt.score) resOpt = res;
		}
	}

	resOpt.score = ResMatch::MAX;
	
	for (int iOrdered = 0; iOrdered < min(params.topKForMultiGSSU, (int)ressOrderedIxs.size()); iOrdered++)
	//for (int i = 0; i < cntHpQ; i++)
	{ // For topKForMultiGSSU best GSSUs from Q (having minimum distance to some GSSU in the DB structure compute the full score
		int i = ressOrderedIxs[iOrdered].i;
		/*
		 * For each of the GSSUs in Q structure identify the most similar GSSU in DB structure 
		 */
		int jMin = 0;
		//for (int j = 1; j < cntHpDB; j++) if (ress[i][j].score < ress[i][jMin].score) jMin = j;
		jMin = ressOrderedIxs[iOrdered].j;
		/*
		 * Align i-th GSSU in Q with minJ-th GSSU in DB and successively align the neighborhood GSSU with penalties corresponding to necessary transformations
		 */

		//int cntAlignedGSSUsInDirection1, cntAlignedGSSUsInDirection2;
		//cntAlignedGSSUsInDirection1 = 1 + min(i, jMin) + min(cntHpQ-i-1, cntHpDB-jMin-1);
		//cntAlignedGSSUsInDirection2 = 1 + min(i, cntHpDB-jMin-1) + min(cntHpQ-i-1, jMin);
		vector< pair<int, int> > alignedGSSUs;
		vector< double > alignedGSSUScores;
		vector< double > rotationPenalties;

		for (int ixReverse = 0; ixReverse < 2; ixReverse++)
		{
			//if (cntAlignedGSSUsInDirection1 > cntAlignedGSSUsInDirection2 && ixReverse == 1) continue;
			//if (cntAlignedGSSUsInDirection1 < cntAlignedGSSUsInDirection2 && ixReverse == 0) continue;
			//ixReverse designates the direction in which the GSSUs should be aligned (solves the situation when
			//last GSSU of Q is initialy alligned with first of DB)

			int cntAlignedGSSUs = 1;
			double scoreAggI = ress[i][jMin].score;

			alignedGSSUs.push_back(pair<int,int>(i, jMin));
			alignedGSSUScores.push_back(scoreAggI);
			int iPrev = i, jPrev = jMin;
		
			for (int l = 0; l < 2; l++)
			{
				//if l==0, GSSUs "upstream" from i,jMin are being aligned, else downstream
				int k;
				l == 0 ? k = 1 : k = -1;
				while (true)
				{
					int iCurrent, jCurrent;
					if (ixReverse == 0)
					{
						
						iCurrent = i + k;
						jCurrent = jMin + k;
					}
					else
					{
						iCurrent = i + k;
						jCurrent = jMin - k;
					}
					if ( iCurrent < 0 || jCurrent < 0 || iCurrent >= cntHpQ || jCurrent >= cntHpDB ) break;			

					//first we need to compute the full score if it has not been computed yet
					//!!! maybe the i,j and not maximum i+k,jMin+k distance should be also recomputed since it might not correspond with the real distance because of the heuristics application !!!
					if (ress[iCurrent][jCurrent].score == ResMatch::MAX) 
					{
						cHairpin hpQ = rnaQ.GetHairpin(iCurrent);
						cHairpin hpDB = rnaDB.GetHairpin(jCurrent);
						sResMatch tempScore;
						ress[iCurrent][jCurrent] = tempScore = MatchHairpins(rnaQ, hpQ, rnaDB, hpDB, params, ResMatch::MAX);
					}

					//rotation penalty is computed as the difference in angles needed to rotate the transposed molecule into i-th state and (i+1)-st state
					double penaltyRot = 1;
					penaltyRot += fabs(norm2((ress[iPrev][jPrev].rot - ress[iCurrent][jCurrent].rot)/2));
					//cout << penaltyRot << "-";
					if (ress[iCurrent][jCurrent].score != ResMatch::MAX)
					{
					//if one of the GSSUs contains only one pair of nts then no score will be assessed to the GSSU pair and this pair should be skipped
						rotationPenalties.push_back(penaltyRot);
						scoreAggI += ress[iCurrent][jCurrent].score;// * penaltyRot * params.rotationPenalty; //the rotation penalty concerns two pairs of GSSU and therefore will be added to the overall score
						alignedGSSUs.push_back(pair<int,int>(iCurrent, jCurrent));
						alignedGSSUScores.push_back(ress[iCurrent][jCurrent].score);
						cntAlignedGSSUs++;
						
						iPrev = iCurrent;
						jPrev = jCurrent;
					}

					// TODO: translation penalty
					//Translation penalty is actually included in the GSSU-GSSU distances, because necessary translation is caused only by
					//different lengths of the GSSUs. The rotation penalty, which is the other part of the penalty is included in previous instruction.
					l == 0 ? k++ : k--;
				}				
			}
			double penalty = 1;
			for (int ixPenalty = 0; ixPenalty < rotationPenalties.size(); ixPenalty++) penalty += rotationPenalties[ixPenalty] * params.rotationPenalty;
			penalty /= (1+rotationPenalties.size());
			//cout << "score: " << scoreAggI << " rotation: " << penalty;
			scoreAggI *= penalty;
			int cntNonAlignedGSSUs = cntHpQ + cntHpDB - 2*cntAlignedGSSUs;
			int cntMaxNonAlignedGSSUs = cntHpQ + cntHpDB - 1;
			int cntMaxAlignedGSSUs = min(cntHpQ, cntHpDB);
			//penalty for different GSSU size		
			//cout << " aligned " << cntAlignedGSSUs << "non-aligned: " << cntNonAlignedGSSUs;
			//scoreAggI = (scoreAggI / cntAlignedGSSUs) * pow(-log(1-(double)cntNonAlignedGSSUs/cntMaxNonAlignedGSSUs)+1, 6);
			//scoreAggI = (scoreAggI / cntAlignedGSSUs) * (-log(cntNonAlignedGSSUs/(double)cntMaxNonAlignedGSSUs)+1) * (-log(1-cntAlignedGSSUs/(double)cntMaxAlignedGSSUs)+1);	
			scoreAggI = (scoreAggI / cntAlignedGSSUs) * (1 + cntNonAlignedGSSUs)*2;
			//cout << " nonaligned penalty:" << pow(-log(1-(double)cntNonAlignedGSSUs/cntMaxNonAlignedGSSUs)+1, 6);
			//cout << "overal: " << scoreAggI;
			if (scoreAggI < resOpt.score)
			{
				//cout << " best ";
				resOpt = ress[i][jMin];
				resOpt.score = scoreAggI;
				resOpt.alignedGSSUs = alignedGSSUs;
				resOpt.alignedScores = alignedGSSUScores;
				resOpt.rotationPenalties = rotationPenalties;
			}
		}
	}

	for (int i = 0; i < cntHpQ; i++) delete[] ress[i];
	delete [] ress;

	if (resOpt.score != ResMatch::MAX)
	{
		//let's recompute the nn-list to refer to the whole structures not only the best GSSU

		vector<tPdbAtom> resQ = rnaQ.GetAllNtResidues();
		vector<tPdbAtom> resDB = rnaDB.GetAllNtResidues();
		vector<tPdbAtom> resDBRot = rnaDB.GetAllNtResidues();

		//transform coord vector of the DB structure
		cRMSD3DCoord coord, coordRot;
		for (int i = 0; i < resDB.size(); i++)
		{
			coord[0] = resDB[i].coords.x; coord[1] = resDB[i].coords.y;	coord[2] = resDB[i].coords.z;
			coordRot = resOpt.rot * coord + resOpt.trans;			
			resDBRot[i].coords.x = coordRot[0]; resDBRot[i].coords.y = coordRot[1]; resDBRot[i].coords.z = coordRot[2];			
		}

		//return resOpt;

		//if the extra step with superposition based on nearest neigbours is requested it has to be done before the nearest neighbors
		//are discovered and inserted into the result since in the result the nearest neighbors appear together with their mutual distances 
		//which will change after the extra RMSD superposition step
		if (extraSuperposition)
		{
			//first let's precompute the distance matrix to avoid duplicate distance computations
			double **dists = new double*[resQ.size()];
			for (int i = 0; i < resQ.size(); i++) dists[i] = new double[resDBRot.size()];
			for (int i = 0 ; i < resQ.size(); i++)
				for (int j = 0 ; j < resDBRot.size(); j++) 
				{
					dists[i][j] = sqrt( pow(resQ[i].coords.x - resDBRot[j].coords.x, 2) + pow(resQ[i].coords.y - resDBRot[j].coords.y, 2) + pow(resQ[i].coords.z - resDBRot[j].coords.z, 2) );
				}
			

			vector< pair<int, int> > nnForRMSD;
			for (int i = 0 ; i < resQ.size(); i++)
			{
				double min_dist = ResMatch::MAX;
				int min_index = -1;
				t3DCoords coordsQ = resQ[i].coords;
				int j;
				for (j = 0 ; j < resDBRot.size(); j++)
				{
					if (min_dist > dists[i][j])
					{
						min_dist = dists[i][j];
						min_index = j;
					}
				}
				if (min_index == -1) continue;				
				
				//check whether resQ[i] is also the NN of resDBRot[min_index], if not do not consider this pair for nnForRMSD
				for (j = 0; j < resQ.size(); j++)
				{	
					if (dists[j][min_index] < min_dist) break;
				}
				if (j < resQ.size()) continue;
				

				nnForRMSD.push_back(make_pair(i, min_index));					
			}

			if (nnForRMSD.size() > 3)
			{
				//global RMSD over all GSSUs
				cRMSDStructExt sQAux, sDBAux;
				sQAux.length = resQ.size();
				sDBAux.length = resDB.size();
				sQAux.coord = new cRMSD3DCoord[sQAux.length];
				sDBAux.coord = new cRMSD3DCoord[sDBAux.length];
				for (int i = 0 ; i < resQ.size(); i++){ sQAux.coord[i][0] = resQ[i].coords.x; sQAux.coord[i][1] = resQ[i].coords.y; sQAux.coord[i][2] = resQ[i].coords.z; };
				for (int i = 0 ; i < resDB.size(); i++){ sDBAux.coord[i][0] = resDB[i].coords.x; sDBAux.coord[i][1] = resDB[i].coords.y; sDBAux.coord[i][2] = resDB[i].coords.z; };
				cRMSDAlign align = cRMSDAlign(nnForRMSD.size());
				for (int i = 0; i < nnForRMSD.size(); i++) { align[i].x = nnForRMSD[i].first; align[i].y = nnForRMSD[i].second; }
				cRMSD3DCoord transQ, transDB;
				cRMSDMatrix rotDB;
				int lenMin = sDBAux.length, lenMax = sQAux.length;
				if (sQAux.length <= sDBAux.length) {lenMin = sQAux.length; lenMax = sDBAux.length; }
				rmsd(&sQAux, &sDBAux, align, rotDB, transQ, transDB); 
				resOpt.rot = rotDB;
				resOpt.trans = transQ-rotDB*transDB;

				//since the rotation and translation matrix probably changed, the coordination of the resulting superposed strucutre also needs to be modified
				//however, we do not modify the resOpt.score simply because the whole SETTER procedure would have to be repeated
				cRMSD3DCoord coord, coordRot;
				for (int i = 0; i < resDB.size(); i++)
				{
					coord[0] = resDB[i].coords.x; coord[1] = resDB[i].coords.y;	coord[2] = resDB[i].coords.z;
					coordRot = resOpt.rot * coord + resOpt.trans;			
					resDBRot[i].coords.x = coordRot[0]; resDBRot[i].coords.y = coordRot[1]; resDBRot[i].coords.z = coordRot[2];			
				}
			}
		}

		//identify NN for each query nt and write it to the result
		
		for (int k = 0; k < 2; k++)
		{
			resOpt.nn[k].clear();
			vector<tPdbAtom> *pRes1, *pRes2;

			if (k == 0) { pRes1 = &resQ; pRes2 = &resDBRot; }
			else { pRes1 = &resDBRot; pRes2 = &resQ; };

			for (int i = 0 ; i < (*pRes1).size(); i++)
			{
				double min_dist = ResMatch::MAX;
				int min_index = -1;
				t3DCoords coordsQ = (*pRes1)[i].coords;
				for (int j = 0 ; j < (*pRes2).size(); j++)
				{
					double dist = sqrt( pow(coordsQ.x - (*pRes2)[j].coords.x, 2) + pow(coordsQ.y - (*pRes2)[j].coords.y, 2) + pow(coordsQ.z - (*pRes2)[j].coords.z, 2) );
					if (min_dist > dist)
					{
						min_dist = dist;
						min_index = j;
					}
				}
				if (min_index == -1) {
					//cout << "el van baszva" << endl;
					continue;
				}
				pair<string, string> nnPair = make_pair(string(1, (*pRes1)[i].chain) + "-" + (*pRes1)[i].residue_num, string(1, (*pRes2)[min_index].chain) + "-" + (*pRes2)[min_index].residue_num);
				resOpt.nn[k].push_back(make_pair(nnPair, min_dist));
			}
		}
		
		
	}

	resOpt.timeStart = tStart;
	resOpt.timeEnd = boost::posix_time::microsec_clock::local_time();
	return resOpt;
}

// resolved
sResMatch Match(
	cRNAStructure rnaQ, 
	cRNAStructure rnaDB, 
	sParams params, 
	sResMatch ***outRess, 
	vector<GSSUPair> &outOptPair, 
	bool nnListNeeded, 
	double bestScore)
{
	int cntHpQ = (int)rnaQ.GetHairpinsCount();
	int cntHpDB = (int)rnaDB.GetHairpinsCount();
	sResMatch **ress = new sResMatch*[cntHpQ];
	vector<sIxsScore> ressOrderedIxs;
	
	for (int i = 0; i < cntHpQ; i++) ress[i] = new sResMatch[cntHpDB];
	sResMatch resOpt;
	resOpt.score = ResMatch::MAX;
		
	for (int i = 0; i < cntHpQ; i++)
	{
		for (int j = 0; j < cntHpDB; j++)
		{	
			cHairpin hpQ = rnaQ.GetHairpin(i);
			cHairpin hpDB = rnaDB.GetHairpin(j);
			sResMatch res;			
			(*outRess)[i][j] = ress[i][j] = res = MatchHairpins(rnaQ, hpQ, rnaDB, hpDB, params, bestScore);
			ressOrderedIxs.push_back(sIxsScore(res.score, i, j));			
			if (ressOrderedIxs.size() == params.topKForMultiGSSU || (ressOrderedIxs.size() > params.topKForMultiGSSU && res.score < ressOrderedIxs[params.topKForMultiGSSU-1].score))
			{
				std::sort(ressOrderedIxs.begin(), ressOrderedIxs.end(), compareIxsScore);
				bestScore = ressOrderedIxs[params.topKForMultiGSSU-1].score;
			}
			//if (res.score < resOpt.score) resOpt = res;
		}
	}

	resOpt.score = ResMatch::MAX;
	
	// the optimal ordered sequence of GSSU pairs
	vector<pair<int, int>> optPairs;

	for (int iOrdered = 0; iOrdered < min(params.topKForMultiGSSU, (int)ressOrderedIxs.size()); iOrdered++)
	//for (int i = 0; i < cntHpQ; i++)
	{ // For topKForMultiGSSU best GSSUs from Q (having minimum distance to some GSSU in the DB structure compute the full score
		vector<GSSUPair> pairs;
		int i = ressOrderedIxs[iOrdered].i;
		/*
		 * For each of the GSSUs in Q structure identify the most similar GSSU in DB structure 
		 */
		int jMin = 0;
		//for (int j = 1; j < cntHpDB; j++) if (ress[i][j].score < ress[i][jMin].score) jMin = j;
		jMin = ressOrderedIxs[iOrdered].j;
		/*
		 * Align i-th GSSU in Q with minJ-th GSSU in DB and successively align the neighborhood GSSU with penalties corresponding to necessary transformations
		 */

		//int cntAlignedGSSUsInDirection1, cntAlignedGSSUsInDirection2;
		//cntAlignedGSSUsInDirection1 = 1 + min(i, jMin) + min(cntHpQ-i-1, cntHpDB-jMin-1);
		//cntAlignedGSSUsInDirection2 = 1 + min(i, cntHpDB-jMin-1) + min(cntHpQ-i-1, jMin);
		vector< pair<int, int> > alignedGSSUs;
		vector< double > alignedGSSUScores;
		vector< double > rotationPenalties;

		for (int ixReverse = 0; ixReverse < 2; ixReverse++)
		{
			//if (cntAlignedGSSUsInDirection1 > cntAlignedGSSUsInDirection2 && ixReverse == 1) continue;
			//if (cntAlignedGSSUsInDirection1 < cntAlignedGSSUsInDirection2 && ixReverse == 0) continue;
			//ixReverse designates the direction in which the GSSUs should be aligned (solves the situation when
			//last GSSU of Q is initialy alligned with first of DB)

			int cntAlignedGSSUs = 1;
			double scoreAggI = ress[i][jMin].score;

			pairs.push_back(GSSUPair(i, jMin));
			alignedGSSUs.push_back(pair<int,int>(i, jMin));
			alignedGSSUScores.push_back(scoreAggI);
			int iPrev = i, jPrev = jMin;

			for (int l = 0; l < 2; l++)
			{
				//if l==0, GSSUs "upstream" from i,jMin are being aligned, else downstream
				int k;
				l == 0 ? k = 1 : k = -1;
				while (true)
				{
					int iCurrent, jCurrent;
					if (ixReverse == 0)
					{
						
						iCurrent = i + k;
						jCurrent = jMin + k;
					}
					else
					{
						iCurrent = i + k;
						jCurrent = jMin - k;
					}
					if ( iCurrent < 0 || jCurrent < 0 || iCurrent >= cntHpQ || jCurrent >= cntHpDB ) break;			

					//first we need to compute the full score if it has not been computed yet
					//!!! maybe the i,j and not maximum i+k,jMin+k distance should be also recomputed since it might not correspond with the real distance because of the heuristics application !!!
					if (ress[iCurrent][jCurrent].score == ResMatch::MAX) 
					{
						cHairpin hpQ = rnaQ.GetHairpin(iCurrent);
						cHairpin hpDB = rnaDB.GetHairpin(jCurrent);
						sResMatch tempScore;
						ress[iCurrent][jCurrent] = tempScore = MatchHairpins(rnaQ, hpQ, rnaDB, hpDB, params, ResMatch::MAX);
					}

					//rotation penalty is computed as the difference in angles needed to rotate the transposed molecule into i-th state and (i+1)-st state
					double penaltyRot = 1;
					penaltyRot += fabs(norm2((ress[iPrev][jPrev].rot - ress[iCurrent][jCurrent].rot)/2));
					//cout << penaltyRot << "-";
					if (ress[iCurrent][jCurrent].score != ResMatch::MAX)
					{
					//if one of the GSSUs contains only one pair of nts then no score will be assessed to the GSSU pair and this pair should be skipped
						pairs.push_back(GSSUPair(iCurrent, jCurrent));
						rotationPenalties.push_back(penaltyRot);
						scoreAggI += ress[iCurrent][jCurrent].score;// * penaltyRot * params.rotationPenalty; //the rotation penalty concerns two pairs of GSSU and therefore will be added to the overall score
						alignedGSSUs.push_back(pair<int,int>(iCurrent, jCurrent));
						alignedGSSUScores.push_back(ress[iCurrent][jCurrent].score);
						cntAlignedGSSUs++;

						iPrev = iCurrent;
						jPrev = jCurrent;
					}

					// TODO: translation penalty
					//Translation penalty is actually included in the GSSU-GSSU distances, because necessary translation is caused only by
					//different lengths of the GSSUs. The rotation penalty, which is the other part of the penalty is included in previous instruction.
					l == 0 ? k++ : k--;
				}				
			}
			double penalty = 1;
			for (int ixPenalty = 0; ixPenalty < rotationPenalties.size(); ixPenalty++) penalty += rotationPenalties[ixPenalty] * params.rotationPenalty;
			penalty /= (1+rotationPenalties.size());
			//cout << "score: " << scoreAggI << " rotation: " << penalty;
			scoreAggI *= penalty;
			int cntNonAlignedGSSUs = cntHpQ + cntHpDB - 2*cntAlignedGSSUs;
			int cntMaxNonAlignedGSSUs = cntHpQ + cntHpDB - 1;
			int cntMaxAlignedGSSUs = min(cntHpQ, cntHpDB);
			//penalty for different GSSU size		
			//cout << " aligned " << cntAlignedGSSUs << "non-aligned: " << cntNonAlignedGSSUs;
			//scoreAggI = (scoreAggI / cntAlignedGSSUs) * pow(-log(1-(double)cntNonAlignedGSSUs/cntMaxNonAlignedGSSUs)+1, 6);
			//scoreAggI = (scoreAggI / cntAlignedGSSUs) * (-log(cntNonAlignedGSSUs/(double)cntMaxNonAlignedGSSUs)+1) * (-log(1-cntAlignedGSSUs/(double)cntMaxAlignedGSSUs)+1);
			//penalty for different GSSU size		
			scoreAggI = (scoreAggI / cntAlignedGSSUs) * (1 + cntNonAlignedGSSUs)*2;
			//cout << " nonaligned penalty:" << pow(-log(1-(double)cntNonAlignedGSSUs/cntMaxNonAlignedGSSUs)+1, 6);
			//cout << "overal: " << scoreAggI;
			if (scoreAggI < resOpt.score)
			{
				outOptPair = pairs;
				resOpt = ress[i][jMin];
				resOpt.score = scoreAggI;
				resOpt.alignedGSSUs = alignedGSSUs;
				resOpt.alignedScores = alignedGSSUScores;
				resOpt.rotationPenalties = rotationPenalties;
			}
		}
	}

	for (int i = 0; i < cntHpQ; i++) delete[] ress[i];
	delete [] ress;

	if (resOpt.score != ResMatch::MAX)
	{
		//let's recompute the nn-list to refer to the whole structures not only the best GSSU

		vector<tPdbAtom> resQ = rnaQ.GetAllNtResidues();
		vector<tPdbAtom> resDB = rnaDB.GetAllNtResidues();
		vector<tPdbAtom> resDBRot = rnaDB.GetAllNtResidues();

		//transform coord vector of the DB structure
		cRMSD3DCoord coord, coordRot;
		for (int i = 0; i < resDB.size(); i++)
		{
			coord[0] = resDB[i].coords.x; coord[1] = resDB[i].coords.y;	coord[2] = resDB[i].coords.z;
			coordRot = resOpt.rot * coord + resOpt.trans;			
			resDBRot[i].coords.x = coordRot[0]; resDBRot[i].coords.y = coordRot[1]; resDBRot[i].coords.z = coordRot[2];			
		}

		return resOpt;

		//identify NN for each query nt and write it to the result
		
		vector< pair<int, int> > nnForRMSD;
		for (int k = 0; k < 2; k++)
		{
			resOpt.nn[k].clear();
			vector<tPdbAtom> *pRes1, *pRes2;

			if (k == 0) { pRes1 = &resQ; pRes2 = &resDBRot; }
			else { pRes1 = &resDBRot; pRes2 = &resQ; };

			for (int i = 0 ; i < (*pRes1).size(); i++)
			{
				double min_dist = ResMatch::MAX;
				int min_index = -1;
				t3DCoords coordsQ = (*pRes1)[i].coords;
				for (int j = 0 ; j < (*pRes2).size(); j++)
				{
					double dist = sqrt( pow(coordsQ.x - (*pRes2)[j].coords.x, 2) + pow(coordsQ.y - (*pRes2)[j].coords.y, 2) + pow(coordsQ.z - (*pRes2)[j].coords.z, 2) );
					if (min_dist > dist)
					{
						min_dist = dist;
						min_index = j;
					}
				}
				if (min_index == -1) {
					continue;
				}
				pair<string, string> nnPair = make_pair(string(1, (*pRes1)[i].chain) + "-" + (*pRes1)[i].residue_num, string(1, (*pRes2)[min_index].chain) + "-" + (*pRes2)[min_index].residue_num);
				resOpt.nn[k].push_back(make_pair(nnPair, min_dist));
				
				//if (min_dist <= 4) nnForRMSD.push_back(make_pair(i, min_index));
			}
		}
	}

	return resOpt;
}