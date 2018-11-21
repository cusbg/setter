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
//#define NDEBUG

//#define VLD_FORCE_ENABLE 1
//#include <vld.h>

#include "algo.h"
#include "dirent.h"
#include "rapidxml-1.13/rapidxml.hpp"
#include "rapidxml-1.13/rapidxml_print.hpp"
#include "input.h"
#include "make_aver_rna.h"
#include "multi_align.h"

#include "parameters.h"
#include "print_to_pdb.h"

#include "tbb/task_scheduler_init.h"

#include <set>
#include <fstream>
#include <cctype>
#include <time.h>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <string.h>

#include <fstream>

#include <boost/serialization/vector.hpp>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp> 
#include <boost/date_time.hpp>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>


using namespace std;
namespace po = boost::program_options;

//! This function is responsible for processing the output datasets.
/*!
	Generates average RNA from family and then alignes all RNAs from dataset.
	If an actually processed RNA from [dataset file] is presented in [input file], the
	RNA is omitted from the average RNA computation (we consider two RNA-s the same, if
	the pdb names are the same and the chains are the same). This procedure is meant
	for testing purposes.

	\param family	reference to input family
	\param dataset	reference to dataset of RNAs
	*/

void CompairDatasetToFamily(Input &family, Input &dataset)
{
	Parameters &param = Parameters::GetInstance();
	vector<int> rnaVectInFamily;     // indices in family
	vector<int> rnaVectNotInFamily;  // indices in dataset

	for (int i = 0; i < dataset.mInputSize; ++i) {
		string dataId = dataset.mOriginalMap[i];
		int j = 0;
		for (; j < family.mInputSize; ++j) {
			string pdbId = family.mOriginalMap[j];
			if (dataId == pdbId) {
				rnaVectInFamily.push_back(j);
				break;
			}
		}

		if (j == family.mInputSize) {
			rnaVectNotInFamily.push_back(i);
		}
	}

	/*
	 * We form a family from all the structures in the family
	 * and align all the remaining structures (dataset) to it one by one
	 * to get the distances.
	 */
	if (rnaVectNotInFamily.size() > 0) {
		MultiAlign *nj = MultiAlign::Factory(param.GetAlgorithm(), family);
		MAResult familyAver = nj->GetResult();

		vector<sResMatch> output(rnaVectNotInFamily.size());
		ParallelMainBlock parallelMainBlock(
			familyAver.rna, rnaVectNotInFamily, dataset.mInputRnaVect, dataset.mOriginalMap, output);

		tbb::parallel_for(
			tbb::blocked_range<size_t>(0, rnaVectNotInFamily.size()), parallelMainBlock);
		parallelMainBlock.PrintOutput();

		delete nj;
	}

	/*
	 * We iteratively remove structures from the family and
	 * each time we build an average structure without that
	 * structure. Then, we measure the distance of that
	 * removed structure to the average. The point is that
	 * measuring distance to the average is correct only
	 * when that structure was not used in the process of building
	 * the average structure.*/
	for (int i = 0; i < rnaVectInFamily.size(); ++i) {
		int pos = rnaVectInFamily[i];
		cRNAStructure &testRNA = *family.mInputRnaVect[pos];

		Input subFamily(family);
		subFamily.RemoveRNA(pos);

		MultiAlign *nj = MultiAlign::Factory(param.GetAlgorithm(), subFamily);
		MAResult aver = nj->GetResult();

		sParams params = GetGlobalParams();
		double score = CalcScore(testRNA, aver.rna, params).score;

		cout << family.mOriginalMap[pos] << " - aver : "
			<< score << endl;

		delete nj;
	}
}

double erf(double x)
{
	// constants
	double a1 = 0.254829592;
	double a2 = -0.284496736;
	double a3 = 1.421413741;
	double a4 = -1.453152027;
	double a5 = 1.061405429;
	double p = 0.3275911;

	// Save the sign of x
	int sign = 1;
	if (x < 0)
		sign = -1;
	x = fabs(x);

	// A&S formula 7.1.26
	double t = 1.0 / (1.0 + p * x);
	double y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * exp(-x * x);

	return sign * y;
}

bool IsStaggeredPosition(vector<vector<int>> alIxs, int ixColumn)
{
	int cntFilled = 0;

	for (int ix = 0; ix < alIxs.size(); ix++)
		if (alIxs[ix][ixColumn] >= 0) cntFilled++;

	return cntFilled > 1 ? false : true;
}

bool IsAssignableIntoGroup(t3DCoords coords, vector<pair<int, t3DCoords>> group, double threshold)
{
	bool assignable = true;
	for (int ixGroup = 0; ixGroup < group.size(); ixGroup++)
	{
		if (group[ixGroup].second.DistL2(coords) > threshold)
		{
			assignable = false;
			break;
		}
	}
	return assignable;
}

t3DCoords GetCoordiantesFromAlignment(vector<cRNAStructure*> input, vector<vector<int>> alIxs, int ix)
{
	int ixNonNullPosition;
	for (ixNonNullPosition = 0; ixNonNullPosition < alIxs.size(); ixNonNullPosition++) if (alIxs[ixNonNullPosition][ix] >= 0) break;

	assert(input[ixNonNullPosition]->GetAllNtResidues().size() > alIxs[ixNonNullPosition][ix]);

	return input[ixNonNullPosition]->GetAllNtResidues()[alIxs[ixNonNullPosition][ix]].coords;
}

void MergeColumns(vector<vector<int>> &alIxs, vector<int> &newColumn, int ixCol)
{
	for (int ix = 0; ix < alIxs.size(); ix++)
	{
		assert(newColumn[ix] == -1 || alIxs[ix][ixCol] == -1);
		if (alIxs[ix][ixCol] != -1) newColumn[ix] = alIxs[ix][ixCol];
	}

	/*cout << "1new: ";
	for (int ix = 0; ix < alIxs.size(); ix++) cout << newColumn[ix];
	cout << endl;
	cout << "old: ";
	for (int ix = 0; ix < alIxs.size(); ix++) cout << alIxs[ix][ixCol];
	cout << endl;*/
}

int GetStaggeredPositionIndex(vector<vector<int>> &alIxs, int ixColumn)
{
	int ix;
	for (ix = 0; ix < alIxs.size(); ix++)
		if (alIxs[ix][ixColumn] > -1) break;
	assert(ix < alIxs.size());
	return ix;
}

//! This function serializes the multiple part (alignment) of the superposition result.
/*!
\param NNInAvg	Stores a vector for each input structure containing for each nucleotide its nearest neighbor in the average structure which is also under given threshold (-1 if none exists)
\param ixMostSimilarInputStructure Index of the input structure which is most similar to the average structure. Sequence of the structure will be used as the one driving the multiple sequence/structure alignment.
\param threshold Threshold specifying how close residues must be to be considered conserved/aligned.
\return A string containing a multiple alignment.
*/
string LogMulti(cRNAStructure *aver, vector<cRNAStructure*> input, vector<vector<int>> NNInAvg, int ixMostSimilarInputStructure, double threshold)
{
	string positionWidth = "5";
	string emptyPosition = "      ";
	vector<vector<tPdbAtom>> atomsInput;
	for (int ix = 0; ix < input.size(); ix++) atomsInput.push_back(input[ix]->GetAllNtResidues());

	vector<vector<int>> positionsTaken; //for each input structure, it stores a vector of ixs of nucleotides which have already been used in the alignment
	//vector<string> al; //stores a string for each line of the resulting alignment	
	vector<vector<int>> alIxs; //stores the alignment as a 2D array, with rows being indexes of the aligned residues
	for (int ix = 0; ix < input.size(); ix++) 
	{
		//al.push_back("");
		alIxs.push_back(vector<int>());
		positionsTaken.push_back(vector<int>());
	}
	for (int ix = 0; ix < NNInAvg[ixMostSimilarInputStructure].size(); ix++)
	{ //iterate over all positions of the most similar structure
		if (NNInAvg[ixMostSimilarInputStructure][ix] > -1)
		{ //if that positions was within given threshold from the average structure
			string resNum = emptyPosition;
			resNum = (boost::format("%|+" + positionWidth + "| ") % (boost::trim_copy(atomsInput[ixMostSimilarInputStructure][ix].residue_num) + atomsInput[ixMostSimilarInputStructure][ix].ins_code)).str();
			//al[ixMostSimilarInputStructure] += resNum;
			if (resNum == emptyPosition) alIxs[ixMostSimilarInputStructure].push_back(-1); else alIxs[ixMostSimilarInputStructure].push_back(ix);
			for (int ixSeq = 0; ixSeq < input.size(); ixSeq++)
				if (ixSeq != ixMostSimilarInputStructure)
				{ //iterate over all input structures different from the most similar structure
					resNum = emptyPosition;
					int ixPos = 0;
					for (; ixPos < NNInAvg[ixSeq].size(); ixPos++)
					{
						if (NNInAvg[ixSeq][ixPos] == NNInAvg[ixMostSimilarInputStructure][ix] && 
							find(positionsTaken[ixSeq].begin(), positionsTaken[ixSeq].end(), ixPos) == positionsTaken[ixSeq].end())
						{ //align, only if there exists a position of the examined structure which is close enough to the positions of the average structure
							//which is being examined from the point of view of the nearest structure ( NNInAvg[ixMostSimilarInputStructure][ix] )
							resNum = (boost::format("%|+" + positionWidth + "| ") % (boost::trim_copy(atomsInput[ixSeq][ixPos].residue_num) + atomsInput[ixSeq][ixPos].ins_code)).str();
							positionsTaken[ixSeq].push_back(ixPos);
							break;
						}						
					}
					//al[ixSeq] += resNum;
					if (resNum == emptyPosition) alIxs[ixSeq].push_back(-1); else alIxs[ixSeq].push_back(ixPos);
				}			
		}
	}

	//impute the residues which are not aligned with the average structure
	for (int ix = 0; ix < input.size(); ix++)
	{
		int ixPos = 0;
		int ixLastNonEmptyPosValue = -1;
		while (ixPos < alIxs[ix].size())
		{
			if (alIxs[ix][ixPos] > -1)
			{	
				if (ixPos == 0 || alIxs[ix][ixPos-1] == -1 || alIxs[ix][ixPos] > alIxs[ix][ixPos-1] + 1)
				{//if we encounter a non-empty position right after an empty position we
					//insert the non-aligned residues which are missing
					for (int ixToInsert = alIxs[ix][ixPos] - 1; ixToInsert > ixLastNonEmptyPosValue; ixToInsert--)
					{
						for (int ixAux = 0; ixAux < input.size(); ixAux++)
						{//we insert the index into the ix sequence and -1 to the rest of the sequences for each of the inserted residues						
							alIxs[ixAux].insert(alIxs[ixAux].begin() + ixPos, -1);
						}
						alIxs[ix][ixPos] = ixToInsert;
					}
				}
				ixLastNonEmptyPosValue = alIxs[ix][ixPos] ;
			}
			ixPos++;
		}
	}

	//Identify staggered (imputed) positions and possibly merge them with neighboring staggered positions
	//First identify ranges of neighboring staggered positions (ranges)
	vector<pair<int, int>> staggeredRanges;
	bool inRange = false;
	for (int ixCol = 0; ixCol < alIxs[0].size(); ixCol++)
	{//for each position
		if (IsStaggeredPosition(alIxs, ixCol))
		{
			if (!inRange) staggeredRanges.push_back(std::pair<int, int>(ixCol, ixCol));
			else staggeredRanges[staggeredRanges.size()-1].second++;
			inRange = true;
		}
		else
			inRange = false;
	}
	/*cout << "staggered: ";
	for (int i = 0; i < staggeredRanges.size(); i++) cout << staggeredRanges[i].first << "-" << staggeredRanges[i].second << " ";*/
	//Go over the ranges backwards
	for (int ix = staggeredRanges.size()-1; ix >= 0; ix--)
	{
		//for each range go position by position and extract groups of residues 
		//in which each pair is within given distance threshold
		vector<pair<int, t3DCoords> > rangeIndexes; //index in the alignment with the coordinates of the respective (staggered) residue
		//extract indexes in the range together with the coordinates
		for (int ixAux = staggeredRanges[ix].first; ixAux <= staggeredRanges[ix].second; ixAux++) rangeIndexes.push_back(std::pair<int, t3DCoords>(ixAux, GetCoordiantesFromAlignment(input, alIxs, ixAux)));
		int ixInRange = 0;		
		vector<vector<pair<int, t3DCoords>>> groups;
		while (rangeIndexes.size() > 0)
		{
			vector<pair<int, t3DCoords>> group;
			//push the last non-assigned index in the range into a new grop
			//group.push_back(rangeIndexes[rangeIndexes.size()-1]);
			//remove the index from the range
			//rangeIndexes.erase(rangeIndexes.end()-1);
			//for the remaining non-assigned residues check the condition and possibly assign them into the new group
			vector<int> alreadyUsedStructures;
			for (int ixRange = rangeIndexes.size()-1; ixRange >= 0; ixRange--)
			{
				//first we have to checke whether the sequence whose residue is staggered has not already been assigned to the group
				//that could happen in situations like:
				// G G                 G
				//     A A U A         G
				//   G         A A U A G
				int staggeredPositionIndex = GetStaggeredPositionIndex(alIxs, rangeIndexes[ixRange].first);
				//cout << "st: " << rangeIndexes[ixRange].first << ": " << staggeredPositionIndex << endl;
				if (find(alreadyUsedStructures.begin(), alreadyUsedStructures.end(), staggeredPositionIndex) != alreadyUsedStructures.end()) continue;
				if (IsAssignableIntoGroup(rangeIndexes[ixRange].second, group, threshold))
				{
					group.push_back(rangeIndexes[ixRange]);
					rangeIndexes.erase(rangeIndexes.begin()+ixRange);
					alreadyUsedStructures.push_back(staggeredPositionIndex);
				}
			}
			groups.push_back(group);
		}
		//when the staggered region is processed, each of the groups will be merged into a single column
		vector<vector<int>> newColumns;
		for (int ixGroup = groups.size()-1; ixGroup >= 0; ixGroup--)
		{
			vector<int> newColumn;
			for (int ixAux = 0; ixAux < input.size(); ixAux++) newColumn.push_back(-1);
			vector<pair<int, t3DCoords>> group = groups[ixGroup];
			for (int ixPosition = 0; ixPosition < group.size(); ixPosition++) MergeColumns(alIxs, newColumn , group[ixPosition].first);
			newColumns.push_back(newColumn);
		}
		//finally, remove the staggered columns and replace them with the merged groups
		for (int ixCol = staggeredRanges[ix].second; ixCol >= staggeredRanges[ix].first; ixCol--)
			for (int ixRow = 0; ixRow < alIxs.size(); ixRow++)
				alIxs[ixRow].erase(alIxs[ixRow].begin() + ixCol);
		for (int ixNewColumns = newColumns.size()-1; ixNewColumns >=0 ; ixNewColumns--)
			for (int ixRow = 0; ixRow < alIxs.size(); ixRow++)
				alIxs[ixRow].insert(alIxs[ixRow].begin() + staggeredRanges[ix].first, newColumns[ixNewColumns][ixRow]);
		
	}

	string alignment = "";
	//for (int ix = 0; ix < input.size(); ix++) alignment += al[ix] + "\n";

	//convert the information about aligned indexes into the alignment string
	for (int ix = 0; ix < alIxs.size(); ix++)
	{//for each structure
		for (int ixPos = 0; ixPos < alIxs[ix].size(); ixPos++)
		{//for each position
			//if the position contains -1, output empty position else output given position
			if (alIxs[ix][ixPos] == -1) alignment += emptyPosition;
			//else alignment += (boost::format("%|+" + positionWidth + "| ") % (boost::trim_copy(atomsInput[ix][alIxs[ix][ixPos]].residue_num) + atomsInput[ix][alIxs[ix][ixPos]].ins_code)).str();
			else alignment += (boost::format("%|+" + positionWidth + "| ") % (boost::trim_copy(atomsInput[ix][alIxs[ix][ixPos]].res_name))).str();
		}
		alignment += "\n";
	}

	return alignment;
}

//! This function serializes the pairwise superposition result into a XML log file .
/*!
	\param logFileName	name of the log file. If the string is empty, the serialized file content will be passed using the return value.
	\param res reference to the pairwise superposition result structure
	\param rna array with two RNA structures
	\param id1 Identifier of the first structure.
	\param id2 Identifier of the first structure.
	\param NNInAvg This aims to be used in multiple superposition as a 2D array of nearest neighboring nt in AVG structure to each of the input structures' nt. 
		Neighbor is a nt which is not farther than the threshold. If there are multiple nts in given distance the nearest is selected.
		The structure is built iteratively as LogResult is being called for every of the structures.

	\return Either an empty string (if the content is serialized into the logFileName) or a XML string containing seriazlied superposition result.
	*/
string LogResult(string logFileName, ResMatch res, cRNAStructure* rna[2], string id1 = "", string id2 = "", string chain1 = "", string chain2 = "", 
	vector<vector<int>> *NNInAvg = new vector<vector<int>>(), double threshold = 0)
{
	using namespace rapidxml;

	xml_document<> outDoc;
	xml_node<> *nodeRoot = outDoc.allocate_node(node_element, "setter");
	xml_node<> *nodeMeta = outDoc.allocate_node(node_element, "meta");
	xml_node<> *nodeResult = outDoc.allocate_node(node_element, "result");
	xml_node<> *nodeAux;
	xml_attribute<> *attrAux;
	char *node_sub_result;

	//meta information serialization
	nodeAux = outDoc.allocate_node(node_element, "time");
	char dateStart[50], dateEnd[50], timeStart[50], timeEnd[50], span[50];
	strcpy(dateStart, boost::gregorian::to_iso_extended_string(res.timeStart.date()).c_str());
	strcpy(dateEnd, boost::gregorian::to_iso_extended_string(res.timeEnd.date()).c_str());
	strcpy(timeStart, boost::posix_time::to_simple_string(res.timeStart.time_of_day()).c_str());
	strcpy(timeEnd, boost::posix_time::to_simple_string(res.timeEnd.time_of_day()).c_str());
	strcpy(span, boost::posix_time::to_simple_string(res.timeEnd - res.timeStart).c_str());
	nodeAux->append_attribute(outDoc.allocate_attribute("start_date", dateStart));
	nodeAux->append_attribute(outDoc.allocate_attribute("end_date", dateEnd));
	nodeAux->append_attribute(outDoc.allocate_attribute("start_time", timeStart));
	nodeAux->append_attribute(outDoc.allocate_attribute("end_time", timeEnd));
	nodeAux->append_attribute(outDoc.allocate_attribute("run_time_ms", span));
	nodeMeta->append_node(nodeAux);

	ostringstream strs; //result information serialization
	char *node_str;

	strs.str(""); strs.clear();	strs << id1; node_str = outDoc.allocate_string(strs.str().c_str());
	nodeAux->append_attribute(outDoc.allocate_attribute("id1", node_str));
	strs.str(""); strs.clear();	strs << id2; node_str = outDoc.allocate_string(strs.str().c_str());
	nodeAux->append_attribute(outDoc.allocate_attribute("id2", node_str));

	nodeAux = outDoc.allocate_node(node_element, "size");
	strs.str(""); strs.clear();	strs << rna[0]->GetLengthNt(); node_str = outDoc.allocate_string(strs.str().c_str());
	nodeAux->append_attribute(outDoc.allocate_attribute("nt_cnt1", node_str));
	strs.str(""); strs.clear();	strs << rna[1]->GetLengthNt(); node_str = outDoc.allocate_string(strs.str().c_str());
	nodeAux->append_attribute(outDoc.allocate_attribute("nt_cnt2", node_str));
	strs.str(""); strs.clear();	strs << rna[0]->GetHairpinsCount(); node_str = outDoc.allocate_string(strs.str().c_str());
	nodeAux->append_attribute(outDoc.allocate_attribute("gssu_cnt1", node_str));
	strs.str(""); strs.clear();	strs << rna[1]->GetHairpinsCount(); node_str = outDoc.allocate_string(strs.str().c_str());
	nodeAux->append_attribute(outDoc.allocate_attribute("gssu_cnt2", node_str));
	nodeMeta->append_node(nodeAux);

	//GSSUs definition
	xml_node<> *nodeGSSUs = outDoc.allocate_node(node_element, "GSSUs");
	char gssus_name[2][50];
	for (int k = 0; k < 2; k++)
	{
		sprintf(gssus_name[k], "GSSUs%d", k + 1);
		xml_node<> *nodeGSSU = outDoc.allocate_node(node_element, gssus_name[k]);
		int cntHp = rna[k]->GetHairpinsCount();
		for (int ixHP = 0; ixHP < cntHp; ixHP++)
		{
			char *strIx = new char[10];
			sprintf(strIx, "GSSU%d", ixHP + 1);
			xml_node<> *nodeGSSUi = outDoc.allocate_node(node_element, strIx);
			cHairpin hp = rna[k]->GetHairpin(ixHP);
			//serialization of the first part of the stem
			strs.str(""); strs.clear();
			for (int ixAux = hp.stem.pairs.size() - 1; ixAux >= 0; ixAux--)
			{
				if (hp.stem.pairs[ixAux].r1.ix_residue != -1)
				{
					string auxChain = hp.stem.pairs[ixAux].r1.chain;
					if (k == 0 && chain1 != "") auxChain = chain1; else if (k == 1 && chain2 != "") auxChain = chain2;
					strs << auxChain << "-" << hp.stem.pairs[ixAux].r1.residue_position << ";";
				}
			}
			string strAux = strs.str();
			if (strAux.find_last_of(";") == strAux.size() - 1) strAux = strAux.substr(0, strAux.size() - 1);
			node_sub_result = outDoc.allocate_string(strAux.c_str());
			nodeGSSUi->append_node(outDoc.allocate_node(node_element, "stem1", node_sub_result));
			//serialization of the head
			strs.str(""); strs.clear();
			for (int ixAux = 0; ixAux < hp.head.residues.size(); ixAux++)
			{
				if (ixAux > 0) strs << ";";
				string auxChain = hp.head.residues[ixAux].r1.chain;
				if (k == 0 && chain1 != "") auxChain = chain1; else if (k == 1 && chain2 != "") auxChain = chain2;
				strs << auxChain << "-" << hp.head.residues[ixAux].r1.residue_position;
			}
			node_sub_result = outDoc.allocate_string(strs.str().c_str());
			nodeGSSUi->append_node(outDoc.allocate_node(node_element, "head", node_sub_result));
			//serialization of the second part of the stem
			strs.str(""); strs.clear();
			bool first = true;
			for (int ixAux = 0; ixAux < hp.stem.pairs.size(); ixAux++)
			{
				if (hp.stem.pairs[ixAux].r2.ix_residue != -1)
				{
					if (!first) strs << ";";
					string auxChain = hp.stem.pairs[ixAux].r2.chain;
					if (k == 0 && chain1 != "") auxChain = chain1; else if (k == 1 && chain2 != "") auxChain = chain2;
					strs << auxChain << "-" << hp.stem.pairs[ixAux].r2.residue_position;
					first = false;
				}
			}
			strAux = strs.str();
			if (strAux.find_last_of(";") == strAux.size() - 1) strAux = strAux.substr(0, strAux.size() - 1);
			node_sub_result = outDoc.allocate_string(strAux.c_str());
			nodeGSSUi->append_node(outDoc.allocate_node(node_element, "stem2", node_sub_result));
			nodeGSSU->append_node(nodeGSSUi);
		}
		nodeGSSUs->append_node(nodeGSSU);
	}
	nodeMeta->append_node(nodeGSSUs);

	nodeRoot->append_node(nodeMeta);

	//score
	strs.str(""); strs.clear();
	strs << res.score;
	node_sub_result = outDoc.allocate_string(strs.str().c_str());
	nodeResult->append_node(outDoc.allocate_node(node_element, "score", node_sub_result));

	//p-value
	double muAux = 0, sigmaAux = 0;
	int minLengthAux = min(rna[0]->GetLengthNt(), rna[1]->GetLengthNt());
	muAux = 0.19397245346996325 * pow(minLengthAux, 0.5) + 25.466345138393969 * 1 / minLengthAux;
	sigmaAux = 8.7053141611673599 * pow(minLengthAux, -0.5) + 0.024897807792426670 * pow(log((double)minLengthAux), 2);
	double pValue = 0.5 + 0.5 * erf((log(res.score) - muAux) / sqrt(2 * sigmaAux * sigmaAux));
	strs.str(""); strs.clear();
	strs << pValue;
	node_sub_result = outDoc.allocate_string(strs.str().c_str());
	nodeResult->append_node(outDoc.allocate_node(node_element, "p-value", node_sub_result));

	//rotation matrix
	strs.str(""); strs.clear();
	strs << "["
		<< "[" << res.rot[0][0] << "," << res.rot[0][1] << "," << res.rot[0][2] << "],"
		<< "[" << res.rot[1][0] << "," << res.rot[1][1] << "," << res.rot[1][2] << "],"
		<< "[" << res.rot[2][0] << "," << res.rot[2][1] << "," << res.rot[2][2] << "]]";
	node_sub_result = outDoc.allocate_string(strs.str().c_str());
	nodeResult->append_node(outDoc.allocate_node(node_element, "rotation", node_sub_result));

	//translation vector
	strs.str(""); strs.clear();
	strs << "[" << res.trans[0] << "," << res.trans[1] << "," << res.trans[2] << "]";
	node_sub_result = outDoc.allocate_string(strs.str().c_str());
	nodeResult->append_node(outDoc.allocate_node(node_element, "translation", node_sub_result));

	//NN set - for each nt in one RNA, it's nearest neighbor is showed	
	char nnset_name[2][10];
	for (int k = 0; k < 2; k++)
	{
		strs.str(""); strs.clear();
		for (int i = 0; i < res.nn[k].size(); i++)
		{
			if (i > 0) strs << ";";
			strs << res.nn[k][i].first.first << ";" << res.nn[k][i].first.second << ";" << res.nn[k][i].second;
		}
		node_sub_result = outDoc.allocate_string(strs.str().c_str());

		sprintf(nnset_name[k], "NNSet%d", k + 1);
		nodeResult->append_node(outDoc.allocate_node(node_element, nnset_name[k], node_sub_result));
	}

	//NN pairing - pairs of mutual NNs (if NN of a nt is also NN of the NN)
	vector<pair<int, int> > nnMutual;
	vector<tPdbAtom> atoms[2];
	vector<vector<double> > nnDists;	
	for (int i = 0; i < 2; i++) atoms[i] = rna[i]->GetAllNtResidues();
	for (int i = 0; i < atoms[1].size(); i++) TranslateCoord(atoms[1][i].coords, res.rot, res.trans);

	vector<vector<int> > nnAux;
	//initialization
	for (int i = 0; i < atoms[0].size(); i++)
	{
		vector<int> nnAux1;
		vector<double> nnDistsAux;
		for (int j = 0; j < atoms[1].size(); j++)
		{
			nnAux1.push_back(0);
			nnDistsAux.push_back(0);
		}

		nnAux.push_back(nnAux1);
		nnDists.push_back(nnDistsAux);
	}
	//NNs for nts in the first structure
	for (int i = 0; i < atoms[0].size(); i++)
	{
		double distMin = 1000000;
		int ixMin = -1;
		for (int j = 0; j < atoms[1].size(); j++)
		{
			double dist = sqrt(pow(atoms[0][i].coords.x - atoms[1][j].coords.x, 2) + pow(atoms[0][i].coords.y - atoms[1][j].coords.y, 2) + pow(atoms[0][i].coords.z - atoms[1][j].coords.z, 2));
			nnDists[i][j] = dist;
			if (dist < distMin)
			{
				distMin = dist;
				ixMin = j;
			}
		}		
		nnAux[i][ixMin]++;
	}
	//NNs for nts in the second structure
	//In the case of MultiSETTER the first structure is the average one
	if (threshold > 0) NNInAvg->push_back(vector<int>());
	for (int i = 0; i < atoms[1].size(); i++)
	{
		double distMin = 100000;
		int ixMin = -1;
		for (int j = 0; j < atoms[0].size(); j++)
		{
			double dist = sqrt(pow(atoms[1][i].coords.x - atoms[0][j].coords.x, 2) + pow(atoms[1][i].coords.y - atoms[0][j].coords.y, 2) + pow(atoms[1][i].coords.z - atoms[0][j].coords.z, 2));
			if (dist < distMin)
			{
				distMin = dist;
				ixMin = j;
			}
		}
		nnAux[ixMin][i]++;
		if (threshold > 0)
		{
			int ixAux = -1;
			if (distMin < threshold) ixAux = ixMin;
			(*NNInAvg)[NNInAvg->size() - 1].push_back(ixAux);
		}
	}
	for (int i = 0; i < atoms[0].size(); i++)
		for (int j = 0; j < atoms[1].size(); j++) if (nnAux[i][j] == 2) nnMutual.push_back(pair<int, int>(i, j));
	strs.str(""); strs.clear();
	const double distThreshold = 4;
	int cntInDistance = 0;
	int cntInDistanceSameType = 0;
	double rmsd = 0;
	for (int i = 0; i < nnMutual.size(); i++)
	{
		if (i > 0) strs << ";";
		string auxChain1 = chain1 == "" ? string(1, atoms[0][nnMutual[i].first].chain) : chain1;
		string auxChain2 = chain2 == "" ? string(1, atoms[1][nnMutual[i].second].chain) : chain2;
		strs << auxChain1 << "-" << atoms[0][nnMutual[i].first].residue_num << "_" << atoms[0][nnMutual[i].first].ins_code << ";"
			<< auxChain2 << "-" << atoms[1][nnMutual[i].second].residue_num << "_" << atoms[1][nnMutual[i].second].ins_code;

		rmsd += atoms[0][nnMutual[i].first].coords.DistL2(atoms[1][nnMutual[i].second].coords);

		if (nnDists[nnMutual[i].first][nnMutual[i].second] <= distThreshold)
		{
			cntInDistance++;
			if (boost::trim_copy(atoms[0][nnMutual[i].first].res_name) == boost::trim_copy(atoms[1][nnMutual[i].second].res_name)) cntInDistanceSameType++;
		}
	}
	rmsd /= (double)nnMutual.size();
	rmsd = sqrt(rmsd);
	node_sub_result = outDoc.allocate_string(strs.str().c_str());
	nodeResult->append_node(outDoc.allocate_node(node_element, "mutualNNSet", node_sub_result));

	//rmsd
	strs.str(""); strs.clear();
	strs << rmsd;
	node_sub_result = outDoc.allocate_string(strs.str().c_str());
	nodeResult->append_node(outDoc.allocate_node(node_element, "RMSD", node_sub_result));

	//number of aligned nucleotides in 4A
	strs.str(""); strs.clear();
	strs << cntInDistance;
	node_sub_result = outDoc.allocate_string(strs.str().c_str());
	nodeResult->append_node(outDoc.allocate_node(node_element, "cnt_aligned_nucleotides", node_sub_result));

	//number of exact base matches in 4A
	strs.str(""); strs.clear();
	strs << cntInDistanceSameType;
	node_sub_result = outDoc.allocate_string(strs.str().c_str());
	nodeResult->append_node(outDoc.allocate_node(node_element, "cnt_exact_based_matches", node_sub_result));

	//PSI
	strs.str(""); strs.clear();
	strs << cntInDistance / (double)min(atoms[0].size(), atoms[1].size());
	node_sub_result = outDoc.allocate_string(strs.str().c_str());
	nodeResult->append_node(outDoc.allocate_node(node_element, "PSI", node_sub_result));

	//PID
	strs.str(""); strs.clear();
	strs << (double)cntInDistanceSameType / min(atoms[0].size(), atoms[1].size());
	node_sub_result = outDoc.allocate_string(strs.str().c_str());
	nodeResult->append_node(outDoc.allocate_node(node_element, "PID", node_sub_result));

	//P atom positions
	char reprpos_name[2][50];
	for (int k = 0; k < 2; k++)
	{
		strs.str(""); strs.clear();
		for (int i = 0; i < atoms[k].size(); i++)
		{
			string auxChain = string(1, atoms[k][i].chain);
			if (k == 0 && chain1 != "") auxChain = chain1; else if (k == 1 && chain2 != "") auxChain = chain2;
			strs << auxChain << "-" << atoms[k][i].residue_num << "-" << boost::trim_copy(atoms[k][i].res_name) << "[" << atoms[k][i].coords.x << ";" << atoms[k][i].coords.y << ";" << atoms[k][i].coords.z << "]";
		}
		node_sub_result = outDoc.allocate_string(strs.str().c_str());

		sprintf(reprpos_name[k], "superposed_representative_coordinates%d", k + 1);
		nodeResult->append_node(outDoc.allocate_node(node_element, reprpos_name[k], node_sub_result));
	}

	//aligned GSSUs
	strs.str(""); strs.clear();
	for (int i = 0; i < res.alignedGSSUs.size(); i++)
	{
		if (i > 0) strs << ";";
		strs << res.alignedGSSUs[i].first + 1 << "-" << res.alignedGSSUs[i].second + 1;
		//<< "(" << res.alignedScores[i] << ")";		
	}

	//strs << " - rotation-penalties: ";
	//for (int i = 0; i < res.rotationPenalties.size(); i++) strs << res.rotationPenalties[i] << "-";

	node_sub_result = outDoc.allocate_string(strs.str().c_str());
	nodeResult->append_node(outDoc.allocate_node(node_element, "aligned_GSSUs", node_sub_result));

	//principal alignment triplet
	strs.str(""); strs.clear();
	for (int i = 0; i < 3; i++)
	{
		try{
			tPdbAtom at[2] = { rna[0]->GetAtom(res.superpositionTriplet[i].r1), rna[1]->GetAtom(res.superpositionTriplet[i].r2) };
			if (i > 0) strs << ";";
			string auxChain1 = chain1 == "" ? string(1, at[0].chain) : chain1;
			string auxChain2 = chain2 == "" ? string(1, at[1].chain) : chain2;
			strs << auxChain1 << "-" << at[0].residue_num << ":" << auxChain2 << "-" << at[1].residue_num;
		}
		catch (const char* err)
		{
			strs << "none";
		}
	}
	node_sub_result = outDoc.allocate_string(strs.str().c_str());
	nodeResult->append_node(outDoc.allocate_node(node_element, "superposition_triplet", node_sub_result));

	nodeRoot->append_node(nodeResult);
	outDoc.append_node(nodeRoot);

	if (boost::trim_copy(logFileName) != "")
	{
		ofstream ofs(logFileName.c_str());
		if (!ofs.is_open()) cerr << "Problem with opening the log file " << logFileName << "." << endl;
		else ofs << outDoc;
		return "";
	}
	else
	{
		std::string s;
		rapidxml::print(std::back_inserter(s), outDoc, 0);
		return s;
	}
}

//! Parses the input xml and outputs the parameters into the BOOST variables map
void ParseInputXML(string xmlFileName, po::variables_map &vm, vector<InputRNA> &rnaInputs)
{
	ifstream ifs(xmlFileName);
	if (ifs.is_open())
	{
		vm.clear();

		using namespace rapidxml;

		string xmlAux, line;
		while (getline(ifs, line)) xmlAux += line;
		char *xml;
		xml = new char[xmlAux.size() + 1];
		strncpy(xml, xmlAux.c_str(), xmlAux.size());
		xml[xmlAux.size()] = 0;

		xml_document<> inDoc;
		inDoc.parse<0>(xml);

		xml_node<> *ndTask = inDoc.first_node()->first_node("task");
		xml_node<> *ndArguments = inDoc.first_node()->first_node("arguments");
		xml_node<> *ndParameters = inDoc.first_node()->first_node("setter-parameters");
		xml_node<> *ndArgument, *ndParameter;

		//reading the type of task to be performed
		string type = string(ndTask->first_attribute("type")->value());
		string prefix;
		if (type == "pairwise") prefix = "p-";
		else if (type == "average") prefix = "a-";
		else if (type == "align") prefix = "l-";
		else if (type == "all-to-all") prefix = "t-";

		vm.insert(std::make_pair(type, boost::program_options::variable_value()));

		string name, value;

		//reading the arguments
		ndArgument = ndArguments->first_node("argument");
		while (ndArgument != NULL)
		{
			name = string(ndArgument->first_attribute("name")->value());
			value = string(ndArgument->first_attribute("value")->value());

			if (name != "configuration-file") name = prefix + name;
			vm.insert(std::make_pair(name, po::variable_value(value, false)));

			ndArgument = ndArgument->next_sibling();
		}

		//reading the params
		sParams params = GetGlobalParams();
		ndParameter = ndParameters->first_node("parameter");
		while (ndParameter != NULL)
		{
			name = string(ndParameter->first_attribute("name")->value());
			value = string(ndParameter->first_attribute("value")->value());

			if (name == "p_neck_shift") params.neckMaxShift = atoi(value.c_str());
			else if (name == "p_identical_pair_nt_type_modifier") params.identicalLetterModificator = atof(value.c_str());
			else if (name == "p_identical_pair_sse_type_modifier") params.identicalPairTypeModificator = atof(value.c_str());
			else if (name == "p_pair_distance_threshold") params.NNDist = atof(value.c_str());
			else if (name == "p_top_k_gssu") params.topKForMultiGSSU = atof(value.c_str());
			else if (name == "p_no_head_step_ratio_divisor") params.noHeadShiftStepRatioDiviser = atof(value.c_str());
			else if (name == "p_early_termination") params.earlyTermination = atoi(value.c_str());
			else if (name == "p_rotation_penalty") params.rotationPenalty = atof(value.c_str());

			ndParameter = ndParameter->next_sibling();
		}
		SetGlobalParams(params);

		if (type == "pairwise" || type == "average" || type == "align")
		{
			InputRNA rnaDesc;
			xml_node<> *ndStructures = inDoc.first_node()->first_node(type.c_str())->first_node("structures");
			xml_node<> *ndStructure;
			ndStructure = ndStructures->first_node("structure");
			while (ndStructure != NULL)
			{
				rnaDesc.chain = string(ndStructure->first_attribute("chain")->value());
				rnaDesc.exclusions = boost::to_upper_copy(string(ndStructure->first_attribute("exclusions")->value()));
				rnaDesc.id = string(ndStructure->first_attribute("id")->value());
				rnaDesc.pdb = string(ndStructure->first_attribute("pdbfile")->value());
				rnaDesc.ranges = boost::to_upper_copy(string(ndStructure->first_attribute("ranges")->value()));
				rnaDesc.x3dna = string(ndStructure->first_attribute("ssfile")->value());
				rnaInputs.push_back(rnaDesc);

				ndStructure = ndStructure->next_sibling();
			}
		}
		if (type == "align")
		{
			xml_node<> *ndStructure = inDoc.first_node()->first_node("align")->first_node("average-structure");
			string avgStructureFileName = string(ndStructure->first_attribute("file")->value());
			vm.insert(std::make_pair(prefix + "input-file", po::variable_value(avgStructureFileName, false)));
		}
		if (type == "align" || type == "average")
		{
			xml_node<> *ndStructure = inDoc.first_node()->first_node(type.c_str())->first_node("align-threshold");
			double threshold = atof(ndStructure->first_attribute("value")->value());
			vm.insert(std::make_pair(prefix + "align-threshold", po::variable_value(threshold, false)));
		}

		delete [] xml;

		po::notify(vm);
	}
	else throw "Unable to open the XML input file.";
}

//! The entry point of the application.
int main(int argc, char *argv[])
{
	SetGlobalParams(GetDefaultParams());

	tbb::task_scheduler_init scheduler;
	scheduler.terminate();
	scheduler.initialize(-1);


	po::options_description odTest("Testing options");
	odTest.add_options()
		("ping", "Prints 'ok' on the console.\n This is used by the Java GUI to test the existence of the program.")
		;

	po::options_description odGeneral("General options");
	odGeneral.add_options()
		("help,h", "Produce help message.")
		("configuration-file,c", po::value<string>()->default_value(""), "File containing the parameters of the application.")
		;

	po::options_description odXML("Task specified in a XML file");
	odXML.add_options()
		("xml,x", po::value<string>(), "XML input file.")
		;

	po::options_description odAllToAll("All to all pairwise superposition options");
	odAllToAll.add_options()
		("all-to-all,t", "Flag indicating all to all pairwise superposition task.")
		("t-input-file", po::value<string>(), "Input file.")
		;

	po::options_description odAverage("Options for building 1) an average RNA from a set of RNA structures 2) superposition of the structures onto the average structure");
	odAverage.add_options()
		("average,a", "Flag indicating the average structure building task.")
		("a-input-file", po::value<string>(), "File containing links to PDB and X3DNA files of the structures to be superposed.")
		("a-output-path", po::value<string>()->default_value(""), "Path for storing results. Note that serialized average structure will be created and saved allways but the PDB files and Jmol script only if 'Script' option in the configuration file is turned on.")
		("a-output-file-prefix", po::value<string>(), "Prefix of the output files.")
		("a-align-threshold", po::value<double>(), "Threshold value which specifies how close a nucleotide needs to be to the average structure to be included in the multiple alignment.")
		("a-log-file", po::value<string>()->default_value(""), "Name of the XML log file. If none is supplied, no log will be created.")
		;

	po::options_description odAlign("Options for aligning a set of RNA structures to an average structure");
	odAlign.add_options()
		("align,l", "Flag indicating the alignment task.")
		("l-input-file", po::value<string>(), "File with a serialized average structure.")
		("l-dataset-file", po::value<string>(), "File with a dataset.")
		("l-output-path", po::value<string>()->default_value(""), "Path for storing results.")
		("l-align-threshold", po::value<double>(), "Threshold value which specifies how close a nucleotide needs to be to the average structure to be included in the multiple alignment.")
		("l-log-file", po::value<string>()->default_value(""), "Name of the XML log file. If none is supplied, no log will be created.")
		;

	po::options_description odCompare("Options for comparison of two average structures.");
	odCompare.add_options()
		("pairwise,p", "Flag indicating pairwise RNA structure superposition task.")
		("p-input-pdb1", po::value<string>(), "First RNA structure (PDB).")
		("p-input-3dna1", po::value<string>(), "First secondary structure (3DNA).")
		("p-input-chain1", po::value<string>(), "Chains of the first structure to be superposed.")
		("p-ranges1", po::value<string>(), "Specification of ranges of residues to be included in the form C1-R1:C2-R2;C3-R3:C4-R4... R being the resiude sequence number and C being the chain id. If non specified, no restriction takes place.")
		("p-exclusion-list1", po::value<string>(), "Exclusion list contains the list of residues to be excluded from the comparison in the form C1-R1;C2-R2;... R being the resiude sequence number and C being the chain id. If non specified, no restriction takes place.")
		("p-input-pdb2", po::value<string>(), "Second RNA structure (PDB).")
		("p-input-3dna2", po::value<string>(), "Second secondary structure (3DNA).")
		("p-input-chain2", po::value<string>(), "Chains of the second structure to be superposed.")
		("p-ranges2", po::value<string>(), "Specification of ranges of residues to be included in the form C1-R1:C2-R2;C3-R3:C4-R4... R being the resiude sequence number and C being the chain id. If non specified, no restriction takes place.")
		("p-exclusion-list2", po::value<string>(), "Exclusion list contains the list of residues to be excluded from the comparison in the form C1-R1;C2-R2;... R being the resiude sequence number and C being the chain id. If non specified, no restriction takes place.")
		("p-output-path", po::value<string>()->default_value(""), "Path for storing the results, i.e., the output PDB files and the resulting JMol script.")
		("p-output-file-prefix", po::value<string>(), "Prefix for the output files (pdb, jmol).")
		("p-log-file", po::value<string>()->default_value(""), "Name of the XML log file. If none is supplied, no log will be created.")
		("p-extra-superposition,s", "Flag allowing to run an extra RMSD superposition after SETTER with alignmnent defined by mutual nearest neighbors of SETTER's structural superposition.")
		;

	po::options_description all("Allowed options");
	all.add(odGeneral).add(odXML).add(odAllToAll).add(odAverage).add(odAlign).add(odCompare);

	po::variables_map vm;
	po::store(parse_command_line(argc, argv, all), vm);
	po::notify(vm);

	vector<InputRNA> rnaInputs; //to be used if the task definition is specified in a XML file
	bool xmlInput = false; //to be used if the task definition is specified in a XML file
	if (vm.count("xml"))
	{
		xmlInput = true;
		//convert the content of the XML into the vm map and continue as if the arguments
		try {
			string xmlFileLocation = vm["xml"].as<string>();
			ParseInputXML(xmlFileLocation, vm, rnaInputs);
			string xmlDir = boost::filesystem::absolute(xmlFileLocation).remove_filename().string();
			boost::filesystem::current_path(xmlDir);
		}
		catch (const char* err)
		{
			cerr << "There was a problem parsing the XML input file: " << err << endl;
			return EXIT_FAILURE;
		}
		catch (string err)
		{
			cerr << "There was a problem parsing the XML input file: " << err << endl;
			return EXIT_FAILURE;
		}
	}

	if (vm.count("help") || argc == 1) {
		cout << all << "\n";
		return EXIT_SUCCESS;
	}
	else if (vm.count("all-to-all"))
	{
		Input input(vm["t-input-file"].as<string>());

		int begin = clock();
		cout << setiosflags(ios::fixed) << setprecision(4);
		for (int i = 0; i < input.mInputSize; ++i) {
			for (int j = 0; j < input.mInputSize; ++j) {

				sParams params = GetGlobalParams();
				double score = Match(*input.mInputRnaVect[i], *input.mInputRnaVect[j],
					params).score;

				cout << score << ";";
			}

			cout << endl;
		}
		int end = clock();

		cout << resetiosflags(ios::fixed);
		cout << "Comparison time: " << (float)(end - begin) / CLOCKS_PER_SEC << "s" << endl;
	}
	else if (vm.count("average"))
	{		
		string inputFileName, prefix, outputPath, pathAndPrefix, logFileName;
		if (!xmlInput) inputFileName = vm["a-input-file"].as<string>();
		prefix = vm["a-output-file-prefix"].as<string>();
		outputPath = vm["a-output-path"].as<string>();
		if (outputPath.size() > 0) outputPath += "/";
		pathAndPrefix = outputPath + prefix;
		string outputFileName = pathAndPrefix + ".aver";

		logFileName = vm["a-log-file"].as<string>();
		if (logFileName != "") logFileName = outputPath + logFileName;

		Parameters::Init(vm["configuration-file"].as<string>());
		Parameters &param = Parameters::GetInstance();
		scheduler.terminate();
		scheduler.initialize(param.GetThreads());

		Input family;
		try{
			cout << "Parsing the input structures..." << endl;
			if (xmlInput) family = Input(rnaInputs);
			else family = Input(inputFileName);
			cout << "Finished parsing the input structures." << endl;
		}
		catch (const char* err)
		{
			cerr << "There was a problem parsing: " << err << endl;
			return EXIT_FAILURE;
		}
		catch (string err)
		{
			cerr << "There was a problem parsing: " << err << endl;
			return EXIT_FAILURE;
		}
		if (family.mInputSize == 0) {
			cerr << "Error: size of the input dataset cannot be 0" << endl;
			return EXIT_FAILURE;
		}

		int begin = clock();
		cout << "Structures parsed. Starting..." << endl;
		MultiAlign *nj = MultiAlign::Factory(param.GetAlgorithm(), family);
		MAResult aver = nj->GetResult();
		int end = clock();
		cout << "Finished. Comparison time: " << (float)(end - begin) / CLOCKS_PER_SEC << "s" << endl;

		ofstream ofs(outputFileName);
		if (!ofs.is_open())
		{
			cerr << "There was a problem storing the average structure into " << outputFileName << ". Does the output directory exist?" << endl;
			return EXIT_FAILURE;
		}
		o_archive oa(ofs);
		oa << aver;

		// Jmol script generator
		ofstream ofsScript(pathAndPrefix + ".jmol");		
		ofsScript << "background White" << endl;
		ofsScript << "appendNew = true" << endl;

		Printer::PrintRNA(aver.rna, pathAndPrefix + ".pdb", 'A');

		ofsScript << "load " << prefix << ".pdb" << endl;

		for (int i = 0; i < family.mInputSize; ++i) {
			ResMatch &res = aver.results[i];
			string prefixName = prefix + "_" + family.mOriginalMap[i] + ".pdb";
			string name = outputPath + prefixName;
			char chainID = 'B' + i;
			Printer::PrintRNA(*family.mInputRnaVect[i], res, name, chainID);
			ofsScript << "load append " << prefixName << endl;
		}

		ofsScript << "select all; trace only;" << endl;

		vector<string> colors = GetRGBColors(family.mInputSize);
		for (int i = 0; i < family.mInputSize; ++i) {
			ofsScript << "select */" << i + 2 << "; color trace " << colors[i] << endl;
		}

		ofsScript << "select */1; color trace Black" << endl;
		ofsScript << endl;
		ofsScript.close();

		/*
		* MultiSETTER logs the same way as SETTER but instead of one pair, 
		* it outputs superposition of all the input structures
		* over the average structure.
		*/
		if (logFileName != "")
		{			
			ofstream ofs(logFileName.c_str());
			if (!ofs.is_open()) cerr << "There was a problem opening the log file " << logFileName << "." << endl;
			else
			{
				vector<vector<int>> NNInAvg;
				double threshold = vm["a-align-threshold"].as<double>();

				cRNAStructure *rna[2];
				rna[0] = &(aver.rna); //when superposing in the ParallelMainBlock the average RNA is allways the first RNA

				ofs << "<multisetter>" << endl;
				double distMin = 1000000;
				int ixMin = 0;
				for (int i = 0; i < family.mInputSize; i++)
				{
					if (aver.dist[i] < distMin)
					{
						distMin = aver.dist[i];
						ixMin = i;
					}
					rna[1] = family.mInputRnaVect[i];
					string id = family.mOriginalMap[i];
					if (xmlInput) id = rnaInputs[i].id;
					ofs << LogResult("", aver.results[i], rna, "average", id, "A", "", &NNInAvg, threshold);
				}
				ofs << "<alignment>" << endl;
				ofs << LogMulti(&(aver.rna), family.mInputRnaVect, NNInAvg, ixMin, threshold);
				ofs << "</alignment>" << endl;
				ofs << "</multisetter>" << endl;
				ofs.close();
			}
		}
	}
	else if (vm.count("align"))
	{
		string inputFileName = vm["l-input-file"].as<string>();

		string outputPath = vm["l-output-path"].as<string>();
		if (outputPath.size() > 0) outputPath += "/";

		string logFileName = vm["l-log-file"].as<string>();
		if (logFileName != "") logFileName = outputPath + logFileName;

		Parameters::Init(vm["configuration-file"].as<string>());
		Parameters &param = Parameters::GetInstance();
		scheduler.terminate();
		scheduler.initialize(param.GetThreads());

		/*
		//Experimental evaluation purposes
		Input family(familyPath);
		if (family.mInputSize == 0) {
		cerr << "ERrror while reading the family..." << endl;
		return 1;
		}
		Input dataset(datasetFileName);
		if (dataset.mInputSize == 0) {
		cerr << "Error while reading the dataset..." << endl;
		return 1;
		}

		CompairDatasetToFamily(family, dataset);
		*/

		cout << "Reading the average structure..." << endl;
		MAResult aver;
		ifstream ifs(inputFileName);
		i_archive ia(ifs);
		ia >> aver;

		cout << "Mean     : " << aver.mean << endl;
		cout << "Deviation: " << aver.deviation << endl;

		cout << "Reading the dataset..." << endl;

		Input dataset;
		try{
			if (xmlInput) dataset = Input(rnaInputs);
			else dataset = Input(vm["l-dataset-file"].as<string>());
		}
		catch (const char* err)
		{
			cerr << "There was a problem when parsing: " << err << endl;
			return EXIT_FAILURE;
		}
		catch (string err)
		{
			cerr << "There was a problem parsing: " << err << endl;
			return EXIT_FAILURE;
		}
		if (dataset.mInputSize == 0) {
			cerr << "Error: size of the input dataset cannot be 0" << endl;
			return EXIT_FAILURE;
		}

		vector<sResMatch> output(dataset.mInputRnaVect.size());
		vector<int> positions(dataset.mInputRnaVect.size());
		for (int i = 0; i < dataset.mInputRnaVect.size(); ++i) {
			positions[i] = i;
		}

		cout << "Comparing the dataset to the average structure..." << endl;

		ParallelMainBlock parallelMainBlock(
			aver.rna, positions, dataset.mInputRnaVect, dataset.mOriginalMap, output);

		tbb::parallel_for(
			tbb::blocked_range<size_t>(0, dataset.mInputRnaVect.size()), parallelMainBlock);
		parallelMainBlock.PrintOutput();


		if (logFileName != "")
		{
			ofstream ofs(logFileName.c_str());
			if (!ofs.is_open()) cerr << "There was a problem opening the log file " << logFileName << "." << endl;
			else
			{
				cRNAStructure *rna[2];
				rna[0] = &(aver.rna); //when superposing in the ParallelMainBlock the average RNA is allways the first RNA

				vector<sResMatch> results = parallelMainBlock.results();
				for (int i = 0; i < results.size(); i++)
				{
					rna[1] = dataset.mInputRnaVect[i];
					ofs << LogResult("", results[i], rna, "", "", "A");
				}
				ofs.close();
			}
		}

	}
	else if (vm.count("pairwise"))
	{
		Parameters::Init(vm["configuration-file"].as<string>());
		Parameters &param = Parameters::GetInstance();

		scheduler.terminate();
		scheduler.initialize(param.GetThreads());		

		string rnaFilename[2], x3dnaFileName[2], chain[2], exclusionList[2] = { "", "" }, ranges[2] = { "", "" }, logFileName;

		bool extraSuperposition = vm.count("p-extra-superposition");

		if (xmlInput)
		{
			for (int i = 0; i < 2; i++)
			{
				rnaFilename[i] = rnaInputs[i].pdb;
				x3dnaFileName[i] = rnaInputs[i].x3dna;
				chain[i] = rnaInputs[i].chain;
				ranges[i] = rnaInputs[i].ranges;
				exclusionList[i] = rnaInputs[i].exclusions;
			}
		}
		else
		{
			rnaFilename[0] = vm["p-input-pdb1"].as<string>();
			x3dnaFileName[0] = vm["p-input-3dna1"].as<string>();
			chain[0] = vm["p-input-chain1"].as<string>();
			if (vm.count("p-ranges1") > 0) ranges[0] = boost::to_upper_copy(vm["p-ranges1"].as<string>());
			if (vm.count("p-exclusion-list1") > 0) exclusionList[0] = boost::to_upper_copy(vm["p-exclusion-list1"].as<string>());
			rnaFilename[1] = vm["p-input-pdb2"].as<string>();
			x3dnaFileName[1] = vm["p-input-3dna2"].as<string>();
			chain[1] = vm["p-input-chain2"].as<string>();
			if (vm.count("p-ranges2") > 0) ranges[1] = boost::to_upper_copy(vm["p-ranges2"].as<string>());
			if (vm.count("p-exclusion-list2") > 0) exclusionList[1] = boost::to_upper_copy(vm["p-exclusion-list2"].as<string>());
		}

		string prefix = vm["p-output-file-prefix"].as<string>();
		string outputPath = vm["p-output-path"].as<string>();
		if (outputPath.size() > 0) outputPath += "/";
		string pathAndPrefix = outputPath + prefix;

		logFileName = vm["p-log-file"].as<string>();
		if (logFileName != "") logFileName = outputPath + logFileName;

		cRNAStructure *rna[2];
		for (int i = 0; i < 2; i++)
			try{
			rna[i] = new cRNAStructure(rnaFilename[i], x3dnaFileName[i], chain[i], exclusionList[i], ranges[i]);
		}
		catch (const char* err)
		{
			cerr << "Problem when parsing the input structre " << rnaFilename[i] << ": " << err << endl;
			return EXIT_FAILURE;
		}
		catch (string err)
		{
			cerr << "Problem when parsing the input structre " << rnaFilename[i] << ": " << err << endl;
			return EXIT_FAILURE;
		}

		sParams params = GetGlobalParams();
		cout << "Starting matching...";
		int begin = clock();
		ResMatch res = Match(*rna[0], *rna[1], params, extraSuperposition);
		int end = clock();
		cout << "Comparison time: " << (float)(end - begin) / CLOCKS_PER_SEC << "s" << endl;
		cout << "Distance: " << res.score << endl;
		string ids[2] = { "", "" };
		if (xmlInput) {
			ids[0] = rnaInputs[0].id;
			ids[1] = rnaInputs[1].id;
		}
		if (logFileName != "") LogResult(logFileName, res, rna, ids[0], ids[1]);

		ofstream ofsScript(pathAndPrefix + ".jmol");

		char chainId = 'X';
		if (chain[0].size() > 0) chainId = chain[0][0];
		Printer::PrintRNA(*rna[0], pathAndPrefix + "_1.pdb", chainId);
		chainId = 'X';
		if (chain[1].size() > 0) chainId = chain[1][0];
		Printer::PrintRNA(*rna[1], res, pathAndPrefix + "_2.pdb", chainId);

		ofsScript << "background White" << endl;
		ofsScript << "appendNew = true" << endl;
		ofsScript << "load " << pathAndPrefix + "_1.pdb" << endl;
		ofsScript << "load append " << pathAndPrefix + "_2.pdb" << endl;

		ofsScript << "select all; trace only;" << endl;
		ofsScript << "select */1" << "; color trace " << "RED" << endl;
		ofsScript << "select */2" << "; color trace " << "BLUE" << endl;
		ofsScript << endl;
	}


	return EXIT_SUCCESS;
}
