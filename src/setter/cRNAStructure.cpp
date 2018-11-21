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
#include "cRNAStructure.h"
#include "common.h"

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <boost/serialization/string.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>

#include <boost/algorithm/string.hpp>
#include <boost/unordered_map.hpp>



PairResiude::PairResiude() {ix_residue = -1;}
PairResiude::PairResiude(const PairResiude &res) 
{
	chain = res.chain;
	residue_position = res.residue_position;		
	residue_label_short = res.residue_label_short;	
	residue_label_long = res.residue_label_long;
	
	ix_residue = res.ix_residue;
	ix_chain = res.ix_chain;
}

PairResiude &PairResiude::operator=(const PairResiude &res) 
{
	chain = res.chain;
	residue_position = res.residue_position;		
	residue_label_short = res.residue_label_short;	
	residue_label_long = res.residue_label_long;
	
	ix_residue = res.ix_residue;
	ix_chain = res.ix_chain;

	return *this;
}

PairResiude::~PairResiude()
{
	this;
	// no-op;
}

unsigned int cHairpin::GetStemLength()
{
	unsigned int cnt = 0;
	for (unsigned int i = 0; i < stem.pairs.size(); i++)
	{
		if (stem.pairs[i].r1.ix_residue >= 0 )
		{
			cnt++;
		}
		if (stem.pairs[i].r2.ix_residue >= 0 )
		{
			cnt++;
		}
	}

	return cnt;
}

cRNAStructure::cRNAStructure(string PdbFileName):cPDBStructure(PdbFileName)
{
	Init();
}

/*cRNAStructure::cRNAStructure(string PdbFileName, string chain):cPDBStructure(PdbFileName)
{
	Init();
}*/

cRNAStructure::cRNAStructure(string PdbFileName, string X3DNAPairingFileName, string chain, string exclusionList, string ranges):cPDBStructure(PdbFileName)
{	
	if (mModels.size() == 0) 
		//parsing of the PDB file did not result in a structure read
		return;
	if (chain != "")
	{
		//the chain can include multiple chains delimited by space, comma or semicolon

		vector<int> ixsChainToRemove;
		for (int i = 0; i < (int)mModels[0].size(); i++)
		{
			if (mModels[0][i].size() == 0 || chain.find(string(1,mModels[0][i][0].chain)) == string::npos)
			{
				ixsChainToRemove.push_back(i);
			}
		}

		if (mPrimarySeqChains.size() > 0)
		{
			for (int i = ixsChainToRemove.size()-1; i >= 0 ; i--)
			{
				for (int j = 0; j < (int)mModels.size(); j++)
				{
					mModels[j].erase(mModels[j].begin() + ixsChainToRemove[i]);
				}
				mPrimarySeqChains.erase(mPrimarySeqChains.begin() + ixsChainToRemove[i]);
			}
		}
	}

	//keep only residues which are specified in the ranges string
	if (ranges != "")
	{
		typedef vector<pair<int, int>> tRangeList;
		typedef boost::unordered_map<string, tRangeList> tRanges; //map of chains and associated vector of ranges
		tRanges mapRanges;
		//ranges are in the form C1-R1:C2-R2;C3-R3:C4-R4... R being the resiude sequence number and C being the chain id.
		vector< boost::iterator_range<string::iterator> > list;
		boost::split(list, ranges, boost::is_any_of(";"), boost::token_compress_on );
		for (int ixList = 0; ixList < list.size(); ixList++)
		{
			std::vector<std::string> range;
			boost::split(range, list[ixList], boost::is_any_of("-:"), boost::token_compress_on );
			if (range.size() != 4) throw "Wrong format ranges specification.";
			if (range[0] != range[2]) throw "Chains in a range must match.";
			if (mapRanges.count(range[0]) == 0) mapRanges.insert(make_pair(range[0], tRangeList()));
			mapRanges[range[0]].push_back(make_pair<int, int>(atoi(CleanseNumber(range[1]).c_str()),atoi(CleanseNumber(range[3]).c_str())));
		}
		for (int ixChain = 0; ixChain < (int)mModels[0].size(); ixChain++)
		{
			vector<int> ixsResToRemove;
			for (int ixRes = 0; ixRes < (int)mModels[0][ixChain].size(); ixRes++)
			{
				bool remove = true;

				string chain = string(1, mModels[0][ixChain][ixRes].chain);
				boost::to_upper(chain);
				int res_num = atoi(CleanseNumber(mModels[0][ixChain][ixRes].residue_num).c_str());

				for (tRangeList::iterator it = mapRanges[chain].begin(); it != mapRanges[chain].end(); ++it)
				{
					if (res_num >= it->first && res_num <= it->second)
					{
						remove = false;
						break;
					}
				}
				
				if (remove) ixsResToRemove.push_back(ixRes);
			}
			for (int ixRes = ixsResToRemove.size()-1; ixRes >= 0 ; ixRes--)
			{
				for (int ixModel = 0; ixModel < (int)mModels.size(); ixModel++)
				{
					mModels[ixModel][ixChain].erase(mModels[ixModel][ixChain].begin() + ixsResToRemove[ixRes]);
				}
			}
		}
	}

	//remove residues present in the exlusion list from the remaining chains 
	if (exclusionList != "")
	{
		exclusionList = ";" + exclusionList + ";";
		//exclusion list contains the list of residues to be excluded from the comparison in the form 
		//C1_R1;C2_R2;... R being the resiude id and C being the chain id		
		for (int ixChain = 0; ixChain < (int)mModels[0].size(); ixChain++)
		{
			vector<int> ixsResToRemove;
			for (int ixRes = 0; ixRes < (int)mModels[0][ixChain].size(); ixRes++)
			{
				string resId = boost::to_upper_copy(string(";") + mModels[0][ixChain][ixRes].chain + "-" + mModels[0][ixChain][ixRes].residue_num + ";");
				if (exclusionList.find(resId) != string::npos) ixsResToRemove.push_back(ixRes);
			}
			for (int ixRes = ixsResToRemove.size()-1; ixRes >= 0 ; ixRes--)
			{
				for (int ixModel = 0; ixModel < (int)mModels.size(); ixModel++)
				{
					mModels[ixModel][ixChain].erase(mModels[ixModel][ixChain].begin() + ixsResToRemove[ixRes]);
				}
			}
		}
	}	

	Init();
	if (mModelsNt.size() == 0) throw string("No nucleotide residues found in ") + PdbFileName + ". Check your specification.";

	//if (GetLengthNt() < 3) throw "Not enough nucleotide residues left for superposition.";

	ifstream ifs(X3DNAPairingFileName.c_str());
	if (ifs)
	{		
		vector<tPair> pairs;
		string line;
		while (!ifs.eof())
		{		
			//ofstream ofs("errAux");
			//ofs << PdbFileName << endl;
			getline(ifs, line);
			if (line.find("Strand I") != string::npos && line.find("Strand II") != string::npos && line.find("Helix") != string::npos)
			{
				getline(ifs, line);	
				//boost::trim(line);
				while (line != "" && line != "\r" && line != "\n" && line != "\r\n" && line[0] != '*')
				{
					//ofs << "a -" << line << "-";
					tPair	p;
					int		temp[10];
					string	str_temp;
					
					temp[0] = line.find_first_of(">");
					temp[1] = line.find_first_of(":");
					temp[2] = line.find_first_of(":", temp[1]+1);

					p.r1.chain = line.substr(temp[0]+1, temp[1]-temp[0]-1);
					//ofs << line << endl; ofs.flush();
					p.r1.residue_label_long = line.substr(29,3);
					//ofs << "xxx" << endl;ofs.flush();
					p.r1.residue_label_short = line.substr(33,1).c_str()[0];
					str_temp = line.substr(temp[1]+1, temp[2]-temp[1]-1);
					replace(str_temp.begin(), str_temp.end(), '.', ' ');
					replace(str_temp.begin(), str_temp.end(), '_', ' ');
					str_temp.resize(std::remove(str_temp.begin(), str_temp.end(), ' ') - str_temp.begin());
					p.r1.residue_position = str_temp;

					//identify position (ix_chain, is_residue) of the atom in the structure
					for (unsigned int i = 0; i < mModelsNt[0].size(); i++)
					{
						if (mModelsNt[0][i][0].chain == p.r1.chain[0])
						{
							p.r1.ix_chain = i;
							for (unsigned int j = 0; j < mModelsNt[0][i].size(); j++)
							{
								if (mModelsNt[0][i][j].residue_num == p.r1.residue_position)
								{
									p.r1.ix_residue = j;
									break;
								}
							}
						}
					}

					temp[0] = line.find_last_of(":");
					temp[1] = line.find_last_of("<");
					temp[2] = line.find_last_of(":", temp[0]-1);

					p.r2.chain = line.substr(temp[0]+1, temp[1]-temp[0]-1);
					p.r2.residue_label_long = line.substr(41,3);
					p.r2.residue_label_short = line.substr(39,1).c_str()[0];
					str_temp = line.substr(temp[2]+1, temp[0]-temp[2]-1);
					replace(str_temp.begin(), str_temp.end(), '.', ' ');
					replace(str_temp.begin(), str_temp.end(), '_', ' ');				
					str_temp.resize(std::remove(str_temp.begin(), str_temp.end(), ' ') - str_temp.begin());
					p.r2.residue_position = str_temp;

					//identify position (ix_chain, is_residue) of the atom in the structure
					for (unsigned int i = 0; i < mModelsNt[0].size(); i++)
					{
						if (mModelsNt[0][i][0].chain == p.r2.chain[0])
						{
							p.r2.ix_chain = i;
							for (unsigned int j = 0; j < mModelsNt[0][i].size(); j++)
							{
								if (mModelsNt[0][i][j].residue_num == p.r2.residue_position)
								{
									p.r2.ix_residue = j;
									break;
								}
							}						
						}
					}

					temp[0] = line.find_first_of("]");
					p.pair_type = line.substr(temp[0]+2, 5);
					p.spatial_information = line[71];

					
					pairs.push_back(p);
					getline(ifs, line);
				}
			}
		}

		for (unsigned int i = 0; i < pairs.size(); i++)
		{
			for (int k = 0; k < 2; k++)
			{
				tPair p;
				if (k ==0)
				{
					p = pairs[i];
				}
				else
				{
					p = SwitchResiduesInPair(pairs[i]);
				}
				
				string chain = p.r1.chain;
				if (mPairsChains[chain].empty())
				{
					mPairsChains[chain].push_back(p);					
				}
				else
				{
					unsigned int j;
					for (j = 0; j < mPairsChains[chain].size(); j++)				
					{
						if (p.r1.ix_residue <= mPairsChains[chain][j].r1.ix_residue)
						{
							mPairsChains[chain].insert(mPairsChains[chain].begin()+j, p);
							break;
						}
					}
					if (j >= mPairsChains[chain].size())
					{
						mPairsChains[chain].push_back(p);
					}
				}
				if (mPairs.empty())
				{
					mPairs.push_back(p);
				}
				else
				{	
					unsigned int j;
					for (j = 0; j < mPairs.size(); j++)				
					{
						if (p.r1.ix_residue <= mPairs[j].r1.ix_residue && p.r1.ix_chain == mPairs[j].r1.ix_chain)
						{
							mPairs.insert(mPairs.begin()+j, p);
							break;
						}
					}
					if (j >= mPairs.size())
					{
						mPairs.push_back(p);
					}
				}
			}
		}
	}
	else
	{
		throw (X3DNAPairingFileName + " could not be opened.").c_str();
	}
	
	ExtractHairpins();
}

tPair cRNAStructure::SwitchResiduesInPair(tPair p)
{
	tPair ps;

	ps = p;
	ps.r1 = p.r2;
	ps.r2 = p.r1;

	return ps;
}

void cRNAStructure::Init()
{
	ExtractSimplified();
}

void cRNAStructure::ExtractSimplified()
{
	for (unsigned int i = 0; i < mModels.size(); i++)
	{
		tPdbModel mNew;
		tPdbModel m = mModels[i];

		for (unsigned int j = 0; j < m.size(); j++)
		{
			tPdbChain chNew;
			tPdbChain ch = m[j];

			for (unsigned int k = 0; k < ch.size(); k++)
			{
				tPdbAtom atom = ch[k];
				//if (atom.name == " C4'")
				if (atom.name == " P  " || atom.name == " C4'")
				{ //we prefer P atoms but if the residue misses a P atom, we take C4
					if (chNew.size() > 0 && boost::trim_copy(chNew[chNew.size() - 1].residue_num) == boost::trim_copy(atom.residue_num) && chNew[chNew.size() - 1].ins_code == atom.ins_code)
					{
						if (atom.name == " P  ")
						{
							chNew.erase(chNew.end() - 1);
							chNew.push_back(atom);
						}

					}
					else
					{
						chNew.push_back(atom);						
					}					
				}
			}

			if (!chNew.empty())
			{
				int l = 0;
				for (; l < (int)mNew.size(); l++)
				{
					if (mNew[l][0].chain == chNew[0].chain)
					{
						break;
					}
				}
				if (l == (int)mNew.size())
				{
					mNew.push_back(chNew);
				}
				else
				{
					for (int m = 0; m < (int)chNew.size(); m++)
					{
						mNew[l].push_back(chNew[m]);
					}
				}
			}
		}

		if (!mNew.empty())
		{
			mModelsNt.push_back(mNew);
		}
	}
}

// if atom is in processed at position r1 or r2, its position in the vector is returned ... else -1 is returned
int AtomInPair(tPdbAtom atom, vector<tPair> pairs, int ixResidue = 1)
{
	for (unsigned int i = 0; i < pairs.size(); i++)
	{
		if (
			(ixResidue == 1 && pairs[i].r1.ix_residue >= 0 && atom.chain == pairs[i].r1.chain[0] && atom.residue_num == pairs[i].r1.residue_position)
			||
			(ixResidue == 2 && pairs[i].r2.ix_residue >= 0 && atom.chain == pairs[i].r2.chain[0] && atom.residue_num == pairs[i].r2.residue_position)
			)
		{
			if (pairs[i].spatial_information == '|' /*&& (pairs[i].pair_type == "-----" || pairs[i].pair_type == "x----" || pairs[i].pair_type == "----x") */)
			{
				return i;
			}
		}
	}
	return -1;	
}

tPairResiude GetPairResidueFromAtom(tPdbAtom atom, int ixChain, int ixResidue)
{
	tPairResiude r;
	r.chain = atom.chain;
	r.ix_chain = ixChain;
	r.ix_residue = ixResidue;
	r.residue_label_short = atom.res_name;
	r.residue_position = atom.residue_num;

	return r;
}

void cRNAStructure::ExtractHairpins(bool GSSUsInSeqOrder)
{
	if (mModelsNt.size() == 0) return;

	tPdbModel m = mModelsNt[0];
	vector<tPair> processed;
	int auxGSSUSeqIx = 0;
	bool in_stem = false;
	cHairpin hp;
	vector<int> postOrderLabeling;
	
	vector<unsigned int> prcsdCounts; // when a stem is finished, number of elements in processed vector is logged 
									//thus we can recognize stems of e.g. cloverleaf
									//morevoer it can help in dealing with some structures such as 2TPK, withou that
									// it would put residues 1, 2, 28-36 into one GSSU
	for (unsigned int i = 0; i < m.size(); i++)
	{
		tPdbChain ch = m[i];
		vector<unsigned int> prcsdCounts; // when a stem is finished, number of elements in processed vector is loged
												//thus we can recognize stems of e.g. cloverleaf
		for (unsigned int j = 0; j < ch.size();  j++)
		{
			tPdbAtom atom = ch[j];
			int ixProcessed = -1;
			//int ixPairs = AtomInPair(atom, mPairs[string(1, atom.chain)]); //if hairpins should contain residues coming from one chain only
			int ixPairs = AtomInPair(atom, mPairs);
			if (ixPairs >= 0)
			{
				ixProcessed = AtomInPair(atom, processed, 2);
			}

			if (in_stem)
			{				
				if (ixPairs >= 0)
				{
					if (ixProcessed >= 0)
					{
						if (prcsdCounts.size() > 0 && processed.size() == prcsdCounts[prcsdCounts.size()-1])
						{
							if (mHairpins.empty() || !GSSUsInSeqOrder || hp.stem.pairs.empty())
							{
								mHairpins.push_back(hp);
								postOrderLabeling.push_back(auxGSSUSeqIx);
							}
							else 
							{
								mHairpins.insert(mHairpins.begin() + hp.stem.pairs[0].aux, hp);
								postOrderLabeling.insert(postOrderLabeling.begin() + hp.stem.pairs[0].aux, auxGSSUSeqIx);
							}
							hp.head.residues.clear();
							hp.stem.pairs.clear();
							prcsdCounts.pop_back();
							auxGSSUSeqIx++;
						}

						for (int k = processed.size()-1; k >= ixProcessed; k--)
						{
							hp.stem.pairs.push_back(processed[k]);
						}						
						processed.erase(processed.begin() + ixProcessed, processed.end());
					}
					else
					{
						if (mHairpins.empty() || !GSSUsInSeqOrder || hp.stem.pairs.empty())
						{
							mHairpins.push_back(hp);
							postOrderLabeling.push_back(auxGSSUSeqIx);
						}
						else
						{
							mHairpins.insert(mHairpins.begin() + hp.stem.pairs[0].aux, hp);
							postOrderLabeling.insert(postOrderLabeling.begin() + hp.stem.pairs[0].aux, auxGSSUSeqIx);
						}
						hp.head.residues.clear();
						hp.stem.pairs.clear();
						in_stem = false;
						auxGSSUSeqIx++;
						if (prcsdCounts.empty() || prcsdCounts[prcsdCounts.size()-1]<processed.size())
						{
							prcsdCounts.push_back(processed.size());
						}
						//processed.push_back(mPairs[string(1,atom.chain)][ixPairs]); //if hairpins should contain residues coming from one chain only
						processed.push_back(mPairs[ixPairs]);
						
					}	
				}
				else
				{
					tPair p;
					p.r2 = GetPairResidueFromAtom(atom, i, j);
					p.r1.ix_residue = -1;

					hp.stem.pairs.push_back(p);
				}
			}
			else
			{
				if (ixPairs >= 0)
				{
					if (ixProcessed >= 0)
					{ //not in stem and paired element discovered whose second element has already been touched => neck (end/start of stem)
						hp.head.residues.insert(hp.head.residues.begin(), processed.begin() + ixProcessed + 1, processed.end());
						hp.stem.pairs.push_back(processed[ixProcessed]);
						processed[processed.size() - 1].aux = auxGSSUSeqIx;
						processed.erase(processed.begin() + ixProcessed , processed.end());
						in_stem = true;
					}
					else
					{
						//processed.push_back(mPairs[string(1,atom.chain)][ixPairs]); //if hairpins should contain residues coming from one chain only
						processed.push_back(mPairs[ixPairs]);
						processed[processed.size() - 1].aux = auxGSSUSeqIx;
					}
				}
				else
				{
					tPair p;
					p.r1 = GetPairResidueFromAtom(atom, i, j);	
					p.r2.ix_residue = -1;
					p.aux = auxGSSUSeqIx;

					processed.push_back(p);
				}
			}
		}

		//this piece of code should be uncommented if hairpins contain reisudes from one chain only
		/*if (hp.head.residues.size() > 0 || hp.stem.pairs.size() > 0)
		{
			mHairpins.push_back(hp);
			hp.head.residues.clear();
			hp.stem.pairs.clear();
		}
		
		if (processed.size() > 0)
		{
			hp.head.residues.clear();
			hp.stem.pairs.clear();
			hp.stem.pairs.insert(hp.stem.pairs.begin(), processed.begin(), processed.end());
			mHairpins.push_back(hp);
			hp.head.residues.clear();
			hp.stem.pairs.clear();

			processed.clear();
		}*/

		
	}
	if (processed.size() > 0)
	{
		/*
		//this works well for pseudoknots but not for the "normal" RNA structure
		//the presumption for the normal RNA structure is that the whole structure forms a sort of hairpin, e.g.,
		//the 5' and 3' ends are approximately at the same position in 2D visualziation but this has not be necesarily
		//true for pseudoknots (2TPK). In such case the beginning and end form one GSSU which looks strange (at least
		//from the view of the 2D image)
		while (!prcsdCounts.empty())
		{
		int ixAux = prcsdCounts[prcsdCounts.size()-1];
		hp.stem.pairs.insert(hp.stem.pairs.begin(), processed.begin() + ixAux, processed.end());
		mHairpins.push_back(hp);
		hp.head.residues.clear();
		hp.stem.pairs.clear();

		processed.erase(processed.begin() + ixAux, processed.end());
		prcsdCounts.pop_back();
		}
		*/
		hp.stem.pairs.insert(hp.stem.pairs.begin(), processed.begin(), processed.end());

		processed.clear();
	}

	if (hp.head.residues.size() > 0 || hp.stem.pairs.size() > 0)
	{
		if (mHairpins.empty() || !GSSUsInSeqOrder || hp.stem.pairs.empty())
		{
			mHairpins.push_back(hp);
			postOrderLabeling.push_back(auxGSSUSeqIx);
		}
		else
		{
			mHairpins.insert(mHairpins.begin() + hp.stem.pairs[0].aux, hp);
			postOrderLabeling.insert(postOrderLabeling.begin() + hp.stem.pairs[0].aux, auxGSSUSeqIx);
		}
		hp.head.residues.clear();
		hp.stem.pairs.clear();
	}

	//mGSSUPostOrderLabeling.insert(pair<std::string, vector<int>>(string(1, ch[0].chain), postOrderLabeling));
	
	

	//extract hairpins with respect to chain
	for (int i = 0; i < (int)mHairpins.size(); i++)
	{
		cHairpin hp = mHairpins[i];
		map<string, int> presentChains;
		for (int j = 0; j < (int)hp.head.residues.size(); j++)
		{
			presentChains[hp.head.residues[j].r1.chain] = 1;
			presentChains[hp.head.residues[j].r2.chain] = 1;
		}
		for (int j = 0; j < (int)hp.stem.pairs.size(); j++)
		{
			presentChains[hp.stem.pairs[j].r1.chain] = 1;
			presentChains[hp.stem.pairs[j].r2.chain] = 1;
		}
		map<string, int>::iterator itPC;
		for (itPC = presentChains.begin(); itPC != presentChains.end(); itPC++)
		{
			if (itPC->first != "")
			{
				mHairpinsChains[itPC->first].push_back(hp);
			}			
		}
	}
}

string cRNAStructure::GetPrimarySeqFromModel(string chain)
{
	string seq;
	
	int ixChain = 0;
	if (chain != "")
	{
		for (int i = 0; i < (int)mModelsNt[0].size(); i++)
		{
			if (string(1,mModelsNt[0][i][0].chain) == chain)
			{
				ixChain = i;
				break;
			}
		}
	}

	for (int i = 0; i < (int)mModelsNt[0][ixChain].size(); i++)
	{
		string auxStr = mModelsNt[0][ixChain][i].res_name;
		std::stringstream trimmer;
		trimmer << auxStr;
		trimmer >> auxStr;
		seq += auxStr;
	}

	return seq;
}

vector<tPdbAtom> cRNAStructure::GetAllNtResidues()
{
	vector<tPdbAtom> residues;

	for (int ixChain = 0; ixChain < mModelsNt[0].size(); ixChain++)
	{ //for each chain
		for (int ixRes = 0; ixRes < mModelsNt[0][ixChain].size();ixRes++) residues.push_back(mModelsNt[0][ixChain][ixRes]);
	}
	return residues;
}

unsigned int cRNAStructure::GetLengthNt()
	{
		unsigned int size = 0;
		if (mModelsNt.size() > 0)
		{
			for (int i = 0; i < (int)mModelsNt[0].size(); i++)
			{
				size += mModelsNt[0][i].size();
			}
		}

		return size;
	}

void cRNAStructure::RenameWithSeqOrder()
{
	/*
	int id = 0;
	for (int ixChain = 0; ixChain <mModelsNt[0].size(); ixChain++)
	{ //for each chain
		string chainName = string(1,mModelsNt[0][ixChain][0].chain);
		for (int ixHp = 0; ixHp < mHairpinsChains[chainName].size(); ixHp++)
		{
			//first process the left side of the stem
			tHPStem stem = mHairpinsChains[chainName][ixHp].stem;
			for (int ixStem = 0; ixStem < stem.pairs.size(); ixStem++)
			{
				tPairResiude pr = stem.pairs[ixStem].r1;
				if (pr.ix_residue >= 0) mModels[0][pr.ix_chain][pr.ix_residue].residue_num = id++;
			}
			//second extract the loop, if present
			tHPHead head = mHairpinsChains[chainName][ixHp].head;
			if (head.residues.size() > 0)
			{
				for (int ixHead = 0; ixHead < head.residues.size(); ixHead++)
				{
					tPairResiude pr = head.residues[ixHead].r1;
					if (pr.ix_residue >= 0) mModels[0][pr.ix_chain][pr.ix_residue].residue_num = id++;
				}

			}
			//if the loop was present, right side of the stem can be processed, 
			//otherwise it needs to be processed only after the GSSUs between left
			//and right side of the stem will be processed
				
			
		}
	}

	//Extract simplified?
	return residues;
	*/

}