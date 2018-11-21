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
#ifndef CRNASTRUCTURE_H_
#define CRNASTRUCTURE_H_

#include "pdb-parser/pdb_parser.h"
#include <map>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/base_object.hpp>

typedef struct PairResiude
{
	template<class Archive>
    void serialize(Archive & ar, const unsigned int version);

	PairResiude();// {ix_residue = -1;}
	PairResiude(const PairResiude &res);
	PairResiude &operator=(const PairResiude &res);

	~PairResiude();

	std::string	chain;
	std::string	residue_position;		//std::string, since eg. position 25A is possible!!
	std::string	residue_label_short;	//should by A, C, G, U but could be e.g. DA, DC, DG, DU
	std::string residue_label_long;
	
	int		ix_residue;	//index in the simplified representation (mModelsNt)
	int		ix_chain;
} tPairResiude;

typedef struct 
{
	template<class Archive>
    void serialize(Archive & ar, const unsigned int version);

	tPairResiude	r1, r2;
	std::string		pair_type;
	char			spatial_information; // | ... normal, + ... isolated BP, x ... helix change
	int				aux; //to be used when building GSSUs, if they should be sorted in sequence order
} tPair;

typedef struct 
{
	template<class Archive>
    void serialize(Archive & ar, const unsigned int version);

	vector<tPair>	pairs;
} tHPStem;

typedef struct
{
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version);

	vector<tPair> residues; //second residue will remain empty

} tHPHead;

class cHairpin{
private:
	friend class boost::serialization::access;

	template<class Archive>
    void serialize(Archive & ar, const unsigned int version);

public:
	tHPStem		stem;
	tHPHead		head;	

	unsigned int GetStemLength();

	inline unsigned int GetHeadLength()
	{
		return head.residues.size();
	}

	//number of residues (does it is different number than |head| + |stem| since stem contains pairs)
	inline unsigned int GetLength()
	{
		return GetHeadLength() + GetStemLength();
	}

	unsigned int GetCntBonds()
	{
		int cnt = 0;
		for (unsigned int i = 0; i < stem.pairs.size(); i++) if (stem.pairs[i].r1.ix_residue >= 0 || stem.pairs[i].r2.ix_residue >= 0) cnt ++;
		return cnt;

	}
};


class cRNAStructure:public cPDBStructure
{
private:
	friend class FormAvg;
	friend class ParallelMakeBackup;
	friend class ParallelFullQVect;
	friend class ParallelPreprocessRNA;
	friend class ParallelMerge;
	friend class ParallelMerge;
	friend class Printer;
	friend class boost::serialization::access;

	template<class Archive>
    void serialize(Archive & ar, const unsigned int version);

	void	Init();
	void	ExtractSimplified();
	tPair	SwitchResiduesInPair(tPair);

protected:

	std::string								mPDBFileName;

	std::string								X3DNAPairingFileName;
	map< std::string, vector<tPair> >		mPairsChains; //one vector for each chain
	vector<tPair>							mPairs; //one vector for each chain

	vector<cHairpin>						mHairpins;
	map< std::string, vector<cHairpin> >	mHairpinsChains; //one vector for each chain
	map< std::string, vector<int> >			mGSSUPostOrderLabeling; //one vector for each chain

	vector<tPdbModel>						mModelsNt;	//simplified representation containing only one atom (C4) per residue
	

	void ExtractHairpins(bool GSSUsInSeqOrder = false);
	void RenameWithSeqOrder();

public:
	
	cRNAStructure(){};
	cRNAStructure(std::string PdbFileName);
	cRNAStructure(std::string PdbFileName, std::string X3DNAPairingFileName, std::string chain = "", std::string exclusionList = "", std::string ranges = "");
	~cRNAStructure(){};

	inline tPdbAtom		GetAtom(tPairResiude r){ 
		if (r.ix_chain > mModelsNt[0].size() - 1) {
			throw "Chain index out of bounds.";
		}
		else if (r.ix_residue > mModelsNt[0][r.ix_chain].size() - 1) {
			cout << "Residue index out of bounds." << endl;
		}
		return mModelsNt[0][r.ix_chain][r.ix_residue]; 
	};
	inline unsigned int GetHairpinsCount(){ return mHairpins.size(); };	
	inline unsigned int GetHairpinsCount(std::string chain){ return mHairpinsChains[chain].size(); };
	inline cHairpin		GetHairpin(int i){ return mHairpins[i]; };
	inline cHairpin		GetHairpin(std::string chain, int i){ return mHairpinsChains[chain][i]; };
	std::string			GetPrimarySeqFromModel(std::string chain = "");

	std::vector<tPdbAtom>	GetAllNtResidues();
	void Relabel(); //used in MultiSETTER when a new average structure arises simply by putting together merged GSSUs
	
	unsigned int GetLengthNt();

	//cRNAStructure *Segregate(char chain, std::string start, std::string end);
};

#endif /*CRNASTRUCTURE_H_*/
