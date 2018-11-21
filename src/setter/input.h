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
#include <string>
#include <vector>
#include "cRNAStructure.h"
#include "common.h"
#include "match_data.h"
#include "parallel.h"

using namespace std;

class MatchData;

//! Structure containing information from input file.
struct InputRNA
{
	//! Id of the structure.
	string id;
	//! Path to pdb file.
	string pdb;
	//! Path to 3dna file.
	string x3dna;
	//! Chain ID
	string chain;
	//! Comments
	string comment;
	//! Ranges of residues to be considered.
	string ranges;
	//! Residues not to be considered.
	string exclusions;
};

//! The calss is used to preprocess and store the input RNA structures
class Input
{
	friend class boost::serialization::access;
	friend class ParallelResult;
public:
	template<class Archive>
	//! Used to BOOST serialization.
    void serialize(Archive & ar, const unsigned int version);
	//! Default constructor
	Input();
	//! Constructs the from a path to input file
	/*!
		\param file	The path to the input file.
	*/
	Input(string file);
	//! Copy constructor.
	Input(Input &in);
	//! Assignment operator.
	Input& operator = (const Input &in);
	//! Constructs from a vector of RNA structures, doesn't have all the informations.
	/*!
		\param inputRnaVect	A vector of RNA structures.
	*/
	Input(vector<cRNAStructure *> &inputRnaVect);
	//! Constructs from a vector of RNA definitions.
	/*!
		\param inputRnaDefVect	A vector of RNA definitions.
	*/
	Input(vector<InputRNA> &inputRnaDefVect);
	~Input();
	//! Removes an RNA from given position in the mInputRnaVect
	/*!
		\param pos The position of RNA in mInputRnaVect, which will be removed.
	*/
	void RemoveRNA(int pos); 

	//! The size of the input (number of RNAs).
	int mInputSize;
	//! Contains RNAs read from input file
	vector<cRNAStructure *> mInputRnaVect;
	//! Contains the names of RNAs read from the file in the same order as in mInputRnaVect.
	vector<string> mOriginalMap;
};

//! Class inherited from Input, used as input for Neighbour-Joinning
class NJInput: public Input
{
	friend class MultiAlign;
	friend class MultiAlignImpl0;
	friend class MultiAlignImpl1;
	friend class MultiAlignImpl2;
	friend class MultiAlignImpl3;
	friend class ParallelResult;

public:
	//! Default constructor.
	NJInput();
	//! Constructor from base class.
	NJInput(Input &in);
	//! Copy constructor.
	NJInput(NJInput &input);
	~NJInput();

	//! Returns the number of RNA structures stored in this data structure.
	int GetSize();
	//! Returns the matchin data created by SETTER.
	/*!
		The returned RNA will be one of the RNAs with indexes from input parameters, 
		which is closer to other RNAs.
		\param i	The index of RNA in mRnaVect which could be the result of method.
		\param j	The index of RNA in mRnaVect which could be the result of method.
		\result		Matching results of RNA, which is closer to other structures in Input
	*/
	MatchData operator()(int i, int j);
	//! Returns a string identificator of rnaQ; it is the one from ith and jth, that is closer to other structures.
	string GetRnaQIdx(int i, int j);
	//! Returns a string identificator of rnaQ; it is the one from ith and jth, that is more far to other structures.
	string GetRnaDBIdx(int i, int j);
	//! Returns the score of matching between RNAs with indexes from  the input parameters.
	/*!
		The result is calculated according the Parameters of algorithm. (min, arithmetical mean, geometrical mean).
	*/
	double CalcScore(int i, int j);
	//! Allmost the same as operator().
	void GetData(int i, int j, MatchData &data);

private:
	//! Calculates an order of RNAs according the distance from other RNAs in the input set.
	void SortData();
	//! Calculates the mean distance between GSSUs
	void CalcMean();
	//! Prints the all-to-all distances
	void Print();
	//! Removes RNA from with given index.
	void RemoveRNA(int idx);
	//! Adds RNA
	/*!
		\param rna		The RNA structure whichwill be added to the input set.
		\param rnaIdx	The name of the RNA.
		\param init		If false, lazy initialization will be used.
	*/
	void AddRNA(
		cRNAStructure *rna, 
		string rnaIdx = "newRNA",
		bool init = true);
	//! Delets data fields.
	void DeleteData();

	//! Number of RNA structures in mRnaVect and mMap.
	int mSize;
	//! At first it's copy of mInputRnaVect, then RNAs will be deleted, added.
	vector<cRNAStructure *> mRnaVect;
	//! Names of RNAs in mRnaVect.
	vector<string> mStringMap;
	//! Map between actual indexes of RNA structures and original indexes.
	vector<int> mIntMap;
	//! Matchin data which will be filled by SETTER.
	MatchData ***mData;

	//! Sorted list of RNAs by distance from each other.
	vector<int> mSortedList;
	//! Mean distance between GSSUs, calculated by CalcMean methode.
	double mMean;
};

//! Static parses input file.
class Parser
{
private:
	//! Private constructor.
	Parser();
public:
	//!Parses inpu file.
	/*!
		\param file		The input file.
		\param vect		The vector of RNAs, this vector will be filled by this methode.
		\param map		The vector of names of RNAs, this vector will be filled by this methode.
	*/
	static void ParseInput(
		string file, 
		vector<cRNAStructure *> &vect,
		vector<string> &map);
	static void ParseInput(string file, cRNAStructure &rna);
};

//! Intel TBB class which reads RNA structures.
class ParallelRead
{
public:
	ParallelRead(
		vector<InputRNA> &inputRNAs,
		vector<cRNAStructure *> &vect,
		vector<string> &map);

	void operator()(const tbb::blocked_range<size_t> &r) const;
	void serial(int size) const;
	
private:
	vector<InputRNA> &mInputRNAs;
	vector<cRNAStructure *> &mVect;
	vector<string> &mMap;
};