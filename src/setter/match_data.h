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

#include <vector>
#include <map>

#include "cRNAStructure.h"
#include "common.h"

//! The class is used to store matching informations.
/*! 
	The class stores the matching results between two RNA and both between GSSUs 
	of two matched RNA structures cretaed by original SETTER algorithm. <br>
	The class supports lazy initialization
*/
class MatchData
{
	//friend class Input;
	friend class FormAvg;
	friend class ParallelMakeBackup;
	friend class ParallelFullQVect;
	friend class ParallelPreprocessRNA;
	friend class ParallelMerge;
	friend class NJInput;
	//friend class ParallelAddRNA;
public:
	//! Default constructor.
	MatchData();

	//! Constructor
	/*!
		\param rnaQ
		\param rnaDB
		\param init	decides if the constructor has to initializes the data fields
	*/
	MatchData(
		cRNAStructure *rnaQ, 
		cRNAStructure *rnaDB,
		bool init = true);

	//! Retrieves the data from a vecor of RNA structures and used positions in the vector.
	/*!
		\param rnaVect	A vector containing RNA structures
		\param rnaQIdx	the id of rnaQ in the vector
		\param rnaDBIdx	the id of rnaDB in the vector
		\param init	decides if the constructor has to initializes the data fields
	*/
	MatchData(
		vector<cRNAStructure *> &rnaVect, 
		int rnaQIdx, 
		int rnaDBIdx,
		bool init = true);
	
	//! Stores the string representation of RNA indexes.
	/*!
		\param rnaQ		the rnaQ
		\param rnaDB	the rnaDB
		\param rnaQIdx	the string representation of rnaQ idenfiier
		\param rnaDBIdx	the string representation of rnaDB idenfiier
		\param init	decides if the constructor has to initializes the data fields
	*/
	MatchData(
		cRNAStructure *rnaQ, 
		cRNAStructure *rnaDB, 
		string rnaQIdx, 
		string rnaDBIdx,
		bool init = true);
	
	//! Copy constructor
	MatchData(MatchData &data);

	//! Assignment operator.
	void operator=(MatchData &data);
	
	//! returns the matching results between GSSUs
	/*!
		\param i	the order number of GSSU in rnaQ
		\param j	the order number of GSSU in rnaDB
		\result		the matching results of GSSUs
	*/
	sResMatch operator()(int i, int j);
	~MatchData();

	//! Returns the matching information between the two RNA input
	sResMatch &GetResult();

	//! the name of the Q RNA
	string sRnaQIdx;
	//! the name of the DB RNA
	string sRnaDBIdx;

private:
	//! calls match for RNAs, sorts the optPairs (from lowes q)
	void Initalize();
	
	//! pointer to query RNA
	cRNAStructure *mRnaQ;
	//! pointer to db RNA
	cRNAStructure *mRnaDB;
	//! number of hairpins in rnaQ
	int mCntHpQ;
	//! number of hairpins in rnaDB
	int mCntHpDB;

	//! default parameters
	static sParams mDefParams;

	//! the result of the RNA match
	sResMatch mResult;
	
	//! the pairing of GSSU pairs (contains positions)
	vector<GSSUPair> mOptPairs;

	//! the result of match each GSSUs
	map<pair<int, int>, sResMatch> mRessMap;

	//! true, if the data is initialized, used at layu initializeation.
	bool mInit;
};

