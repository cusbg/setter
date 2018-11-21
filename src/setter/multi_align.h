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
#include "input.h"

using namespace std;
//! The result structure of Multipe Alignment.
struct MAResult {

	template<class Archive>
    //! Serialization
	void serialize(Archive & ar, const unsigned int version);

	//! the average RNA created from a set of RNA structures.
	cRNAStructure rna;
	//! Contains the pdb identificators and chains of input RNA structures..
	vector<string> map;
	//! Distances of input RNAs from the average RNA.
	vector<double> dist;
	//! Contains Matching Results with originals
	vector<ResMatch> results;
	
	//! Mean distance of input RNAs from the average RNA.
	double mean;
	//! Deviation of distances between input RNAs and the average RNA.
	double deviation;

	//! Deafult constructor.
	MAResult ();
	//! Constructor from building blocks.
	MAResult (
		cRNAStructure &rna, 
		vector<string> &map, 
		vector<double> &dist,
		vector<ResMatch> &results);
};

//! Base class for Multiple Alignment algorithms.
/*!
	Implemeted Multiple Alignment algorithms differs in the implementation of 
	Neighbour Joinning algorithm used to calculate the guide tree.
	<br>
	The implementation uses strategy pattern and factory methode.
 */
class MultiAlign
{
public:
	//! virutal funcon which has to implement the inherited classes.
	virtual MAResult GetResult() = 0;
	//! Returns the score of the RNA merging. NOT USED
	double GetScore();
	//! Factory method returning the concrete implementation of the MultipleAlignment.
	/*!	\param id the id of the algorithm <Parameter.Algo>
		\param the input data
		\result Concrete implementation of the Multiple Alignment.
	*/
	static MultiAlign *Factory(int id, Input &in);

protected:
	//! Constructor which stores reference for input data.
	MultiAlign(Input &in);
	//! Reference for input data.
	Input &mIn;

	int cnt;
	double score;
};

//! Multiple Alignment class using standard neighbor joinning method.
class MultiAlignImpl0 : public MultiAlign
{
public:
	MultiAlignImpl0(Input &in);
	MAResult GetResult();
};

//! Multiple Alignment class using modified neighbor joinning method.
/*!
	All the merged nodes are added to the set of input structures and the distance matrix is recalculated.
*/
class MultiAlignImpl1 : public MultiAlign
{
public:
	MultiAlignImpl1(Input &in);
	MAResult GetResult();
};

//! Multiple Alignment class using modified neighbor joinning method, NOT USED
class MultiAlignImpl2 : public MultiAlign
{
public:
	MultiAlignImpl2(Input &in);
	MAResult GetResult();
};

//! Class for distance matrix calculation
class DistanceMatrix
{
	friend class QMatrix;
public:
	//! Default constructor
	DistanceMatrix();
	//! Constructs distance matrix from the set input structures
	DistanceMatrix(NJInput &input);
	//! Copy constructor
	DistanceMatrix(DistanceMatrix &ref);
	~DistanceMatrix();

	//! Returns the distance between from the given cell.
	double operator()(int i, int j);

	//! Calculates the philogenetic distance between two node
	void CalcDist(
		pair<int, int> min, 
		int &cnt, 
		double &score);

	//! Removes a node from distance matrix.
	void RemoveRNA(int idx);
	
	//! Adds new node to the distance matrix.
	/*! Respectively pairs the nodes i, j and adds a new node to the matrix.
	*/
	void AddRna(
		DistanceMatrix &copyDM, 
		pair<int, int> min);

private:
	int mSize;
	NJInput mInput;
	double **matrix;
};

//! class for Q matrix for reducing the effects of long edges
class QMatrix
{
public:
	QMatrix();
	QMatrix(DistanceMatrix &distMatrix);
	~QMatrix();

	pair<int, int> GetMinPair();

	double **matrix;

private:
	void Calculate();
	double Calculate(int i, int j);

	int mSize;
	DistanceMatrix mDistMatrix;
	pair<int, int> mMinPair;
};

//! Intel TBB class to compute the distances between the parent RNAs and the newly merged RNA
class ParallelResult
{
public:
	ParallelResult(
		NJInput &input,
		cRNAStructure &newRna,
		vector<ResMatch> &resultVector,
		vector<double> &dist);

	void operator()(const tbb::blocked_range<size_t> &r) const;
	void serial(int size) const;

private:
	NJInput &mInput;
	cRNAStructure &mNewRna;
	vector<ResMatch> &mResultVector;
	vector<double> &mDist;
};