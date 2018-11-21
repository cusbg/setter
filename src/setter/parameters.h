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
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <string>

#ifndef TBB
#define TBB 1
#endif

#if TBB

//! Paralelly calculates the distance matrix.
#ifndef PARALLEL_INPUT
#define PARALLEL_INPUT 1
#endif

//! Parallel SETTER
#ifndef PARALLEL
#define PARALLEL 1
#endif

#endif //TBB

typedef boost::archive::text_iarchive i_archive;
typedef boost::archive::text_oarchive o_archive;

using namespace std;

//! Singleton class containing the algorithms parameters.
class Parameters 
{
public:
	//! The Neighbor Joining algorithm implementations
	/*!
		Three Neighbour Joining algorithms are implemented:
	*/
	enum Algo {
		STANDARD,	/*!< The standard Neighbour Joining algorithm. */
		RECALC,		/*!< Allways merges the best candidate and recalculates the guide tree*/
		ALL 		/*!< Merges as many pairs, as posible, then puts back average RNAs to the
							input set and recalculates the guide tree -> probaby biologically incorrect.*/
	};

	//! Determines how to calculate the final score from two scores.
	/*!
		The SETTER algorithm is not symetric (match(A, B) != match(B, A)).<br>
		This enum defines how to sole this problem.<br>
		There are three solutions implemeted:
	*/
	enum Mean {
		MIN,			/*!< Selects the minimum from two input */
		ARITHM,			/*!< Arithmetic mean */
		GEOM			/*!< Geometric mean	*/
	};

	enum Merge {
		ALTERNATING,
		PATTERN
	};

private:
	//! Pattern in [parameters file] to parse the number of threads.
	static const string THREADS;
	//! Pattern in [parameters file] to parse the whether half rotation of both RNAs or full rotation of single RNA is used..
	static const string HALF_ROTATION;
	//! Pattern in [parameters file] to parse the neighbour joining algorithm type.
	static const string ALGO;

	static const string MERGE;
	//! Pattern in [parameters file] to parse the mean type.
	static const string MEAN;
	//! Pattern in [parameters file] to parse the parameter of allowed divergence from mean distance between GSSUs.
	static const string DIVER_MEAN;
	//! Pattern in [parameters file] to parse the parameter of allowed divergence from the parent GSSUs.
    static const string DIVER_PARENT;
	//! Pattern in [parameters file] to parse the parameter of JAVA highliting.
	static const string HIGHLIGHT;
	//! pdb printer enabled/disabled
	static const string SCRIPT;

	//! The only instance of the singleton.
	static Parameters mInstance;
	//! Determines whether the singleton is initialized.
	static bool mInitalized;

	//! The number of working threads.
	int mThreads;
	//! Whether use half rotation or full rotation.
	bool   mHalfRotation;
	//! The type of neighbour joining algorithm
	Algo   mAlgo;
	//! The type of mean.
	Mean   mMean;
	//! The type of the merging algorithm
	Merge mMerge;
	//! Multiplicator for allowed divergence from mean distance  between GSSUs.
	double mDiverMean;
	//! Multiplicator for allowed divergence from the parent GSSUs.
	double mDiverParent;
	//! The highlither string to Java.
	string mHighlight;
	//! PDB print enabled/disabled
	bool mScriptOutput;

	Parameters();
public:
	//! Initializes the algorithms parameters from input file.
	/*!
		\param file The file containing the parameters of the algorithm.
	*/
	static void Init(string file);

	//! Returns the instance of the singleton.
	static Parameters &GetInstance();

	//! Returns the number of working threads.
	int GetThreads();
	//! Returns whether use half rotation on both RNA or full rotation on single RNA
	bool GetHalfRotation();
	//! Returns the type of the Neighbour Joining algorithm.
	Algo GetAlgorithm();
	//! Returns the type of the Merging algorithm.
	Merge GetMergeAlgorithm();
	//! Returns the way, how to solve the unsymmetricity of SETTER algorithm (match(A, B) != match(B, A)).
	Mean GetMean();
	//! Returns the highliter string for Java GUI.
	/*!
		The Java GUI parses the output of the extern process and shows the lines, 
		which are highlited (which starts) by this string.
	*/
	string GetHighlight();
	//! Returns the multiplicator for allowed divergence from mean distance between GSSUs.
	double GetDiverMean();
	//! Returns the multiplicator for allowed divergence from parentGSSUs.
	double GetDiverParent();
	//! Returns wether the pdb printer is enabled or disabled
	bool GetScriptOutput();
};