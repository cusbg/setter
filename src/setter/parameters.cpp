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
#include "parameters.h"
#include "common.h"

#include <boost/algorithm/string.hpp>

using namespace std;

const string Parameters::THREADS = "Threads:";
const string Parameters::HALF_ROTATION = "HalfRotation:";
const string Parameters::ALGO = "Algo:";
const string Parameters::MERGE = "Merge:";
const string Parameters::MEAN = "Mean:";
const string Parameters::DIVER_MEAN = "DiverMean:";
const string Parameters::DIVER_PARENT = "DiverParent:";
const string Parameters::HIGHLIGHT = "Highlight:";
const string Parameters::SCRIPT = "Script:";

bool Parameters::mInitalized = false;
Parameters Parameters::mInstance;

Parameters::Parameters()
{
	if (mInitalized) {
		return;
	}

	mInitalized = true;

	mInstance.mHalfRotation = false;
	mInstance.mAlgo = Algo::STANDARD;
	mInstance.mMerge = Merge::ALTERNATING;
	mInstance.mMean = Mean::MIN;
	mInstance.mHighlight = "";
	mInstance.mDiverMean = 0.9;
	mInstance.mDiverParent = 2.0;
	mInstance.mScriptOutput = false;
	mThreads = -1;
}

inline double StringToDouble(const string &s)
{
	std::istringstream i(s);
	double x;
	if (!(i >> x))
		return 0;
	return x;
}

inline int StringToInt(const string &s)
{
	std::istringstream i(s);
	int x;
	if (!(i >> x))
		return 0;
	return x;
}

void Parameters::Init(string file)
{
	if (!FileExists(file)) {
		return;
	}

	mInitalized = true;

	ifstream ifs(file);

	mInstance.mHalfRotation = false;

	string line;
	while (!ifs.eof()) {
		getline(ifs, line);

		string option, param;
		option = line.substr(0, line.find(":") + 1);
		param = line.substr(line.find(":") + 1);
		boost::algorithm::to_upper(param);
		boost::algorithm::trim(param);

		if (option == THREADS) {
			mInstance.mThreads = StringToInt(param);
		}
		else if (option == HALF_ROTATION) {
			// no-op
		}
		else if (option == ALGO) {
			if (param == "STANDARD") {
				mInstance.mAlgo = Algo::STANDARD;
			}
			else if (param == "RECALC"){
				mInstance.mAlgo = Algo::RECALC;
			}
			else if (param == "ALL"){
				mInstance.mAlgo = Algo::ALL;
			}
		}
		else if (option == MERGE) {
			if (param == "ALTERNATING") {
				mInstance.mMerge = Merge::ALTERNATING;
			}
			else if (param == "PATTERN"){
				mInstance.mMerge = Merge::PATTERN;
			}
		}
		else if (option == MEAN) {
			if (param == "MIN") {
				mInstance.mMean = Mean::MIN;
			}
			else if (param == "ARITHM"){
				mInstance.mMean = Mean::ARITHM;
			}
			else if (param == "GEOM"){
				mInstance.mMean = Mean::GEOM;
			}
		}
		else if (option == DIVER_MEAN) {
			mInstance.mDiverMean = StringToDouble(param);

		}
		else if (option == DIVER_PARENT) {
			mInstance.mDiverParent = StringToDouble(param);
		}
		else if (option == HIGHLIGHT) {
			mInstance.mHighlight = param;
		}
		else if (option == SCRIPT) {
			if (param[0] == '1') {
				mInstance.mScriptOutput = true;
			}
			else {
				mInstance.mScriptOutput = false;
			}
		}
	}
}

Parameters &Parameters::GetInstance()
{
	return mInstance;
}

int Parameters::GetThreads()
{
	return mThreads;
}

bool Parameters::GetHalfRotation()
{
	return mHalfRotation;
}

Parameters::Algo Parameters::GetAlgorithm()
{
	return mAlgo;
}

Parameters::Merge Parameters::GetMergeAlgorithm()
{
	return mMerge;
}

Parameters::Mean Parameters::GetMean()
{
	return mMean;
}

string Parameters::GetHighlight()
{
	return mHighlight;
}

double Parameters::GetDiverMean()
{
	return mDiverMean;
}

double Parameters::GetDiverParent()
{
	return mDiverParent;
}

bool Parameters::GetScriptOutput()
{
	return mScriptOutput;
}