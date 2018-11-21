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
#include <fstream>
#include <time.h>
#include <math.h>

#include "tbb/task_scheduler_init.h"

#include "input.h"
#include "algo.h"
#include "parameters.h"

#include <boost/algorithm/string.hpp>

using namespace std;

Input::Input()
{
	mInputSize = 0;
}

Input::Input(string file)
{
	Parser::ParseInput(file, mInputRnaVect, mOriginalMap);
	mInputSize = mInputRnaVect.size();
}

Input::Input(Input &in)
{
	mInputSize = in.mInputSize;
	for (int i = 0; i < mInputSize; ++i) {
		mInputRnaVect.push_back(new cRNAStructure(*in.mInputRnaVect[i]));
	}
	mOriginalMap = in.mOriginalMap;
}

Input& Input::operator = (const Input &in)
{
	if (this != &in) // protect against invalid self-assignment
    {
		mInputSize = in.mInputSize;
		for (int i = 0; i < mInputSize; ++i) {
			mInputRnaVect.push_back(new cRNAStructure(*in.mInputRnaVect[i]));
		}
		mOriginalMap = in.mOriginalMap;
	}

	return *this;
}

Input::Input(vector<cRNAStructure *> &inputRnaVect)
{
	for (int i = 0; i < inputRnaVect.size(); ++i) {
		mInputRnaVect.push_back(new cRNAStructure(*inputRnaVect[i]));
	}

	mInputSize = mInputRnaVect.size();
}

Input::Input(vector<InputRNA> &inputRnaDefVect)
{
	mInputRnaVect.resize(inputRnaDefVect.size());
	mOriginalMap.resize(inputRnaDefVect.size());

	ParallelRead parallelRead(inputRnaDefVect, mInputRnaVect, mOriginalMap);
	//parallelRead.serial(inputRnaDefVect.size());
	tbb::parallel_for(
		tbb::blocked_range<size_t>(0, inputRnaDefVect.size()), parallelRead);

	mInputSize = mInputRnaVect.size();
}

Input::~Input()
{
	for (int i = 0; i < mInputSize; i++) {
		delete mInputRnaVect[i];
	}
}

void Input::RemoveRNA(int pos)
{
	mInputRnaVect.erase(mInputRnaVect.begin() + pos);
	mOriginalMap.erase(mOriginalMap.begin() + pos);
	--mInputSize;
}

NJInput::NJInput(): Input()
{
	mSize = mInputSize;
	mData = NULL;
}

NJInput::NJInput(Input &in): Input(in)
{
	mSize = mInputSize;
	for (int i = 0; i < mInputSize; ++i) {
		mRnaVect.push_back(new cRNAStructure(*mInputRnaVect[i]));
	}

	stringstream sstream;
	mData = new MatchData **[mSize];
	for (unsigned int i = 0; i < mSize; ++i) {
		sstream.str(string());
		sstream.clear();
		sstream << i;
		mStringMap.push_back(sstream.str());
		mData[i] = new MatchData *[mSize];
	}
	cout << "Computing distances matrix for NJ..." << endl;
#if PARALLEL_INPUT
	
	int start = clock();

	ParallelInput parallelInput(mRnaVect, &mData);

	tbb::parallel_for(
		tbb::blocked_range2d<size_t>(0, mSize, 0, mSize), parallelInput);

	int end = clock();

	cout << "Distance matrix computation runtime: " << (float)(end - start)/CLOCKS_PER_SEC << endl; 
#else
	for (unsigned int i = 0; i < mSize; ++i) {
		for (unsigned int j = 0; j < mSize; ++j) {
			mData[i][j] = new MatchData(mRnaVect, i, j);
		}
	}
#endif

	//Print();
	//Clean();
	//Print();
	SortData();
	CalcMean();
}

double round(double number)
{
    return number < 0.0 ? ceil(number - 0.5) : floor(number + 0.5);
}

void NJInput::CalcMean() {
	
	
	vector<double> distances;
	
	for (int i = 0; i < mSize; ++i) {
		for (int j = 0; j < mSize; ++j) {
			if (i == j) continue;

			MatchData &md = *mData[i][j];

			for (int qIdx = 0; qIdx < md.mCntHpQ; ++qIdx) {
				for (int dbIdx = 0; dbIdx < md.mCntHpDB; ++dbIdx) {
					if (md(qIdx, dbIdx).score != sResMatch::MAX) {
						distances.push_back(md(qIdx, dbIdx).score);
					}
				}
			}
		}
	}


	sort(distances.begin(), distances.end());
	
	double p = 0.8; //seems ideal
	
	int begin =  round(distances.size() * (0.5 - p / 2.0));
	int end	  =  round(distances.size() * (0.5 + p / 2.0));
	double sum = 0.0;
	for (int i = begin; i < end; ++i) {
		sum += distances[i];
	}
	mMean = sum / (end - begin);
}

NJInput::NJInput(NJInput &input)
{
	mInputSize = input.mInputSize;
	for (int i = 0; i < input.mInputRnaVect.size(); ++i) {
		mInputRnaVect.push_back(new cRNAStructure(*input.mInputRnaVect[i]));
	}
	for (int i = 0; i < input.mRnaVect.size(); ++i) {
		mRnaVect.push_back(new cRNAStructure(*input.mRnaVect[i]));
	}

	mStringMap = input.mStringMap;
	mSize = input.mSize;

	mData = new MatchData **[mSize];
	for (unsigned int i = 0; i < mSize; ++i) {
		mData[i] = new MatchData *[mSize];
	}

	for (unsigned int i = 0; i < mSize; ++i) {
		for (unsigned int j = 0; j < mSize; ++j) {
			mData[i][j] = new MatchData(*input.mData[i][j]); // modified
		}
	}
}

NJInput::~NJInput()
{
	for (int i = 0; i < mRnaVect.size(); ++i) {
		delete mRnaVect[i];
	}

	DeleteData();
}

int NJInput::GetSize()
{
	return mSize;
}

MatchData NJInput::operator()(int i, int j)
{
	if (i < 0 || j < 0 || i > mSize - 1|| j > mSize - 1 ) {
		MatchData md;
		return md;
	}
	MatchData result;

	if (mSortedList[i] < mSortedList[j]) {
		result = *mData[i][j];
	} else {
		result = *mData[j][i];
	}

	result.GetResult().score = CalcScore(i, j);

	return result;
}

void NJInput::GetData(int i, int j, MatchData &data)
{
	//data = *mData[i][j];

	if (mSortedList[i] < mSortedList[j]) {
		data = *mData[i][j];
	} else {
		data = *mData[j][i];
	}

	data.GetResult().score = CalcScore(i, j);

}

string NJInput::GetRnaQIdx(int i, int j)
{
	return mData[i][j]->sRnaQIdx;
}
string NJInput::GetRnaDBIdx(int i, int j)
{
	return mData[i][j]->sRnaDBIdx;
}

bool ReverseCompareFirst(const pair<int, double> x1, const pair<int, double> x2 ) 
{
	return x1.first > x2.first;
}

bool CompareSecond(const pair<int, double> x1, const pair<int, double> x2 ) 
{
	return x1.second < x2.second;
}

void NJInput::SortData()
{
	typedef pair<int, double> map;
	mSortedList = vector<int>(mSize);
	vector<map> distances(mSize);
	for (int i = 0; i < mSize; ++i) {
		double sum = 0.0;
		for (int j = 0; j < mSize; ++j) {
			if (i == j) continue;
			sum += mData[i][j]->GetResult().score;
			sum += mData[j][i]->GetResult().score;
		}
		distances[i] = make_pair(i, sum);
	}

	//sort(distances.begin(), distances.end(), [](map a, map b){ return a.second < b.second; });
	sort(distances.begin(), distances.end(), CompareSecond);

	for (int i = 0; i < mSize; ++i) {
		mSortedList[distances[i].first] = i;
	}
}

void NJInput::Print()
{
	cout << setiosflags(ios::fixed) << setprecision(4);
	cout << endl;
	for (int i = 0; i < mSize; ++i) {
		//cout << i << ":" << endl;
		for (int j = 0; j < mSize; ++j) {
			cout << mData[i][j]->GetResult().score << ";";
		}
		cout << '\b' << ' ' << endl;
	}
	cout << endl;
	cout << resetiosflags(ios::fixed);
}


double NJInput::CalcScore(int i, int j)
{
	return mData[i][j]->GetResult().score;
}

void NJInput::RemoveRNA(int idx)
{
	unsigned int size = mSize - 1;
	MatchData ***data;

	data = new MatchData **[size];
	for (unsigned int i = 0; i < size; ++i) {
		data[i] = new MatchData *[size];
	}

	unsigned int ii = 0;
	for (unsigned int i = 0; i < mSize; ++i) {
		if (i == idx) {
			continue;
		}
		unsigned int jj = 0;
		for (unsigned int j = 0; j < mSize; ++j) {
			if (j == idx) {
				continue;
			}
			data[ii][jj] = new MatchData(*mData[i][j]);
			++jj;
		}
		++ii;
	}

	delete mRnaVect[idx];
	mRnaVect.erase(mRnaVect.begin() + idx);
	mStringMap.erase(mStringMap.begin() + idx);
	DeleteData();
	
	mSize = size;
	mData = data;
}

void NJInput::AddRNA(
	cRNAStructure *rna, 
	string rnaIdx,
	bool init)
{
	unsigned int size = mSize + 1;
	MatchData ***data;

	data = new MatchData **[size];
	for (unsigned int i = 0; i < size; ++i) {
		data[i] = new MatchData *[size];
	}

	for (unsigned int i = 0; i < mSize; ++i) {
		for (unsigned int j = 0; j < mSize; ++j) {
			data[i][j] = new MatchData(*mData[i][j]);
		}
	}

#if PARALLEL_INPUT
	int start = clock();

	ParallelAddRNA parallelAddRna(mSize, *rna, rnaIdx, mRnaVect, mStringMap, &data, init);
	tbb::parallel_for(
		tbb::blocked_range<size_t>(0, mSize), parallelAddRna);

	int end = clock();
	//cout << "add rna: " << end - start << endl;

#else
	for (int i = 0; i < mSize; ++i) {
		MatchData m1(mRnaVect[i], rna, mMap[i], rnaIdx, init);
		MatchData m2(rna, mRnaVect[i], rnaIdx, mMap[i], init);
		if (init) {
			if (m1.GetResult().score < m2.GetResult().score) {
				data[i][mSize] = &m1;
			} else {
				data[i][mSize] = &m2;
			}
		} else {
			data[i][mSize] = &m1;
			data[mSize][i] = &m2;
		}
	}
#endif
	
	data[mSize][mSize] = new MatchData(rna, rna, rnaIdx, rnaIdx, false);
	data[mSize][mSize]->GetResult().score = 0.0;

	mRnaVect.push_back(rna);
	mStringMap.push_back(rnaIdx);
	
	DeleteData();
	mSize = size;
	mData = data;

	SortData();
}

void NJInput::DeleteData()
{
	if (mSize == 0) {
		return;
	}

	for (int i = 0; i < mSize; ++i) {
		for (int j = 0; j < mSize; ++j) {
			delete mData[i][j];
		}
		delete [] mData[i];
	}
	delete [] mData;
}


Parser::Parser()
{
	// no-op
}

void Parser::ParseInput(
	string file, 
	vector<cRNAStructure *> &vect,
	vector<string> &map)
{
	if (!FileExists(file)) {
		throw "Error: the input file " + file + " could not be opened.";
	}

	ifstream ifs(file);

	vector<InputRNA> inputRNAs;
	InputRNA inRNA;

	cout << "Parsing the input structures..." << endl;

	string line, pdb, x3dna, chain, comment;
	int cntTried = 0;
	while (!ifs.eof()) {
		line = pdb = x3dna = chain = comment = "";
		for (int i = 0; ; ++i) {
			getline(ifs, line);
			boost::algorithm::trim(line);
			if (line == "") {
				break;
			};
			switch (i) {
			case 0: 
				inRNA.pdb = line;
				cntTried++; 
				break;
			case 1:
				inRNA.x3dna = line;
				break;
			case 2: 
				inRNA.chain = line;
				break;
			case 3:
				inRNA.comment = line;
				break;
			}
		}
		inputRNAs.push_back(inRNA);

	}

	int start = clock();

	vect.resize(inputRNAs.size());
	map.resize(inputRNAs.size());

	ParallelRead parallelRead(inputRNAs, vect, map);
	tbb::parallel_for(
		tbb::blocked_range<size_t>(0, inputRNAs.size()), parallelRead);

	int end = clock();
	cout << "Parse time: " << (float)(end - start)/CLOCKS_PER_SEC << "s" << endl;
}

void Parser::ParseInput(string file, cRNAStructure &rna)
{
	ifstream ifs(file);
	if (ifs) {
		string line, pdb, x3dna, chain, comment;
	
		line = pdb = x3dna = chain = comment = "";
		for (int i = 0; ; ++i) {
			getline(ifs, line);
			if (line == "") {
				break;
			};
			switch (i) {
			case 0: 
				pdb = line; 
				break;
			case 1:
				x3dna = line;
				break;
			case 2: 
				chain = line;
				break;
			case 3:
				comment = line;
				break;
			}
		}
		try {
			rna = cRNAStructure(pdb, x3dna, chain);
			cout << "query: " << pdb << ", " << comment << " was read" << endl;
		} catch (exception &e) {
			cout << "problem with " << pdb << ", " << comment << endl;
		}
	}
}

ParallelRead::ParallelRead(
		vector<InputRNA> &inputRNAs,
		vector<cRNAStructure *> &vect,
		vector<string> &map) :
		mInputRNAs(inputRNAs),
		mVect(vect),
		mMap(map)
{
	// no-op
}


void ParallelRead::operator()(const tbb::blocked_range<size_t> &r) const
{
	for( size_t i = r.begin(); i != r.end(); ++i) {
		cout << "Parsing " << mInputRNAs[i].pdb << endl;
		cRNAStructure *rna = new cRNAStructure(
			mInputRNAs[i].pdb, mInputRNAs[i].x3dna, mInputRNAs[i].chain, mInputRNAs[i].exclusions, mInputRNAs[i].ranges);
		cout << "Finished parsing " << mInputRNAs[i].pdb << endl;
		mVect[i] = rna;
		mMap[i] = GetPdbId(*rna) + "-" + mInputRNAs[i].chain;
	}
}

void ParallelRead::serial(int size) const
{
	for( size_t i = 0; i < size; ++i) {
		cout << "Parsing " << mInputRNAs[i].pdb << endl;
		cRNAStructure *rna = new cRNAStructure(
			mInputRNAs[i].pdb, mInputRNAs[i].x3dna, mInputRNAs[i].chain, mInputRNAs[i].exclusions, mInputRNAs[i].ranges);
		cout << "Finished parsing " << mInputRNAs[i].pdb << endl;
		mVect[i] = rna;
		mMap[i] = GetPdbId(*rna) + "-" + mInputRNAs[i].chain;
	}
}