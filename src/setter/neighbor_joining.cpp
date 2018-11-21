#include <limits>
#include <math.h>
#include "neighbor_joining.h"
#include "make_aver_rna.h"
#include "algo.h"


NJResult::NJResult() {
	cout << "NJ Best" << endl;
	// no-op
}

NJResult::NJResult (
	cRNAStructure &rna, 
	vector<string> &map, 
	vector<double> &dist) 
	: 
	rna(rna),
	map(map),
	dist(dist)
{
	mean = 0.0;
	for (int i = 0; i < dist.size(); ++i) {
		mean += dist[i];
	}

	mean /= dist.size();

	deviation = 0.0;
	for (int i = 0; i < dist.size(); ++i) {
		deviation += (dist[i] - mean) * (dist[i] - mean);
	}

	deviation /= dist.size();
	deviation = sqrt(deviation);
}

NeighborJoining::NeighborJoining(Input &in) : mIn(in)
{
	cnt = 0;
	score = 0;
}

double NeighborJoining::GetScore()
{
	return score / cnt;
}

NeighborJoining *NeighborJoining::Factory(int id, Input &in)
{
	NeighborJoining *result;
	if (id == 0) {
		result = new NeighborJoiningImpl0(in);
	} else if (id == 1) {
		result = new NeighborJoiningImpl1(in);
	} else if (id == 2) {
		result = new NeighborJoiningImpl2(in);
	} /*else if (id == 4) {
		result = new NeighborJoiningImpl4(file);
	}*/ else {
		result = new NeighborJoiningImpl0(in);
	}

	return result;
}

NeighborJoining *NeighborJoining::Factory(string id, Input &in)
{
	NeighborJoining *result;
	if (id == "1") {
		result = new NeighborJoiningImpl0(in);
	}/* else if (id == "2") {
		result = new NeighborJoiningImpl2(file);
	} else if (id == "3") {
		result = new NeighborJoiningImpl3(file);
	} else if (id == "4") {
		result = new NeighborJoiningImpl4(file);
	}*/ else {
		result = new NeighborJoiningImpl0(in);
	}

	return result;
}

DistanceMatrix::DistanceMatrix() : mInput(NJInput())
{
	mSize = 0;
	matrix = NULL;
}

DistanceMatrix::DistanceMatrix(NJInput &input) : mInput(input)
{
	mSize = mInput.GetSize();
	matrix = new double * [mSize];
	for (int i = 0; i < mSize; ++i) {
		matrix[i] = new double [mSize];
	}

	for (int i = 0; i< mSize; ++i) {
		for (int j = 0; j < mSize; ++j) {
			matrix[i][j] = mInput.CalcScore(i, j);
		}
 	}
}

DistanceMatrix::DistanceMatrix(DistanceMatrix &ref) :
	mSize(ref.mSize),
	mInput(ref.mInput)
{
	matrix = new double * [mSize];
	for (int i = 0; i < mSize; ++i) {
		matrix[i] = new double [mSize];
	}

	for (int i = 0; i< mSize; ++i) {
		for (int j = 0; j < mSize; ++j) {
			matrix[i][j] = ref.matrix[i][j];
		}
 	}
}

DistanceMatrix::~DistanceMatrix()
{
	for (int i = 0; i < mSize; ++i) {
		delete [] matrix[i];
	}
	delete [] matrix;
}

double DistanceMatrix::operator()(int i, int j)
{
	if (i < 0 || j < 0 || i > mSize - 1|| j > mSize - 1 ) {
		return sResMatch::MAX;
	}
	return matrix[i][j];
}

void DistanceMatrix::CalcDist(
	pair<int, int> min, 
	int &cnt, 
	double &score)
{
	int f = min.first;
	int g = min.second;
	int r = mSize;

	if (r > 2) {
		double res1 = 0.0;
		for (int k = 0; k < r; ++k) {
			res1 += (matrix[f][k] - matrix[g][k]);
		}
		res1 /= 2 * (r - 2);
		res1 += 0.5 * matrix[f][g];

		double res2 = matrix[f][g] - res1; 
		
		cout << matrix[f][g] << endl;
		cout << mInput.GetRnaQIdx(f, g) << ":" << res1 << endl;
		cout << mInput.GetRnaQIdx(f, g) << ":" << res2 << endl;

	} else {
		double res = matrix[f][g] / 2.0;
		cout << matrix[f][g] << endl;
		cout << mInput.GetRnaQIdx(f, g) << ":" << res << endl;
		cout << mInput.GetRnaQIdx(f, g) << ":" << res << endl;
	}

	++cnt;
	score += matrix[f][g];
}

void DistanceMatrix::RemoveRNA(int idx)
{
	unsigned int size = mSize - 1;
	double **newMatrix;

	newMatrix = new double *[size];
	for (unsigned int i = 0; i < size; ++i) {
		newMatrix[i] = new double [size];
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
			newMatrix[ii][jj] = matrix[i][j];
			++jj;
		}
		++ii;
	}

	for (int i = 0; i < mSize; ++i) {
		delete [] matrix[i];
	}
	delete [] matrix;
	
	
	mSize = size;
	matrix = newMatrix;
}

void DistanceMatrix::AddRna(
	DistanceMatrix &copyDM, 
	pair<int, int> min)
{
	unsigned int size = mSize + 1;
	double **newMatrix;

	newMatrix = new double *[size];
	for (unsigned int i = 0; i < size; ++i) {
		newMatrix[i] = new double [size];
	}

	for (unsigned int i = 0; i < mSize; ++i) {
		for (unsigned int j = 0; j < mSize; ++j) {
			newMatrix[i][j] = matrix[i][j];
		}
	}

	int f = min.first;
	int g = min.second;
	for (int k = 0; k < mSize; ++k) {
		double dist = (copyDM(f, k) + copyDM(g, k) - copyDM(f, g)) / 2.0;
		newMatrix[k][mSize] = newMatrix[mSize][k] = abs(dist);
	}
	newMatrix[mSize][mSize] = 0.0;

	for (int i = 0; i < mSize; ++i) {
		delete [] matrix[i];
	}
	delete [] matrix;
	
	
	mSize = size;
	matrix = newMatrix;
}

QMatrix::QMatrix() : mDistMatrix(DistanceMatrix())
{
	mSize = 0;
	matrix = NULL;
}

QMatrix::QMatrix(DistanceMatrix &distMatrix) : mDistMatrix(distMatrix)
{
	mSize = distMatrix.mSize;
	matrix = new double * [mSize];
	for (int i = 0; i < mSize; ++i) {
		matrix[i] = new double [mSize];
	}
	Calculate();
}

QMatrix::~QMatrix()
{
	for (int i = 0; i < mSize; ++i) {
		delete [] matrix[i];
	}
	delete [] matrix;
}

void QMatrix::Calculate()
{
	double min = numeric_limits<double>::max();

	//cout << "min dist matrix" << endl;
	for (int i = 0; i < mSize; ++i) {
		for (int j = 0; j < mSize; ++j) {
			matrix[i][j] = Calculate(i, j);
			//cout << matrix[i][j] << ", ";

			if (i != j && matrix[i][j] < min) {
				min = matrix[i][j];
				mMinPair = make_pair(i, j);
			}
		}
		//cout << endl;
	}
}

double QMatrix::Calculate(int i, int j)
{
	if (i == j) {
		return 0.0;
	}
	int r = mSize;

	double result = (r - 2) * mDistMatrix(i, j);
	for (int k = 0; k < r; ++k) {
		result -= mDistMatrix(i, k);
		result -= mDistMatrix(j, k);
	}

	return result;
}

pair<int, int> QMatrix::GetMinPair()
{
	return mMinPair;
}