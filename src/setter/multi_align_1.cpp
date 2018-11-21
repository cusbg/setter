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
#include <limits>
#include "multi_align.h"
#include "make_aver_rna.h"
#include "algo.h"
#include "parameters.h"

MultiAlignImpl1::MultiAlignImpl1(Input &in) :
MultiAlign(in)
{
	cout << "NJ - recalculate" << endl;
}

MAResult MultiAlignImpl1::GetResult()
{
	NJInput input(mIn);

	while (input.mSize != 2) {
		DistanceMatrix dm(input);
		QMatrix qm(dm);
		pair<int, int> min = qm.GetMinPair();
		MatchData minData;
		input.GetData(min.first, min.second, minData);

		string rnaIdx = "(" + minData.sRnaQIdx + "-" + minData.sRnaDBIdx + ")";
		cout << "Matching: " << rnaIdx << endl;
		try {
			FormAvg ma(input.mMean, minData, cnt, score);
			cRNAStructure *rna = new cRNAStructure(ma.getResult());


			input.RemoveRNA(min.first);
			if (min.first < min.second) {
				--min.second;
			}
			input.RemoveRNA(min.second);
			input.AddRNA(rna, rnaIdx);
		}
		catch (exception e) {
			cout << e.what() << endl;
		}
	}
	MatchData minData;
	input.GetData(0, 1, minData);

	string rnaIdx = "(" + minData.sRnaQIdx + "-" + minData.sRnaDBIdx + ")";
	cout << "Matching: " << rnaIdx << endl;
	FormAvg ma(input.mMean, minData, cnt, score);
	cRNAStructure *rna = new cRNAStructure(ma.getResult());

	string sequence = rnaIdx;

	cout << /*Parameters::GetHighlight() <<*/ "the sequence: " << sequence << endl;
	cout << /*Parameters::GetHighlight() <<*/ "score: " << GetScore() << endl;
	cout << /*Parameters::GetHighlight() <<*/ "match with originals: " << endl;

	vector<double> dist(input.mInputSize);
	vector<ResMatch> results(input.mInputSize);

	int start = clock();
	ParallelResult parallelResult(input, *rna, results, dist);
	tbb::parallel_for(
		tbb::blocked_range<size_t>(0, input.mInputRnaVect.size()), parallelResult);
	int end = clock();
	cout << "compute time " << end - start << endl;

	for (int i = 0; i < input.mInputRnaVect.size(); ++i) {
		cout << mIn.mOriginalMap[i] << ": " << dist[i] << endl;
	}

	MAResult result(*rna, input.mOriginalMap, dist, results);

	cout << "mean     : " << result.mean << endl;
	cout << "deviation: " << result.deviation << endl;

	return result;
}
