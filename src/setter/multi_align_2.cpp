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


MultiAlignImpl2::MultiAlignImpl2(Input &in) :
MultiAlign(in)
{
	cout << "NJ, recalculate if necessary" << endl;
	// no-op
}

MAResult MultiAlignImpl2::GetResult()
{
	NJInput input(mIn);

	cRNAStructure *rna;
	string rnaIdx;
	MatchData data;
	while (true) {
		vector<cRNAStructure *> vect;
		vector<string> map;
		while (input.mSize > 1) {
			DistanceMatrix dm(input);
			QMatrix qm(dm);
			pair<int, int> min = qm.GetMinPair();
			MatchData minData;
			input.GetData(min.first, min.second, minData);
			data = minData;

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
				vect.push_back(rna);
				map.push_back(rnaIdx);
			}
			catch (exception e) {
				cout << e.what() << endl;
			}
		}

		if (vect.size() == 1 && input.mSize == 0) {
			rna = vect[0];
			rnaIdx = map[0];
			break;
		}

		for (int i = 0; i < vect.size(); ++i) {
			input.AddRNA(vect[i], map[i]);
		}
	}
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

	return result;
}
