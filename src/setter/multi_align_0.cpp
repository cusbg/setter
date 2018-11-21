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

MultiAlignImpl0::MultiAlignImpl0(Input &in) :
MultiAlign(in)
{
	cout << "Using the standard NJ algorithm..." << endl;
	// no-op
}

MAResult MultiAlignImpl0::GetResult()
{
	Parameters &param = Parameters::GetInstance();
	NJInput input(mIn);
	DistanceMatrix dm(input);

	cout << "NJ mergiging ( " << input.mSize << " structures) to merge: ";
	while (input.mSize != 2) {
		/*
		 * Iteratively identify the pair of structures with minimum distance,
		 * merge them into one structure, remove respective records from the
		 * all to all distance matrix, add the merged structure and modify
		 * the distances in the distance matrix (as the NJ algorithm prescribes).
		 */
		DistanceMatrix copyDM(dm);
		QMatrix qm(dm);
		pair<int, int> min = qm.GetMinPair();
		MatchData minData;
		input.GetData(min.first, min.second, minData);

		string rnaIdx = "(" + minData.sRnaQIdx + "-" + minData.sRnaDBIdx + ")";
		cout << rnaIdx << ".";


		try {
			FormAvg ma(input.mMean, minData, cnt, score);
			cRNAStructure *rna = new cRNAStructure(ma.getResult());

			input.RemoveRNA(min.first);
			dm.RemoveRNA(min.first);
			if (min.first < min.second) {
				--min.second;
			}
			input.RemoveRNA(min.second);
			dm.RemoveRNA(min.second);
			input.AddRNA(rna, rnaIdx, false);
			dm.AddRna(copyDM, min);
		}
		catch (exception e) {
			cout << "Exception while merging: " << e.what() << endl;
		}
	}

	MatchData minData;
	input.GetData(0, 1, minData);

	string rnaIdx = "(" + minData.sRnaQIdx + "-" + minData.sRnaDBIdx + ")";
	cout << rnaIdx << ".";

	FormAvg ma(input.mMean, minData, cnt, score);
	cRNAStructure *rna = new cRNAStructure(ma.getResult());

	string sequence = rnaIdx;

	cout << /*Parameters::GetHighlight() <<*/ endl << "Input-Average distances: " << endl;

	vector<double> dist(input.mInputSize);
	vector<ResMatch> results(input.mInputSize);

	int start = clock();
	//Align the input RNA structures with the new average structure
	ParallelResult parallelResult(input, *rna, results, dist);
	tbb::parallel_for(
		tbb::blocked_range<size_t>(0, input.mInputRnaVect.size()), parallelResult);
	int end = clock();
	cout << "Runtime: " << (float)(end - start) / CLOCKS_PER_SEC << endl;

	for (int i = 0; i < input.mInputRnaVect.size(); ++i) {
		cout << mIn.mOriginalMap[i] << ": " << dist[i] << endl;
	}

	//Put all the details about the merging (average structure, distances to the input structures, ...)
	MAResult result(*rna, input.mOriginalMap, dist, results);

	cout << "Mean     : " << result.mean << endl;
	cout << "Deviation: " << result.deviation << endl;

	return result;
}
