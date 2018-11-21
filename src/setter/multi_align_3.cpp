//#include <limits>
//#include "neighbor_joining.h"
//#include "make_aver_rna.h"
//#include "algo.h"
//
//
//MultiAlignImpl4::MultiAlignImpl4(string file) : 
//	MultiAlign(file)
//{
//	// no-op
//}
//
//	// will be imporoved...
//cRNAStructure MultiAlignImpl4::GetResult()
//{
//	Input input(mFile);
//	DistanceMatrix dm(input);
//
//	while(input.mSize != 2) {
//		DistanceMatrix copyDM(dm);
//		QMatrix qm(dm);
//		pair<int, int> min = qm.GetMinPair();
//
//		MatchData &minData = *input.mData[min.first][min.second];
//
//		string rnaIdx = "(" + minData.sRnaQIdx + "-" + minData.sRnaDBIdx + ")";
//		cout << "Matching: " << rnaIdx << endl;
//		dm.CalcDist(min, cnt, score);
//		
//		input.RemoveRNA(min.first);
//		dm.RemoveRNA(min.first);
//		if (min.first < min.second) {
//			--min.second;
//		}
//		input.RemoveRNA(min.second);
//		dm.RemoveRNA(min.second);
//		input.AddRNA(new cRNAStructure(), rnaIdx, false);
//		dm.AddRna(copyDM, min);
//	}
//
//	MatchData &minData = *input.mData[0][1];
//	
//	string rnaIdx = "(" + minData.sRnaQIdx + "-" + minData.sRnaDBIdx + ")";
//	cout << "Matching: " << rnaIdx << endl;
//	dm.CalcDist(make_pair(0, 1), cnt, score);
//	cRNAStructure *rna = new cRNAStructure();
//
//	cout << "the sequence: " << rnaIdx << endl;
//	cout << "score: " << GetScore() << endl;
//	/*cout << "match with originals: " << endl;
//	
//	for (int i = 0; i < input.mInputRnaVect.size(); ++i) {
//		cout << i << ": ";
//		cout << Match(*rna, *input.mInputRnaVect[i], GetDefaultParams()).score << endl;
//	}*/
//
//	return *rna;
//}