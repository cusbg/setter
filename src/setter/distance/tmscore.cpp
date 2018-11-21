/*
 Copyright (c) 2013 David Hoksza

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
#include <cmath>
#include <algorithm>
#include "rmsd.hpp"
#include "protein.hpp"
#include "tmscore.hpp"

using namespace std;



/**
 * The context shared by TM-score functions.
 */
struct context
{
	double _d0;
	double _d0_search;
	double _score_max;  // the actually best nn-TM-score
	cRMSD3DCoord _tx_max;     // the actually best query protein translation
	cRMSD3DCoord _ty_max;     // the actually best database protein translation
	cRMSDMatrix _u_max;      // the actually best database protein rotation (after the translation)
	cRMSDAlign  _align_max;  // the actually best alignment
	cRMSDAlign  _cut_max;    // the actually best cut (of the actually best alignment)
};



/**
 * Shortcuts
 */
#define d0        (c._d0)
#define d0_search (c._d0_search)
#define score_max (c._score_max)
#define u_max     (c._u_max)
#define tx_max    (c._tx_max)
#define ty_max    (c._ty_max)
#define align_max (c._align_max)
#define cut_max   (c._cut_max)



/**
 * The core of the TM-Score algorithm.
 */
void
tmscore_core(context &c, cRMSDStruct *qprot, cRMSDStruct *dprot, cRMSDAlign &align, cRMSDAlign &cut)
{
	/** parameter */
	static int max_iteration = 20;      // maximum number of the refinement of the cut

	bool stable = false;                // flag indicating whether the refinement (cut) becomes stable (i.e. whether a new cut is equal to previous one)


	/**
	 * performs the refinement, until the cut becomes stable or the maximum number of the iteration is reached.
	 */
	for (int it = 0; it < max_iteration && !stable; it++)
	{
		/**
		 * calculates RMSD transformation for the actual cut.
		 */
		cRMSDMatrix u;
		cRMSD3DCoord tx;
		cRMSD3DCoord ty;
		rmsd(qprot, dprot, cut, u, tx, ty);


		/**
		 * calculates nn-TM-score for the RMSD transformation and creates the new cut (on the place of the actual one).
		 */
		double d_tmp = it == 0 ? d0_search - 1 : d0_search + 1; // a cut off distance
		double score_sum;                                       // nn-TM-score
		int old_cut_len = cut.length();                         // the length of the previous cut

		/// sets the length of the cut to zero.
		cut.setlength(0);

		do
		{
			/**
			 * calculates nn-TM-score for the actual cut and creates the new cut.
			 */
			{
				// initial values
				stable = true;
				score_sum = 0;

				/// for each aligned residue-pair ...
				for (int i = 0; i < align.length(); i++)
				{
					/// the distance vektor of the residue-pair (after performing the RMSD transformation)
					cRMSD3DCoord diff = (qprot->coord[align[i].x] - tx) - u * (dprot->coord[align[i].y] - ty);

					/// the distance of the residue-pair (after performing the RMSD transformation)
#ifdef NOWEIGHT
					double dis = sqrt(diff * diff);
#else
					double dis = sqrt(diff * diff) * align[i].w;
#endif /*NOWEIGHT*/
					//cout << " align[i].w:" << align[i].w;

					/// add the score of the resid pair
					score_sum += 1 / (1 + ((dis / d0) * (dis / d0)));

					/// the new cut includes all residue-pairs, that have distance less than d_tmp
					if (dis < d_tmp)
					{
						int len = cut.length();
						cut.setlength(len + 1);

						/** checks for the difference between the constructed cut and the previous cut (on the actual position).
						 *      cut[len]  ... previous value
						 *      align[i]  ... new value
						 */
						if (len >= old_cut_len || cut[len] != align[i])
							stable = false;

						// adds residue-pair into the cut.
						cut[len] = align[i];
					}
				}
			}

			/** checks for the difference between the constructed cut and the previous cut */
			if (cut.length() != old_cut_len)
				stable = false;

			d_tmp = d_tmp + 0.5;
		} while (cut.length() < 3 && align.length() > 3); // if the new cut is too short, it is recalculated for the bigger cut off distance (d_tmp).


		/**
		 * checks, whether the obtain nn-TM-score is better than the actual best nn-TM-score
		 */
		if (score_max <=/*<*/ score_sum)
		{
			score_max = score_sum;

			u_max = u;
			tx_max = tx;
			ty_max = ty;

			align_max = align;
			cut_max = cut.copy();
		}
	}
}



/**
 * The TM-Score algorithm.
 */
void
tmscore(context &c, tmmode_t mode, cRMSDStruct *qprot, cRMSDStruct *dprot, cRMSDAlign &align)
{
	if (align.length() == 0)
		return;

	cRMSDAlign cut(std::min(qprot->length, dprot->length));      // cut (yet unfilled) with the maximal possible length

	/**
	 * Generate initials cuts (according to the specified mode)
	 */
	switch (mode)
	{
	case TM_FULL:
	{
		/** parameter */
		static int iterations_max = 6;                  // maximum number of different initial lengths of cuts

		int cut_length_min = min(align.length(), 4);    // minimal possible length of initial cuts
		int cut_length = align.length();            // initial length of initial cuts

		for (int i = 0; i < iterations_max; i++)         // each iteration works with the different length of cuts
		{
			// for each possible start positions ...
			for (int iL = 0; iL <= align.length() - cut_length; iL++)
			{
				// create a cut and refine it 
				cut.set(align, iL, iL + cut_length);
				tmscore_core(c, qprot, dprot, align, cut);
			}

			// exit loop, if the minimal length of initial cuts is reached.
			if (cut_length == cut_length_min)
				break;

			// set the new length of initial cuts
			cut_length = max(cut_length >> 1, cut_length_min);
		}

		break;
	}

	case TM_FAST_SSE:
	{
		int start = 0;      // a possible start possition of the cut

		for (int i = 1; i <= align.length(); i++)
		{                                                           // close the cut ...
			if (i == align.length() ||                                   // at the end of the alignment
				align[i].x != align[i - 1].x + 1 ||                        // on a gap in the query protein
				align[i].y != align[i - 1].y + 1 ||                        // on a gap in the database protein
				qprot->sse[align[i].x] != qprot->sse[align[i - 1].x] ||    // on a change in secondary structure of the query protein                          
				dprot->sse[align[i].y] != dprot->sse[align[i - 1].y] ||    // on a change in secondary structure of the database protein
				qprot->sse[align[i].x] != dprot->sse[align[i].y] ||      // on a difference in secondary structures
				qprot->sse[align[i].x] == 'T')                           // on a turn in secondary structure(s)
			{
				if (i - start > 4)
				{
					cut.set(align, start, i);
					tmscore_core(c, qprot, dprot, align, cut);
				}

				start = i;
			}
		}

		// at the end use the full alignment
		cut.set(align, 0, align.length());
		tmscore_core(c, qprot, dprot, align, cut);
		break;
	}

	case TM_FAST_ESSE:
	{
		int start = 0;
		int pre_start[2] = { -1, -1 };
		int pre_stop[2] = { -1, -1 };

		for (int i = 1; i <= align.length(); i++)
		{
			if (i == align.length() ||
				align[i].x != align[i - 1].x + 1 ||
				align[i].y != align[i - 1].y + 1 ||
				qprot->sse[align[i].x] != qprot->sse[align[i - 1].x] ||
				dprot->sse[align[i].y] != dprot->sse[align[i - 1].y] ||
				qprot->sse[align[i].x] != dprot->sse[align[i].y] ||
				qprot->sse[align[i].x] == 'T')
			{
				if (i - start > 4)
				{
					cut.set(align, start, i);
					tmscore_core(c, qprot, dprot, align, cut);

					if (pre_start[1] != -1)
					{
						cut.set(align, pre_start[1], pre_stop[1]);
						cut.add(align, start, i);
						tmscore_core(c, qprot, dprot, align, cut);
					}

					if (pre_start[0] != -1)
					{
						cut.set(align, pre_start[0], pre_stop[0]);
						cut.add(align, pre_start[1], pre_stop[1]);
						cut.add(align, start, i);
						tmscore_core(c, qprot, dprot, align, cut);
					}

					pre_start[0] = pre_start[1];
					pre_stop[0] = pre_stop[1];
					pre_start[1] = start;
					pre_stop[1] = i;
				}

				start = i;
			}
		}

		cut.set(align, 0, align.length());
		tmscore_core(c, qprot, dprot, align, cut);
		break;
	}

	case TM_FAST_eSSE:
	{
		int start = 0;
		int pre_start = -1;
		int pre_stop = -1;

		for (int i = 1; i <= align.length(); i++)
		{
			if (i == align.length() ||
				align[i].x != align[i - 1].x + 1 ||
				align[i].y != align[i - 1].y + 1 ||
				qprot->sse[align[i].x] != qprot->sse[align[i - 1].x] ||
				dprot->sse[align[i].y] != dprot->sse[align[i - 1].y] ||
				qprot->sse[align[i].x] != dprot->sse[align[i].y] ||
				qprot->sse[align[i].x] == 'T')
			{
				if (i - start > 4)
				{
					cut.set(align, start, i);
					tmscore_core(c, qprot, dprot, align, cut);

					if (pre_start != -1)
					{
						cut.set(align, pre_start, pre_stop);
						cut.add(align, start, i);
						tmscore_core(c, qprot, dprot, align, cut);
					}

					pre_start = start;
					pre_stop = i;
				}

				start = i;
			}
		}

		cut.set(align, 0, align.length());
		tmscore_core(c, qprot, dprot, align, cut);
		break;
	}

	case TM_FAST:
		// use only the whole alignment as the initial cut
		cut.set(align, 0, align.length());
		tmscore_core(c, qprot, dprot, align, cut);
		break;
	}
}



/**
 * The iterative TM-Score algorithm.
 */
double
itmscore(int dpiterations, int dpwidth, tmmode_t mode, double normParam, cRMSDStruct *qprot, cRMSDStruct *dprot, cRMSDAlign &align, std::ostream *stream)
{
	/**
	 * Create context
	 */
	context c;
	d0 = normParam > 21 ? 1.24 * pow(normParam - 15, 1 / 3.0) - 1.8 : 0.5;
	d0_search = max(min(d0, 8.0), 4.5);
	score_max = -1;


	/**
	 * TM-score for the original alignment
	 */
	tmscore(c, mode, qprot, dprot, align);


	/**
	 * DP iterations
	 */
	for (int i = 0; i < dpiterations; i++)
	{
		cRMSDStruct prot(dprot, u_max, tx_max - u_max * ty_max);

		cRMSDAlign dpalign = align.dp(qprot, &prot, dpwidth, d0);

		if (dpalign == align)
			break;

		align = dpalign;

		tmscore(c, mode, qprot, dprot, align);
	}


	/**
	 * VMD
	 */
	if (stream)
		generate_vmd(*stream, qprot, dprot, align, align_max, cut_max, u_max, tx_max, ty_max);


	/**
	 * Returns
	 */
	align = align_max;
	return 1 - score_max / normParam;
}
