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

#include <iostream>
#include <cmath>
#include <algorithm>
#include "rmsd.hpp"
#include "protein.hpp"
#include "tmscore.hpp"
#include "export.hpp"

using namespace std;



/**
 *  Generate VMD skript.
 *
 * \param qprot         query protein
 * \param dprot         database protein
 * \param align         original alignment
 * \param align_max     best alignment
 * \param cut_max       best cut
 * \param u             rotation of the database protein (after its transposition)
 * \param xt            transposition of the query protein
 * \param yt            transposition of the database protein
 */
void
generate_vmd(ostream &stream, cRMSDStruct *qprot, cRMSDStruct *dprot, cRMSDAlign &align, cRMSDAlign &align_max, cRMSDAlign &cut_max, cRMSDMatrix &bu, cRMSD3DCoord &bxt, cRMSD3DCoord &byt)
{
	stream << "set Q [ mol new " << QUERY_PATH << qprot->name.c_str() << " last 0 ]" << endl;

	stream << "set D [ mol new " << DATA_PATH << dprot->name.c_str() << " last 0 ]" << endl;
	stream << "set sel [atomselect $D \"all\"]" << endl;
	stream << "set M { ";

	cRMSD3DCoord mt = bxt - bu*byt;

	for (int i = 0; i < 3; i++)
	{
		stream << "{ ";

		for (int j = 0; j < 3; j++)
			stream << bu[i][j] << " ";

		stream << mt[i] << "} ";
	}

	stream << "{0 0 0 1}}" << endl;
	stream << "$sel move $M" << endl;


	stream << "set G [ mol new ]" << endl;

	for (int i = 0; i < align.length(); i++)
	{
		cRMSD3DCoord c = qprot->coord[align[i].x];
		cRMSD3DCoord d = bu*dprot->coord[align[i].y] + mt;

		stream << "graphics $G color " << ALIGN_LINE_COLOR << endl;
		stream << "graphics $G line " << c << " " << d << " width " << LINE_WIDTH << endl;

		stream << "graphics $G color " << QUERY_ALIGN_COLOR << endl;
		stream << "graphics $G sphere " << c << " radius " << BALL_RADIUS << endl;

		stream << "graphics $G color " << DATA_ALIGN_COLOR << endl;
		stream << "graphics $G sphere " << d << " radius " << BALL_RADIUS << endl;
	}



	stream << "set N [ mol new ]" << endl;

	for (int i = 0; i < align_max.length(); i++)
	{
		cRMSD3DCoord c = qprot->coord[align_max[i].x];
		cRMSD3DCoord d = bu*dprot->coord[align_max[i].y] + mt;

		stream << "graphics $N color " << MAXALIGN_LINE_COLOR << endl;
		stream << "graphics $N line " << c << " " << d << " width " << MAXLINE_WIDTH << endl;

		stream << "graphics $N color " << QUERY_ALIGN_COLOR << endl;
		stream << "graphics $N sphere " << c << " radius " << MAXBALL_RADIUS << endl;

		stream << "graphics $N color " << DATA_ALIGN_COLOR << endl;
		stream << "graphics $N sphere " << d << " radius " << MAXBALL_RADIUS << endl;
	}



	stream << "set S [ mol new ]" << endl;

	for (int i = 0; i < cut_max.length(); i++)
	{
		cRMSD3DCoord c = qprot->coord[cut_max[i].x];
		cRMSD3DCoord d = bu*dprot->coord[cut_max[i].y] + mt;

		stream << "graphics $S color " << CUTALIGN_LINE_COLOR << endl;
		stream << "graphics $S line " << c << " " << d << " width " << CUTLINE_WIDTH << endl;

		stream << "graphics $S color " << QUERY_ALIGN_COLOR << endl;
		stream << "graphics $S sphere " << c << " radius " << CUTBALL_RADIUS << endl;

		stream << "graphics $S color " << DATA_ALIGN_COLOR << endl;
		stream << "graphics $S sphere " << d << " radius " << CUTBALL_RADIUS << endl;
	}


	stream << "mol modselect 0 $Q \"resid ";

	for (int i = 0; i < align.length(); i++)
		stream << align[i].x + qprot->resid_offset << " ";

	stream << "\"" << endl;
	stream << "mol modcolor 0 $Q ColorID " << QUERY_ALIGN_COLOR << endl;
	stream << "mol modstyle 0 $Q " << QUERY_ALIGN_STYLE << endl;

	stream << "mol addrep $Q" << endl;
	stream << "mol modselect 1 $Q \"not resid ";

	for (int i = 0; i < align.length(); i++)
		stream << align[i].x + qprot->resid_offset << " ";

	stream << "\"" << endl;
	stream << "mol modcolor 1 $Q ColorID " << QUERY_NOTALIGN_COLOR << endl;
	stream << "mol modstyle 1 $Q " << QUERY_NOTALIGN_STYLE << endl;



	stream << "mol modselect 0 $D \"resid ";

	for (int i = 0; i < align.length(); i++)
		stream << align[i].y + dprot->resid_offset << " ";

	stream << "\"" << endl;
	stream << "mol modcolor 0 $D ColorID " << DATA_ALIGN_COLOR << endl;
	stream << "mol modstyle 0 $D " << DATA_ALIGN_STYLE << endl;

	stream << "mol addrep $D" << endl;
	stream << "mol modselect 1 $D \"not resid ";

	for (int i = 0; i < align.length(); i++)
		stream << align[i].y + dprot->resid_offset << " ";

	stream << "\"" << endl;
	stream << "mol modcolor 1 $D ColorID " << DATA_NOTALIGN_COLOR << endl;
	stream << "mol modstyle 1 $D " << DATA_NOTALIGN_STYLE << endl;
}
