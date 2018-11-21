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
#include "print_to_pdb.h"

#include <boost/algorithm/string.hpp>

void Printer::PrintRNA(cRNAStructure &rna, string path, char chainId, bool append)
{
	FILE *file;
	if (append) {
		file = fopen(path.c_str(), "a");
	}
	else {
		file = fopen(path.c_str(), "w");
	}

	tPdbChain &ch = rna.mModels[0][0];
	tPdbAtom atom;
	//cout << "0         1         2         3         4         5         6         7         ";
	//cout << "01234567890123456789012345678901234567890123456789012345678901234567890123456789";

	for (int i = 0; i < ch.size(); ++i) {
		PrintAtom(ch, i, file, chainId);
	}

	fprintf(file, "TER   ");
	fprintf(file, "%5d", ch.size() + 1);
	fprintf(file, "          ");
	fprintf(file, "%c", chainId);
	fprintf(file, "\n");

	fclose(file);
}

void Printer::PrintRNA(cRNAStructure &rna, ResMatch &res, string path, char chainId, bool append)
{
	FILE *file;
	if (append) {
		file = fopen(path.c_str(), "a");
	}
	else {
		file = fopen(path.c_str(), "w");
	}

	tPdbChain &ch = rna.mModels[0][0];
	tPdbAtom atom;
	//cout << "0         1         2         3         4         5         6         7         ";
	//cout << "01234567890123456789012345678901234567890123456789012345678901234567890123456789";

	RMatrix matrix(res.rot);

	for (int i = 0; i < ch.size(); ++i) {
		PrintAtom(ch, res.trans, matrix, i, file, chainId);
	}

	fprintf(file, "TER   ");
	fprintf(file, "%5d", ch.size() + 1);
	fprintf(file, "          ");
	fprintf(file, "%c", chainId);
	fprintf(file, "\n");

	fclose(file);
}

void Printer::PrintAtom(tPdbChain &ch, int i, FILE *file, char chainId)
{
	tPdbAtom &atom = ch[i];


	fprintf(file, "ATOM  ");
	fprintf(file, "%5d", i + 1);
	fprintf(file, " ");
	fprintf(file, "%4s", boost::trim_copy(atom.name).c_str());
	fprintf(file, " "); //Alternate location indicator
	fprintf(file, "%3s", boost::trim_copy(atom.res_name).c_str());
	fprintf(file, " "); //Position 21 not used
	fprintf(file, "%c", chainId);			// chId
	fprintf(file, "%4s", boost::trim_copy(atom.residue_num).c_str());
	fprintf(file, " ");			// insertion code
	fprintf(file, "   ");
	fprintf(file, "%8.3f", atom.coords.x);
	fprintf(file, "%8.3f", atom.coords.y);
	fprintf(file, "%8.3f", atom.coords.z);
	fprintf(file, "  1.00");	// occupancy
	fprintf(file, "  0.00");	// tempFactor
	fprintf(file, "      ");
	fprintf(file, "    ");	// recID
	fprintf(file, "%s", atom.at_name.c_str());
	fprintf(file, "  ");
	fprintf(file, "\n");
}


void Printer::PrintAtom(tPdbChain &ch, cRMSD3DCoord &trans, RMatrix &matrix, int i, FILE *file, char chainId)
{
	tPdbAtom &atom = ch[i];


	cRMSD3DCoord coord;
	coord[0] = atom.coords.x;
	coord[1] = atom.coords.y;
	coord[2] = atom.coords.z;

	coord = matrix * coord;
	coord = coord + trans;
	//coord = matrix * coord;


	fprintf(file, "ATOM  ");
	fprintf(file, "%5d", i + 1);
	fprintf(file, " ");
	fprintf(file, "%4s", boost::trim_copy(atom.name).c_str());
	fprintf(file, " "); //Alternate location indicator
	fprintf(file, "%3s", boost::trim_copy(atom.res_name).c_str());
	fprintf(file, " ");
	fprintf(file, "%c", chainId);			// chId
	fprintf(file, "%4s", boost::trim_copy(atom.residue_num).c_str());
	fprintf(file, " ");			// insertion code
	fprintf(file, "   ");
	fprintf(file, "%8.3f", coord[0]);
	fprintf(file, "%8.3f", coord[1]);
	fprintf(file, "%8.3f", coord[2]);
	fprintf(file, "  1.00");	// occupancy
	fprintf(file, "  0.00");	// tempFactor
	fprintf(file, "      ");
	fprintf(file, "    ");	// recID
	fprintf(file, "%s", atom.at_name.c_str());
	fprintf(file, "  ");
	fprintf(file, "\n");
}
