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
#include <string>

#include "cRNAStructure.h"
#include "common.h"
#include "rotation.h"

using namespace std;

//! Static class allowing printing to pdb files.
class Printer
{
public:
	//! Static methode which prints RNA structure to pdb file.
	/*!
	\param rna		reference to RNA structure
	\param path		the path to th eoutput file
	\param chainId	the chian ID of the RNA structure in th epdb file
	\param append	the output pdb file will be apended, if this parameter is true
	*/
	static void PrintRNA(cRNAStructure &rna, string path, char chainId = 'A', bool append = false);
	//! Static methode which rotates and translates the RNA structure and prints it to pdb file.
	/*!
	\param rna		reference to RNA structure
	\param res		the source of the rotation matrix and translation vector
	\param path		the path to th eoutput file
	\param chainId	the chian ID of the RNA structure in th epdb file
	\param append	the output pdb file will be apended, if this parameter is true
	*/
	static void PrintRNA(cRNAStructure &rna, ResMatch &res, string path, char chainId = 'A', bool append = false);
private:
	static void PrintAtom(tPdbChain &ch, int i, FILE *file, char chainId);
	static void PrintAtom(tPdbChain &ch, cRMSD3DCoord &trans, RMatrix &matrix, int i, FILE *file, char chainId);
};