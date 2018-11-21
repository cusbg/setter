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
#include "pdb_parser.h"
#include <errno.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <map>

#include <boost/algorithm/string.hpp>

#undef max
//#undef min

//----- cPDBStructure methods ------

cPDBStructure::cPDBStructure(string PdbFileName)
{
	this->mSSETypesCount = 3; //2 means that alpha, beta and turnes are used
	ParsePDBFile(PdbFileName);
}

cPDBStructure::~cPDBStructure()
{
}

void cPDBStructure::Clear()
{
	mPdbID = "";
	mModels.clear();
	mPrimarySeqChains.clear();
	m3DCoords.clear();
	mHelices.clear();
	mSheets.clear();
	mTurns.clear();
}

/*
 * Returns sum of letters in all primary sequence chains.
 * Beware - in PDB files this does not need be equal to getLinear3DCoords()!!!
 */
string cPDBStructure::get1LetterPrimerySeq(char chain)
{
	vector<string> seq = getPrimarySeq(chain);
	string rv;
	for (unsigned int i = 0; i < seq.size(); i++)
	{
		string aa = seq[i];
		if (aa == "ALA")
			rv += "A";
		else if (aa == "ARG")
			rv += "R";
		else if (aa == "ASN")
			rv += "N";
		else if (aa == "ASP")
			rv += "D";
		else if (aa == "CYS")
			rv += "C";
		else if (aa == "GLN")
			rv += "Q";
		else if (aa == "GLU")
			rv += "E";
		else if (aa == "GLY")
			rv += "G";
		else if (aa == "HIS")
			rv += "H";
		else if (aa == "ILE")
			rv += "I";
		else if (aa == "LEU")
			rv += "L";
		else if (aa == "LYS")
			rv += "K";
		else if (aa == "MET")
			rv += "M";
		else if (aa == "PHE")
			rv += "F";
		else if (aa == "PRO")
			rv += "P";
		else if (aa == "SER")
			rv += "S";
		else if (aa == "THR")
			rv += "T";
		else if (aa == "TRP")
			rv += "W";
		else if (aa == "TYR")
			rv += "Y";
		else if (aa == "VAL")
			rv += "V";
		else if (aa.find(' ') != string::npos)
		{
			//aa is not amino acid but nucleotide			
			rv += aa.substr(2, 1);
		}
	}

	if (rv.size() == 0)
		cout << endl << "Some kind of error in function get1LetterPrimerySeq." << endl;

	return rv;
}

vector<string> cPDBStructure::getPrimarySeq(char chain)
{
	vector<string> seq = vector<string>();
	vector<string>::iterator it;
	for (unsigned int i = 0; i < mPrimarySeqChains.size(); i++)
		if (chain == ' ' || toupper(mModels[0][i][0].chain) == toupper(chain))
		{
			seq.insert(seq.end(), mPrimarySeqChains[i].begin(), mPrimarySeqChains[i].end());
		}
	return seq;
}

/*
 * Deserializes MODEL section of a PDB file.
 */
void cPDBStructure::ReadModel(string line, ifstream &ifs, tPdbModel &m)
{
	t3DCoordsInfo i3DInfo;
	int psChains = mPrimarySeqChains.size();
	//There exist molecules (e.g., 483D) having alternative locations where some atoms are present multiple
	//times, since their exact location is not clear. Therefore we use alternateLocationIdentifier which is 
	//fixed on the first location identifier found and only atoms with this identifier are considered
	const char alternateLocationIdentifierNaN = '_';
	char alternateLocationIdentifier = alternateLocationIdentifierNaN;

	//the old version of parser relied on the PDB sequence information which is not ok since
	//some PDB file can contain only the coordinates sections
	//for (int i = 0; i < psChains; i++)
	//{
	tPdbAtom	atom;
	tPdbChain	ch;

	map<char, tPdbChain> chains;
	map<char, vector<string>> chainsPrimarySequence;

	string processedResidues = "";
	do
	{
		//i.e. 1msh is missing second chain in it's last model
		if (line.find("ENDMDL") == 0)
		{
			cout << endl << "WATCH OUT: " << mPdbID << " has wrong number of chains in one of it's models!" << endl;
			break;
		}

		//ATOM records are sometimes iterlaid by SIGATM, ANISOU and SIGUIJ records
		if (line.find("ATOM") != 0 && line.find("HETATM") != 0)
		{
			getline(ifs, line);
			continue;
		}

		atom.name = line.substr(12, 4);
		atom.alt_loc = line.c_str()[16];
		//set alternateLocationIdentifier to the first alternate location identifier met
		if (alternateLocationIdentifier == alternateLocationIdentifierNaN) alternateLocationIdentifier = atom.alt_loc;
		//if alt_loc is ' ' then there does not exist alternate location of this atom and moreover it can
		//occur after some alterenate atoms location, therefore alternateLocationIdentifier has to be reinitialized			
		if (atom.alt_loc == ' ') alternateLocationIdentifier = alternateLocationIdentifierNaN;
		//if the identifier does not match the first one met, skip this atom
		else if (alternateLocationIdentifier != atom.alt_loc)
		{
			getline(ifs, line);
			continue;
		}
		atom.res_name = boost::algorithm::trim_copy(line.substr(17, 4));
		//atom.residue_num =	atoi(line.substr(22,4).c_str());
		atom.residue_num = line.substr(22, 4).c_str();
		atom.residue_num.resize(std::remove(atom.residue_num.begin(), atom.residue_num.end(), ' ') - atom.residue_num.begin());
		atom.chain = line.c_str()[21];
		if (atom.chain == ' ')
			atom.chain = 'A';
		atom.ins_code = line.c_str()[26];
		atom.occupancy = atof(line.substr(54, 6).c_str());
		if (line.size() >= 78) atom.at_name = line.substr(76, 2);

		atom.coords.x = atof(line.substr(30, 8).c_str());
		atom.coords.y = atof(line.substr(38, 8).c_str());
		atom.coords.z = atof(line.substr(46, 8).c_str());

		line.find("ATOM") == 0 ? atom.regular = true : atom.regular = false;

		//atom.sse_type = getAtomSSE(atom.residue_num, atom.chain);

		if (chains.find(atom.chain) == chains.end()) {
			chains[atom.chain] = tPdbChain();
			chainsPrimarySequence[atom.chain] = vector<string>();
		}
		chains[atom.chain].push_back(atom);
		string resId = ";" + atom.residue_num + "-" + atom.chain + ";";
		if (processedResidues.find(resId) == string::npos)
		{
			chainsPrimarySequence[atom.chain].push_back(atom.res_name);
			processedResidues += resId;
		}
		//ch.push_back(atom);

		getline(ifs, line);
		//} while (!ifs.eof() && line.find(string("TER")) != 0);
	} while (!ifs.eof());

	for (map<char, tPdbChain>::iterator it = chains.begin(); it != chains.end(); it++) m.push_back(it->second);
	for (map<char, vector<string>>::iterator it = chainsPrimarySequence.begin(); it != chainsPrimarySequence.end(); it++) mPrimarySeqChains.push_back(it->second);
	//m.push_back(ch);
	//}
}

void cPDBStructure::ProcessPDBSSE(string line)
{
	if (line.find(string("HELIX")) == 0)
	{
		tSSEHelix helix;

		helix.id = line.substr(11, 3);
		helix.init_res = line.substr(15, 3);
		helix.init_chain = line.substr(19, 1)[0];
		if (helix.init_chain == ' ')
			helix.init_chain = 'A';
		helix.init_pos = atoi(line.substr(21, 4).c_str());
		helix.term_res = line.substr(27, 3);
		helix.term_chain = line.substr(31, 1)[0];
		if (helix.term_chain == ' ')
			helix.term_chain = 'A';
		helix.term_pos = atoi(line.substr(33, 4).c_str());
		helix.hel_class = atoi(line.substr(38, 2).c_str());

		mHelices.push_back(helix);
	}
	if (line.find(string("SHEET")) == 0)
	{
		tSSESheet sheet;

		sheet.id = line.substr(11, 3);
		sheet.strands_num = atoi(line.substr(14, 2).c_str());
		sheet.init_res = line.substr(17, 3);
		sheet.init_chain = line.substr(21, 1)[0];
		if (sheet.init_chain == ' ')
			sheet.init_chain = 'A';
		sheet.init_pos = atoi(line.substr(22, 4).c_str());
		sheet.term_res = line.substr(28, 3);
		sheet.term_chain = line.substr(32, 1)[0];
		if (sheet.term_chain == ' ')
			sheet.term_chain = 'A';
		sheet.term_pos = atoi(line.substr(33, 4).c_str());
		sheet.sense = atoi(line.substr(38, 2).c_str());

		mSheets.push_back(sheet);
	}
	if (line.find(string("TURN")) == 0)
	{
		tSSETurn turn;

		turn.id = line.substr(11, 3);
		turn.init_res = line.substr(15, 3);
		turn.init_chain = line.substr(19, 1)[0];
		if (turn.init_chain == ' ')
			turn.init_chain = 'A';
		turn.init_pos = atoi(line.substr(20, 4).c_str());
		turn.term_res = line.substr(26, 3);
		turn.term_chain = line.substr(30, 1)[0];
		if (turn.term_chain == ' ')
			turn.term_chain = 'A';
		turn.term_pos = atoi(line.substr(31, 4).c_str());

		mTurns.push_back(turn);
	}
}

void cPDBStructure::ProcessPDBCONECT(string line)
{
	pair< int, vector<int> > conect;
	string temp;

	int k;
	conect.first = atoi(line.substr(6, 10).c_str());
	temp = line.substr(11, 5); k = temp.find_first_not_of(' '); if (k != -1) conect.second.push_back(atoi(temp.substr(k).c_str()));
	temp = line.substr(16, 5); k = temp.find_first_not_of(' '); if (k != -1) conect.second.push_back(atoi(temp.substr(k).c_str()));
	temp = line.substr(21, 5); k = temp.find_first_not_of(' '); if (k != -1) conect.second.push_back(atoi(temp.substr(k).c_str()));
	temp = line.substr(26, 5); k = temp.find_first_not_of(' '); if (k != -1) conect.second.push_back(atoi(temp.substr(k).c_str()));

	unsigned int i = 0;
	for (i = 0; i < mConect.size(); i++)
	{
		if (mConect[i].first == conect.first)
		{
			break;
		}
	}

	if (i < mConect.size())
	{
		for (unsigned int j = 0; j < conect.second.size(); j++)
		{
			mConect[i].second.push_back(conect.second[j]);
		}
	}
	else
	{
		mConect.push_back(conect);
	}
}

void cPDBStructure::ParsePDB(ifstream &ifs)
{
	string line;
	bool AtomRecProcessed = false;	//i.e. pdb1abw has more TER records than chains in SEQRES,
	//hence after processing the designated chains no more reading
	//of ATOM records occures (it looks like some short chains)
	while (!ifs.eof())
	{
		//read line from pdb file
		getline(ifs, line);
		//read PDB ID from HEADER
		if (line.find("HEADER") == 0 && line.size()>=66)
			mPdbID = line.substr(62, 4);

		//read primary sequence in form of multiple chains from SEQRES
		//MOD: The primary sequence is now determined from the coordinates section
		/*
		if (line.find(string("SEQRES")) == 0)
		{
		bool first_line = true;
		string chain;
		do
		{
		if ( atoi(line.substr(8,2).c_str()) == 1 && !first_line) //line.c_str()[11] != 'A' ) not all chains start from A (as in pdb12e8)
		{
		mPrimarySeqChains.push_back(chain);
		chain = "";
		}
		first_line = false;
		for (int i = 0; i < 13; i++)
		{
		if (line.substr(19+i*4, 4) == "    ")
		break;
		chain.append(line.substr(19+i*4, 4));
		}
		getline(ifs, line);
		} while (line.find("SEQRES") == 0);
		if (chain != "")
		mPrimarySeqChains.push_back(chain);
		}
		*/

		//read tertiary structre (1 MODEL)
		//if ((line.find(string("ATOM")) == 0 && line.find(string("HETATM")) == 0 && !AtomRecProcessed) || line.find(string("MODEL")) == 0)
		//if (line.find(string("ATOM")) == 0 || line.find(string("HETATM")) == 0 || line.find(string("MODEL")) == 0)
		if ((line.find(string("ATOM")) == 0 && !AtomRecProcessed) || line.find(string("MODEL")) == 0)
		{
			tPdbModel m;

			if (line.find(string("MODEL")) == 0)
				getline(ifs, line);
			ReadModel(line, ifs, m);
			mModels.push_back(m);
			if (line.find(string("ATOM")) == 0 || line.find(string("HETATM")) == 0)
				AtomRecProcessed = true;
		}

		if (line.find(string("HELIX")) == 0 || line.find(string("SHEET")) == 0 || line.find(string("TURN")) == 0)
			ProcessPDBSSE(line);

		if (line.find(string("CONECT")) == 0)
		{
			ProcessPDBCONECT(line);
		}

	}
}


void cPDBStructure::AddSSEInfo(ifstream &ifs)
{
	string line;
	while (!ifs.eof())
	{
		//read line from pdb file
		getline(ifs, line);

		if (line.find(string("HELIX")) == 0 || line.find(string("SHEET")) == 0 || line.find(string("TURN")) == 0)
			ProcessPDBSSE(line);
	}
}


void cPDBStructure::ParsePDBFile(string ptName)
{
	Clear();
	string fileName;
	string dir;

	fileName = ptName;
	dir = PDB_LIB_DIR;
	ifstream ifs((dir + fileName).c_str());
	if (!ifs)
	{
		throw dir + fileName + " could not be opened.";
	}

	ParsePDB(ifs);

	if (mPdbID == "") mPdbID = ptName;
}

string SerializeAtom(tPdbAtom a)
{
	ostringstream str;
	str << a.coords.x << ";" << a.coords.y << ";" << a.coords.z << ";";
	str << a.name << ";" << a.res_name << ";" << a.chain << ";"
		<< a.at_name << ";" << a.alt_loc << ";" << a.residue_num
		<< ";" << a.ins_code << ";" << a.occupancy << ";";

	return str.str();
}

tPdbAtom DesrializeAtom(istringstream &istr)
{
	tPdbAtom a;
	string token;

	getline(istr, token, ';');
	a.coords.x = atof(token.c_str());
	getline(istr, token, ';');
	a.coords.y = atof(token.c_str());
	getline(istr, token, ';');
	a.coords.z = atof(token.c_str());
	getline(istr, token, ';');
	a.name = token;
	getline(istr, token, ';');
	a.res_name = token;
	getline(istr, token, ';');
	a.chain = token[0];
	getline(istr, token, ';');
	a.at_name = token;
	getline(istr, token, ';');
	a.alt_loc = token[0];
	getline(istr, token, ';');
	a.residue_num = token.c_str();
	getline(istr, token, ';');
	a.ins_code = token[0];
	getline(istr, token, ';');
	a.occupancy = atof(token.c_str());

	return a;
}

void cPDBStructure::Serialize(ofstream &ofs)
{
	ofs << mPdbID << endl;
	ofs << mPrimarySeqChains.size() << endl;
	for (unsigned int i = 0; i < mPrimarySeqChains.size(); i++)
	{
		for (unsigned int j = 0; j < mPrimarySeqChains[i].size(); ofs << mPrimarySeqChains[i][j++] + " ");
		ofs << endl;
	}
	ofs << mModels.size() << endl;
	//for each mode
	for (unsigned int i = 0; i < mModels.size(); i++)
		//for each chain in each model
		for (unsigned int j = 0; j < mPrimarySeqChains.size(); j++)
		{
			//for each atom
			for (unsigned int k = 0; k < mModels[i][j].size(); k++)
				ofs << SerializeAtom(mModels[i][j][k]);
			ofs << endl;
		}
}

void cPDBStructure::Serialize(string fileName)
{
	if (fileName == "")
		fileName = EXP_PDB_LIB_DIR + mPdbID + ".spdb";
	else
		fileName = EXP_PDB_LIB_DIR + fileName;

	ofstream ofs(fileName.c_str());
	if (!ofs)
	{
		cerr << "Unable to open file " << fileName << endl;
		exit(1);
	}

	Serialize(ofs);
}

void cPDBStructure::Deserialize(ifstream &ifs)
{
	string line;
	int size;
	getline(ifs, mPdbID);
	getline(ifs, line);
	size = atoi(line.c_str());
	for (int i = 0; i < size; i++)
	{
		getline(ifs, line);
		boost::algorithm::trim(line);
		vector<string> strs;
		boost::split(strs, line, boost::is_any_of(" "));
		mPrimarySeqChains.push_back(strs);
	}
	getline(ifs, line);
	size = atoi(line.c_str());
	for (int i = 0; i < size; i++)
	{
		tPdbModel	m;
		for (unsigned int j = 0; j < mPrimarySeqChains.size(); j++)
		{
			getline(ifs, line);
			istringstream istr(line);
			tPdbChain	ch;
			t3DCoordsInfo coordsInfo;
			while (!istr.eof())
			{
				tPdbAtom a = DesrializeAtom(istr);
				ch.push_back(a);
				char c = (char)istr.get();
				if (c != -1)
					istr.putback(c);
			}
			m.push_back(ch);
		}
		mModels.push_back(m);
	}
}

void cPDBStructure::Deserialize(string fileName)
{
	ifstream ifs((EXP_PDB_LIB_DIR + fileName).c_str());
	if (!ifs)
	{
		cerr << "Unable to open file " << fileName << endl;
		exit(1);
	}
	Deserialize(ifs);
}

void cPDBStructure::RenameChains(char startWithCharacter)
{
	for (int ixChain = 0; ixChain < mModels[0].size(); ixChain++)
	{
		for (int ixRes = 0; ixRes < mModels[0].size(); ixRes++) mModels[0][ixChain][ixRes].chain = startWithCharacter + ixChain;
	}
}
