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
#ifndef PDB_PARSER_H_
#define PDB_PARSER_H_

/*
 * Notes to ATOM PDB record description - http://bmerc-www.bu.edu/needle-doc/latest/atom-format.html
*/
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <algorithm>
#include <limits>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/base_object.hpp>

//#define PDB_LIB_DIR "e:\\biological_data\\PDB\\PDBFormat\\unzipped\\"
//#define PDB_LIB_DIR "..\\data\\pdb\\"
#define PDB_LIB_DIR ""
#define EXP_PDB_LIB_DIR "e:\\biological_data\\PDB\\PDBFormat\\experiments\\"
#define PSIST_COORD_DIR "PSIST\\PDB_coord\\"

#define ASTRAL_1_65_LIB_DIR "e:\\biological_data\\PDB\\ASTRAL\\1.65\\"
#define ASTRAL_1_67_LIB_DIR "e:\\biological_data\\PDB\\ASTRAL\\1.67\\"
#define ASTRAL_1_73_LIB_DIR "e:\\biological_data\\PDB\\ASTRAL\\1.73\\"

#define MAX_SEQ_LEN 2000

using namespace std;

typedef struct t3DCoords{
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version);

	double	x, y, z;	

	t3DCoords():x(0),y(0),z(0){};
	t3DCoords(double _x, double _y, double _z): x(_x), y(_x), z(_z){};

	t3DCoords	operator+=(const t3DCoords& b) { x += b.x; y += b.y; z += b.z; return *this; }
	t3DCoords	operator+=(const double& add) { x += add; y += add; z += add; return *this; }
	t3DCoords	operator/=(const double& divisor) { x /= divisor; y /= divisor; z /= divisor; return *this; }
	double &    operator[](int ix) { if (ix == 0) return x; else if (ix == 1) return y; else return z; }
	
	inline double DistL2(t3DCoords coords){return pow(pow(x-coords.x, 2)+pow(y-coords.y, 2)+pow(z-coords.z, 2),0.5); };

};

typedef struct {
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version);

	t3DCoords	coords;
	char		sse_type; //see tPdbAtom definition for explanation
	char		chain;
	string		aa;
} t3DCoordsInfo; //3D coordinations with usefull info for VPT generation a.s.o.

struct sDist
	{
	public:

		string str;
		double number;
		bool operator<(const sDist& b)const { return number < b.number; }
		sDist(){str=""; number=0;};
		sDist(const sDist& b){ number = b.number; str = b.str; }
		~sDist() {};
		sDist operator=(const sDist& b) {number = b.number; str = b.str; return b; }
	};	


/*struct sortAscending3D
{
	bool operator()(t3DCoords*& a, t3DCoords*& b)
     {
          return rpStart->GetTimestamp() < rpEnd->GetTimestamp();
     }
};*/

enum eSource{
	ePDB = 0,
	ePSIST, 
	eAstral165,
	eAstral167
};

typedef struct {
	template<class Archive>
    void serialize(Archive & ar, const unsigned int version);

	t3DCoords		coords;
	string			name;			//name with alpha, beta, gama indicator (can by nucleid acid) - 4letters
	string			res_name;		//residue name - code of the AA (VAL, ALA, ...) of NA (DA, DG, DC, DT) - 4 letters
	char			chain;
	string			at_name;		//atom name (C, O, H, ...)	- 2 letters
	char			alt_loc;		//alternate location
	string			residue_num;	//number of AA in which the atome resides
	char			ins_code;		//insertion code
	double			occupancy;
	char			sse_type;			//A-alpha (possibly A1, A2, A3, ...), B, T (turn)
	bool			regular;		//whether it is ATOM or HETATM

} tPdbAtom;

typedef struct {
	template<class Archive>
    void serialize(Archive & ar, const unsigned int version);

	string	id;
	string	init_res;
	int		init_pos;
	char	init_chain;
	string	term_res;
	int		term_pos;
	char	term_chain;
	int		hel_class;

} tSSEHelix;

typedef struct {
	template<class Archive>
    void serialize(Archive & ar, const unsigned int version);

	string	id;
	int		strands_num;
	string	init_res;
	int		init_pos;
	char	init_chain;
	string	term_res;
	int		term_pos;
	char	term_chain;
	int		sense; //Sense of strand with respect to previous strand in the sheet. 0 if first strand, 1 if parallel, -1 if anti-parallel.

} tSSESheet;

typedef struct {
	template<class Archive>
    void serialize(Archive & ar, const unsigned int version);

	string	id;
	string	init_res;
	int		init_pos;
	char	init_chain;
	string	term_res;
	int		term_pos;
	char	term_chain;
} tSSETurn;

typedef vector<tPdbAtom>		tPdbChain;
typedef vector<tPdbChain>		tPdbModel;

class cPDBStructure
{
	/*
	 * The class is meant to be "pure" PDB reader. All operations not directly connected with the PDB
	 * format (such as C-alpha atoms extraction, distance computations, modification to be able to read
	 * ASTRAL files, etc.) should be done in derived classes.
	*/
private:
	friend class boost::serialization::access;
	template<class Archive>
    void serialize(Archive & ar, const unsigned int version);

public:
	string					mPdbID;
	vector<vector<string>>	mPrimarySeqChains;
	vector<tPdbModel>		mModels;
	vector<t3DCoordsInfo>	m3DCoords;		

	vector<tSSEHelix>		mHelices;
	vector<tSSESheet>		mSheets;
	vector<tSSETurn>		mTurns;

	vector< pair< int, vector<int> > >	mConect; //connection between units (amino acids, nucleic acids) - each unit in connection contains list
									//of units to which it is connected

	int						mSSETypesCount; //generally 3, if we do not differntiate types of alpha helices
	
public:

	cPDBStructure() {};
	cPDBStructure(string PdbFileName);	

	~cPDBStructure();

	void				Clear();

	void				ParsePDBFile(string ptName);
	void				ParsePDB(ifstream &ifs);
	void				ParsePSIST(ifstream &ifs);
	void				ParseAstral(ifstream &ifs);
	//char				getAtomSSE(string atom_pos, char chain);
	void				ProcessPDBSSE(string line);
	void				ProcessPDBCONECT(string line);

	tPdbAtom			getPdbAtom(unsigned int i);
	vector<string>		getPrimarySeq(char chain = ' ');
	string				get1LetterPrimerySeq(char chain = ' ');
	vector<t3DCoords>	getLinear3DCoords(char chain = ' ');
	vector<t3DCoordsInfo> getEnhancedLinear3DCoords(char chain = ' ');

	void				Serialize(string fileName = "");
	void				Serialize(ofstream &ofs);
	void				Deserialize(string fileName);
	void				Deserialize(ifstream &ifs);

	string				getId(){ return mPdbID; };
	unsigned int		getPrimarySeqLen() { return get1LetterPrimerySeq().length(); };
	unsigned int		getChainsCount() { return mPrimarySeqChains.size(); };
	unsigned int		getModelCount() { return mModels.size(); };

	void				RenameChains(char startWithCharacter);

protected:	
	void	ReadModel(string line, ifstream &ifs, tPdbModel &m);
	void	AddSSEInfo(ifstream &ifs);
	
};

#endif /*PDB_PARSER_H_*/
