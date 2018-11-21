#include "../../pdb_parser/pdb_parser.h"


template<class Archive>
void tPdbAtom::serialize(Archive & ar, const unsigned int version)
{
	ar & coords;
	ar & name;
	ar & res_name;
	ar & chain;
	ar & at_name;
	ar & alt_loc;
	ar & residue_num;
	ar & ins_code;
	ar & occupancy;
	ar & sse_type;
	ar & regular;
}