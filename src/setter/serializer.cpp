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
#include "pdb-parser/pdb_parser.h"
#include "cRNAStructure.h"
#include "input.h"
#include "multi_align.h"

#include "parameters.h"

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <boost/serialization/string.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/base_object.hpp>

template<class Archive>
void t3DCoords::serialize(Archive & ar, const unsigned int version)
{
	ar & x;
	ar & y;
	ar & z;
}

template void t3DCoords::serialize<i_archive>(i_archive &, unsigned int);
template void t3DCoords::serialize<o_archive>(o_archive &, unsigned int);

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

template void tPdbAtom::serialize<i_archive>(i_archive &, unsigned int);
template void tPdbAtom::serialize<o_archive>(o_archive &, unsigned int);

template<class Archive>
void tPairResiude::serialize(Archive & ar, const unsigned int version)
{
	ar & chain;
	ar & residue_position;
	ar & residue_label_short;
	ar & residue_label_long;
	ar & ix_residue;
	ar & ix_chain;
}

template void tPairResiude::serialize<i_archive>(i_archive &, unsigned int);
template void tPairResiude::serialize<o_archive>(o_archive &, unsigned int);

template<class Archive>
void tPair::serialize(Archive & ar, const unsigned int version)
{
	ar & r1;
	ar & r2;
	ar & pair_type;
	ar & spatial_information;
}

template void tPair::serialize<i_archive>(i_archive &, unsigned int);
template void tPair::serialize<o_archive>(o_archive &, unsigned int);

template<class Archive>
void tHPStem::serialize(Archive & ar, const unsigned int version)
{
	ar & pairs;
}

template void tHPStem::serialize<i_archive>(i_archive &, unsigned int);
template void tHPStem::serialize<o_archive>(o_archive &, unsigned int);

template<class Archive>
void tHPHead::serialize(Archive & ar, const unsigned int version)
{
	ar & residues;
}

template void tHPHead::serialize<i_archive>(i_archive &, unsigned int);
template void tHPHead::serialize<o_archive>(o_archive &, unsigned int);

template<class Archive>
void cHairpin::serialize(Archive & ar, const unsigned int version)
{
	ar & stem;
	ar & head;
}

template void cHairpin::serialize<i_archive>(i_archive &, unsigned int);
template void cHairpin::serialize<o_archive>(o_archive &, unsigned int);

template<class Archive>
void t3DCoordsInfo::serialize(Archive & ar, const unsigned int version)
{
	ar & aa;
	ar & chain;
	ar & coords;
	ar & sse_type;
}

template void t3DCoordsInfo::serialize<i_archive>(i_archive &, unsigned int);
template void t3DCoordsInfo::serialize<o_archive>(o_archive &, unsigned int);

template<class Archive>
void tSSEHelix::serialize(Archive & ar, const unsigned int version)
{
	ar & hel_class;
	ar & id;
	ar & init_chain;
	ar & init_pos;
	ar & init_res;
	ar & term_chain;
	ar & term_pos;
	ar & term_res;
}

template void tSSEHelix::serialize<i_archive>(i_archive &, unsigned int);
template void tSSEHelix::serialize<o_archive>(o_archive &, unsigned int);

template<class Archive>
void tSSESheet::serialize(Archive & ar, const unsigned int version)
{
	ar & id;
	ar & strands_num;
	ar & init_res;
	ar & init_pos;
	ar & init_chain;
	ar & term_res;
	ar & term_pos;
	ar & term_chain;
	ar & sense;
}

template void tSSESheet::serialize<i_archive>(i_archive &, unsigned int);
template void tSSESheet::serialize<o_archive>(o_archive &, unsigned int);



template<class Archive>
void tSSETurn::serialize(Archive & ar, const unsigned int version)
{
	ar & id;
	ar & init_res;
	ar & init_pos;
	ar & init_chain;
	ar & term_res;
	ar & term_pos;
	ar & term_chain;
}

template void tSSETurn::serialize<i_archive>(i_archive &, unsigned int);
template void tSSETurn::serialize<o_archive>(o_archive &, unsigned int);

template<class Archive>
void cPDBStructure::serialize(Archive & ar, const unsigned int version)
{
	//ar & m3DCoords;
	//ar & mConect;
	//ar & mHelices;
	ar & mModels;
	//ar & mPdbID;
	//ar & mPrimarySeqChains;
	//ar & mSheets;
	//ar & mSSETypesCount;
	//ar & mTurns;
}

template void cPDBStructure::serialize<i_archive>(i_archive &, unsigned int);
template void cPDBStructure::serialize<o_archive>(o_archive &, unsigned int);

template<class Archive>
void cRNAStructure::serialize(Archive & ar, const unsigned int version)
{
	//ar & boost::serialization::base_object<cPDBStructure>(*this);
	ar & mModels;
	ar & mPDBFileName;
	ar & X3DNAPairingFileName;
	ar & mPairsChains;
	ar & mPairs;
	ar & mHairpins;
	ar & mHairpinsChains;
	ar & mModelsNt;
}

template void cRNAStructure::serialize<i_archive>(i_archive &, unsigned int);
template void cRNAStructure::serialize<o_archive>(o_archive &, unsigned int);

template<class Archive>
void Input::serialize(Archive & ar, const unsigned int version)
{
	ar & mInputSize;
	ar & mInputRnaVect;
	ar & mOriginalMap;
}

template void Input::serialize<i_archive>(i_archive &, unsigned int);
template void Input::serialize<o_archive>(o_archive &, unsigned int);

template<class Archive>
void MAResult::serialize(Archive & ar, const unsigned int version)
{
	ar & rna;
	ar & map;
	ar & dist;
	ar & mean;
	ar & deviation;
}

template void MAResult::serialize<i_archive>(i_archive &, unsigned int);
template void MAResult::serialize<o_archive>(o_archive &, unsigned int);


