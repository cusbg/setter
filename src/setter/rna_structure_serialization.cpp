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
#include "cRNAStructure.h"

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <boost/serialization/string.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>

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

template void tPairResiude::serialize<class boost::archive::text_iarchive>(class boost::archive::text_iarchive &, unsigned int);
template void tPairResiude::serialize<class boost::archive::text_oarchive>(class boost::archive::text_oarchive &, unsigned int);

template<class Archive>
void tPair::serialize(Archive & ar, const unsigned int version)
{
	ar & r1;
	ar & r2;
	ar & pair_type;
	ar & spatial_information;
}

template void tPair::serialize<class boost::archive::text_iarchive>(class boost::archive::text_iarchive &, unsigned int);
template void tPair::serialize<class boost::archive::text_oarchive>(class boost::archive::text_oarchive &, unsigned int);

template<class Archive>
void tHPStem::serialize(Archive & ar, const unsigned int version)
{
	ar & pairs;
}

template void tHPStem::serialize<class boost::archive::text_iarchive>(class boost::archive::text_iarchive &, unsigned int);
template void tHPStem::serialize<class boost::archive::text_oarchive>(class boost::archive::text_oarchive &, unsigned int);

template<class Archive>
void tHPHead::serialize(Archive & ar, const unsigned int version)
{
	ar & residues;
}

template void tHPHead::serialize<class boost::archive::text_iarchive>(class boost::archive::text_iarchive &, unsigned int);
template void tHPHead::serialize<class boost::archive::text_oarchive>(class boost::archive::text_oarchive &, unsigned int);

template<class Archive>
void cHairpin::serialize(Archive & ar, const unsigned int version)
{
	ar & stem;
	ar & head;
}

template void cHairpin::serialize<class boost::archive::text_iarchive>(class boost::archive::text_iarchive &, unsigned int);
template void cHairpin::serialize<class boost::archive::text_oarchive>(class boost::archive::text_oarchive &, unsigned int);

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

template void cRNAStructure::serialize<class boost::archive::text_iarchive>(class boost::archive::text_iarchive &, unsigned int);
template void cRNAStructure::serialize<class boost::archive::text_oarchive>(class boost::archive::text_oarchive &, unsigned int);