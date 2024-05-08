
/***********************************************************************
* Copyright 1998-2020 CERN for the benefit of the EvtGen authors       *
*                                                                      *
* This file is part of EvtGen.                                         *
*                                                                      *
* EvtGen is free software: you can redistribute it and/or modify       *
* it under the terms of the GNU General Public License as published by *
* the Free Software Foundation, either version 3 of the License, or    *
* (at your option) any later version.                                  *
*                                                                      *
* EvtGen is distributed in the hope that it will be useful,            *
* but WITHOUT ANY WARRANTY; without even the implied warranty of       *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        *
* GNU General Public License for more details.                         *
*                                                                      *
* You should have received a copy of the GNU General Public License    *
* along with EvtGen.  If not, see <https://www.gnu.org/licenses/>.     *
***********************************************************************/

#include "EvtGenBase/EvtDalitzCoord.hh"

#include "EvtGenBase/EvtPatches.hh"

#include <assert.h>
#include <iostream>
using EvtCyclic3::Pair;
using std::endl;
using std::ostream;

// For coordinates it's good to alway have a
// default ctor. Initialize to something invalid.

EvtDalitzCoord::EvtDalitzCoord() :
    m_i1( EvtCyclic3::AB ), m_i2( EvtCyclic3::BC ), m_q1( -1. ), m_q2( -1. )
{
}

EvtDalitzCoord::EvtDalitzCoord( const EvtDalitzCoord& other ) :
    m_i1( other.m_i1 ), m_i2( other.m_i2 ), m_q1( other.m_q1 ), m_q2( other.m_q2 )
{
}

EvtDalitzCoord::EvtDalitzCoord( Pair i1, double q1, Pair i2, double q2 ) :
    m_i1( i1 ), m_i2( i2 ), m_q1( q1 ), m_q2( q2 )
{
}

bool EvtDalitzCoord::operator==( const EvtDalitzCoord& other ) const
{
    return ( m_i1 == other.m_i1 && m_i2 == other.m_i2 && m_q1 == other.m_q1 &&
             m_q2 == other.m_q2 );
}

void EvtDalitzCoord::print( ostream& os ) const
{
    os << m_i1 << " " << m_q1 << endl;
    os << m_i2 << " " << m_q2 << endl;
}

ostream& operator<<( ostream& os, const EvtDalitzCoord& p )
{
    p.print( os );
    return os;
}
