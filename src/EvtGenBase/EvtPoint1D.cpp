
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

#include "EvtGenBase/EvtPoint1D.hh"

#include <cstdio>

EvtPoint1D::EvtPoint1D() :
    m_min( 0. ), m_max( -1. ), m_value( 0. ), m_valid( false )
{
}

EvtPoint1D::EvtPoint1D( double value ) :
    m_min( 0. ), m_max( -1. ), m_value( value ), m_valid( true )
{
}

EvtPoint1D::EvtPoint1D( double min, double max, double value ) :
    m_min( min ),
    m_max( max ),
    m_value( value ),
    m_valid( ( m_min <= m_value && m_value <= m_max ) ? true : false )
{
}

void EvtPoint1D::print() const
{
    printf( "%f (%f : %f)\n", m_value, m_min, m_max );
}
