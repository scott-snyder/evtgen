
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

#include "EvtGenBase/EvtIntervalFlatPdf.hh"

#include "EvtGenBase/EvtRandom.hh"

#include <cassert>

EvtIntervalFlatPdf::EvtIntervalFlatPdf( double min, double max ) :
    EvtPdf<EvtPoint1D>(), m_min( min ), m_max( max )
{
    assert( max >= min );
}

EvtIntervalFlatPdf::EvtIntervalFlatPdf( const EvtIntervalFlatPdf& other ) :
    EvtPdf<EvtPoint1D>( other ), m_min( other.m_min ), m_max( other.m_max )
{
}

EvtPdf<EvtPoint1D>* EvtIntervalFlatPdf::clone() const
{
    return new EvtIntervalFlatPdf( *this );
}

double EvtIntervalFlatPdf::pdf( const EvtPoint1D& ) const
{
    return 1.;
}

EvtValError EvtIntervalFlatPdf::compute_integral() const
{
    return EvtValError( m_max - m_min, 0. );
}

EvtPoint1D EvtIntervalFlatPdf::randomPoint()
{
    return EvtPoint1D( m_min, m_max, EvtRandom::Flat( m_min, m_max ) );
}
