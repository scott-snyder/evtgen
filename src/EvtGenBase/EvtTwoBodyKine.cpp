
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

#include "EvtGenBase/EvtTwoBodyKine.hh"

#include "EvtGenBase/EvtReport.hh"

#include <assert.h>
#include <iostream>
#include <math.h>
using std::endl;
using std::ostream;

EvtTwoBodyKine::EvtTwoBodyKine() : m_mA( 0. ), m_mB( 0. ), m_mAB( 0. )
{
}

EvtTwoBodyKine::EvtTwoBodyKine( const double mA, const double mB,
                                const double mAB ) :
    m_mA( mA ), m_mB( mB ), m_mAB( mAB )
{
    if ( mAB < mA + mB ) {
        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << mAB << " < " << mA << " + " << mB << endl;
        assert( 0 );
    }
}

double EvtTwoBodyKine::m( Index i ) const
{
    double ret = m_mAB;
    if ( A == i )
        ret = m_mA;
    else if ( B == i )
        ret = m_mB;

    return ret;
}

double EvtTwoBodyKine::p( Index i ) const
{
    double p0 = 0.;

    if ( i == AB ) {
        double x = m_mAB * m_mAB - m_mA * m_mA - m_mB * m_mB;
        double y = 2 * m_mA * m_mB;
        p0 = sqrt( x * x - y * y ) / 2. / m_mAB;
    } else if ( i == A ) {
        double x = m_mA * m_mA - m_mAB * m_mAB - m_mB * m_mB;
        double y = 2 * m_mAB * m_mB;
        p0 = sqrt( x * x - y * y ) / 2. / m_mA;
    } else {
        double x = m_mB * m_mB - m_mAB * m_mAB - m_mA * m_mA;
        double y = 2 * m_mAB * m_mA;
        p0 = sqrt( x * x - y * y ) / 2. / m_mB;
    }

    return p0;
}

double EvtTwoBodyKine::e( Index i, Index j ) const
{
    double ret = m( i );
    if ( i != j ) {
        double pD = p( j );
        ret = sqrt( ret * ret + pD * pD );
    }
    return ret;
}

void EvtTwoBodyKine::print( ostream& os ) const
{
    os << " mA = " << m_mA << endl;
    os << " mB = " << m_mB << endl;
    os << "mAB = " << m_mAB << endl;
}

ostream& operator<<( ostream& os, const EvtTwoBodyKine& p )
{
    p.print( os );
    return os;
}
