
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

#include "EvtGenBase/EvtTwoBodyVertex.hh"

#include "EvtGenBase/EvtMacros.hh"

#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>

using std::endl;
using std::ostream;

// Default ctor can sometimes be useful

EvtTwoBodyVertex::EvtTwoBodyVertex() : m_LL( 0 ), m_p0( 0 )
{
}

EvtTwoBodyVertex::EvtTwoBodyVertex( double mA, double mB, double mAB, int L ) :
    m_kine(), m_LL( L ), m_p0( 0 )
{
    // Kinematics is initialized only if the decay is above threshold

    if ( mAB > mA + mB ) {
        m_kine = EvtTwoBodyKine( mA, mB, mAB );
        m_p0 = m_kine.p();
    }
}

EvtTwoBodyVertex::EvtTwoBodyVertex( const EvtTwoBodyVertex& other ) :
    m_kine( other.m_kine ),
    m_LL( other.m_LL ),
    m_p0( other.m_p0 ),
    m_f( ( other.m_f ) ? new EvtBlattWeisskopf( *other.m_f ) : nullptr )
{
}

EvtTwoBodyVertex& EvtTwoBodyVertex::operator=( const EvtTwoBodyVertex& other )
{
    m_kine = other.m_kine;
    m_LL = other.m_LL;
    m_p0 = other.m_p0;
    m_f.reset( other.m_f ? new EvtBlattWeisskopf( *other.m_f ) : nullptr );
    return *this;
}

void EvtTwoBodyVertex::set_f( double R )
{
    m_f = std::make_unique<EvtBlattWeisskopf>( m_LL, R, m_p0 );
}

double EvtTwoBodyVertex::widthFactor( EvtTwoBodyKine x ) const
{
    assert( m_p0 > 0. );

    double p1 = x.p();
    double ff = formFactor( x );
    double factor = pow( p1 / m_p0, 2 * m_LL + 1 ) * mAB() / x.mAB() * ff * ff;

    return factor;
}

double EvtTwoBodyVertex::phaseSpaceFactor( EvtTwoBodyKine x,
                                           EvtTwoBodyKine::Index i ) const
{
    double p1 = x.p( i );
    double factor = pow( p1, m_LL );
    return factor;
}

double EvtTwoBodyVertex::formFactor( EvtTwoBodyKine x ) const
{
    double ff = 1.;

    if ( m_f ) {
        double p1 = x.p();
        ff = ( *m_f )( p1 );
    }

    return ff;
}

void EvtTwoBodyVertex::print( ostream& os ) const
{
    os << " mA = " << mA() << endl;
    os << " mB = " << mB() << endl;
    os << "mAB = " << mAB() << endl;
    os << "  L = " << m_LL << endl;
    os << " p0 = " << m_p0 << endl;
}

ostream& operator<<( ostream& os, const EvtTwoBodyVertex& v )
{
    v.print( os );
    return os;
}
