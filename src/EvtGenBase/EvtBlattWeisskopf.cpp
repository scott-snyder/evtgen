
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

#include "EvtGenBase/EvtBlattWeisskopf.hh"

#include "EvtGenBase/EvtReport.hh"

#include <cassert>
#include <cmath>
#include <iostream>

using std::endl;

EvtBlattWeisskopf::EvtBlattWeisskopf( int LL, double R, double p0 ) :
    m_LL( LL ), m_radial( R ), m_p0( p0 )
{
    if ( R < 0 ) {
        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << "Radius " << R << " negative" << endl;
        assert( 0 );
    }

    m_radial = R;

    // compute formula for nominal momentum

    m_F0 = compute( m_p0 );
    if ( m_F0 <= 0 ) {
        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << "Invalid nominal form factor computed " << m_F0 << endl;
        assert( 0 );
    }
}

EvtBlattWeisskopf::EvtBlattWeisskopf( const EvtBlattWeisskopf& other ) :
    m_LL( other.m_LL ),
    m_radial( other.m_radial ),
    m_p0( other.m_p0 ),
    m_F0( other.m_F0 )
{
}

double EvtBlattWeisskopf::operator()( double p ) const
{
    double ret = compute( p ) / m_F0;
    //  EvtGenReport(EVTGEN_INFO,"EvtGen") << p << " " << m_p0 << " " << m_F0 << " " << m_LL << " " << m_radial << " " << ret << endl;
    return ret;
}

// Blatt-Weisskopf form factors
// see e.g. hep-ex/0011065
// Dalitz Analysis of the Decay D0->K-pi+pi0 (CLEO)
//
// p   - momentum of either daugher in the meson rest frame,
//       the mass of the meson is used
// pAB - momentum of either daughter in the candidate rest frame
//       the mass of the candidate is used
// R - meson radial parameter
//
// In the CLEO paper R=5 GeV-1 for D0, R=1.5 for intermediate resonances

double EvtBlattWeisskopf::compute( double p ) const
{
    double value( 1.0 );

    double z = p * m_radial;
    double zSq = z * z;

    if ( m_LL == 0 ) {
        value = 1.0;
    } else if ( m_LL == 1 ) {
        value = sqrt( 1.0 / ( 1.0 + zSq ) );
    } else if ( m_LL == 2 ) {
        value = sqrt( 1.0 / ( zSq * ( zSq + 3.0 ) + 9.0 ) );
    } else if ( m_LL == 3 ) {
        double denom = zSq * ( zSq * ( zSq + 6.0 ) + 45.0 ) + 225.0;
        value = sqrt( 1.0 / denom );
    } else if ( m_LL == 4 ) {
        double denom = zSq * ( zSq * ( zSq * ( zSq + 10.0 ) + 135.0 ) + 1575.0 ) +
                       11025.0;
        value = sqrt( 1.0 / denom );
    } else if ( m_LL == 5 ) {
        double denom =
            zSq * ( zSq * ( zSq * ( zSq * ( zSq + 15.0 ) + 315.0 ) + 6300.0 ) +
                    99225.0 ) +
            893025.0;
        value = sqrt( 1.0 / denom );
    }

    return value;
}
