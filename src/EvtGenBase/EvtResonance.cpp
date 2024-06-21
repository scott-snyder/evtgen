
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

#include "EvtGenBase/EvtResonance.hh"

#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtConst.hh"
#include "EvtGenBase/EvtKine.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtVector4R.hh"

#include <math.h>
using std::endl;

EvtResonance& EvtResonance::operator=( const EvtResonance& n )
{
    if ( &n == this )
        return *this;
    m_p4_p = n.m_p4_p;
    m_p4_d1 = n.m_p4_d1;
    m_p4_d2 = n.m_p4_d2;
    m_ampl = n.m_ampl;
    m_theta = n.m_theta;
    m_gamma = n.m_gamma;
    m_spin = n.m_spin;
    m_bwm = n.m_bwm;
    return *this;
}

EvtResonance::EvtResonance( const EvtVector4R& p4_p, const EvtVector4R& p4_d1,
                            const EvtVector4R& p4_d2, double ampl, double theta,
                            double gamma, double bwm, int spin ) :
    m_p4_p( p4_p ),
    m_p4_d1( p4_d1 ),
    m_p4_d2( p4_d2 ),
    m_ampl( ampl ),
    m_theta( theta ),
    m_gamma( gamma ),
    m_bwm( bwm ),
    m_spin( spin )
{
}

EvtComplex EvtResonance::resAmpl()
{
    double pi180inv = 1.0 / EvtConst::radToDegrees;

    EvtComplex ampl;
    //EvtVector4R  _p4_d3 = m_p4_p-m_p4_d1-m_p4_d2;

    //get cos of the angle between the daughters from their 4-momenta
    //and the 4-momentum of the parent

    //in general, EvtDecayAngle(parent, part1+part2, part1) gives the angle
    //the missing particle (not listed in the arguments) makes
    //with part2 in the rest frame of both
    //listed particles (12)

    //angle 3 makes with 2 in rest frame of 12 (CS3)
    double cos_phi_0 = EvtDecayAngle( m_p4_p, m_p4_d1 + m_p4_d2, m_p4_d1 );
    //angle 3 makes with 1 in 12 is, of course, -cos_phi_0

    switch ( m_spin ) {
        case 0:
            ampl = ( m_ampl *
                     EvtComplex( cos( m_theta * pi180inv ),
                                 sin( m_theta * pi180inv ) ) *
                     sqrt( m_gamma / EvtConst::twoPi ) *
                     ( 1.0 / ( ( m_p4_d1 + m_p4_d2 ).mass() - m_bwm -
                               EvtComplex( 0.0, 0.5 * m_gamma ) ) ) );
            break;

        case 1:
            ampl = ( m_ampl *
                     EvtComplex( cos( m_theta * pi180inv ),
                                 sin( m_theta * pi180inv ) ) *
                     sqrt( m_gamma / EvtConst::twoPi ) *
                     ( cos_phi_0 / ( ( m_p4_d1 + m_p4_d2 ).mass() - m_bwm -
                                     EvtComplex( 0.0, 0.5 * m_gamma ) ) ) );
            break;

        case 2:
            ampl = ( m_ampl *
                     EvtComplex( cos( m_theta * pi180inv ),
                                 sin( m_theta * pi180inv ) ) *
                     sqrt( m_gamma / EvtConst::twoPi ) *
                     ( ( 1.5 * cos_phi_0 * cos_phi_0 - 0.5 ) /
                       ( ( m_p4_d1 + m_p4_d2 ).mass() - m_bwm -
                         EvtComplex( 0.0, 0.5 * m_gamma ) ) ) );
            break;

        case 3:
            ampl = ( m_ampl *
                     EvtComplex( cos( m_theta * pi180inv ),
                                 sin( m_theta * pi180inv ) ) *
                     sqrt( m_gamma / EvtConst::twoPi ) *
                     ( ( 2.5 * cos_phi_0 * cos_phi_0 * cos_phi_0 -
                         1.5 * cos_phi_0 ) /
                       ( ( m_p4_d1 + m_p4_d2 ).mass() - m_bwm -
                         EvtComplex( 0.0, 0.5 * m_gamma ) ) ) );
            break;

        default:
            EvtGenReport( EVTGEN_DEBUG, "EvtGen" )
                << "EvtGen: wrong spin in EvtResonance" << endl;
            ampl = EvtComplex( 0.0 );
            break;
    }

    return ampl;
}

EvtComplex EvtResonance::relBrWig( int i )
{
    //this function returns relativistic Breit-Wigner amplitude
    //for a given resonance (for P-wave decays of scalars only at the moment!)

    EvtComplex BW;
    EvtVector4R _p4_d3 = m_p4_p - m_p4_d1 - m_p4_d2;
    EvtVector4R _p4_12 = m_p4_d1 + m_p4_d2;

    double msq13 = ( m_p4_d1 + _p4_d3 ).mass2();
    double msq23 = ( m_p4_d2 + _p4_d3 ).mass2();
    double msqParent = m_p4_p.mass2();
    double msq1 = m_p4_d1.mass2();
    double msq2 = m_p4_d2.mass2();
    double msq3 = _p4_d3.mass2();

    double M;

    double p2 =
        sqrt( ( _p4_12.mass2() - ( m_p4_d1.mass() + m_p4_d2.mass() ) *
                                     ( m_p4_d1.mass() + m_p4_d2.mass() ) ) *
              ( _p4_12.mass2() - ( m_p4_d1.mass() - m_p4_d2.mass() ) *
                                     ( m_p4_d1.mass() - m_p4_d2.mass() ) ) ) /
        ( 2.0 * _p4_12.mass() );

    double p2R =
        sqrt( ( m_bwm * m_bwm - ( m_p4_d1.mass() + m_p4_d2.mass() ) *
                                    ( m_p4_d1.mass() + m_p4_d2.mass() ) ) *
              ( m_bwm * m_bwm - ( m_p4_d1.mass() - m_p4_d2.mass() ) *
                                    ( m_p4_d1.mass() - m_p4_d2.mass() ) ) ) /
        ( 2.0 * m_bwm );

    double gam, R;

    if ( i == 1 ) {
        R = 2.0 / ( 0.197 );

    } else
        R = 5.0 / ( 0.197 );

    gam = m_gamma * ( m_bwm / _p4_12.mass() ) * ( p2 / p2R ) * ( p2 / p2R ) *
          ( p2 / p2R ) * ( ( 1 + R * R * p2R * p2R ) / ( 1 + R * R * p2 * p2 ) );
    M = ( msq13 - msq23 -
          ( msqParent - msq3 ) * ( msq1 - msq2 ) / ( m_bwm * m_bwm ) ) *
        sqrt( ( 1 + R * R * p2R * p2R ) / ( 1 + R * R * p2 * p2 ) );

    BW = sqrt( m_gamma ) * M /
         ( ( m_bwm * m_bwm - _p4_12.mass2() ) -
           EvtComplex( 0.0, 1.0 ) * gam * m_bwm );

    return BW;
}
