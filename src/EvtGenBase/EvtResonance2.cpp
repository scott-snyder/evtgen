
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

#include "EvtGenBase/EvtResonance2.hh"

#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtConst.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"

#include <cmath>

EvtResonance2& EvtResonance2::operator=( const EvtResonance2& n )
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
    m_invmass_angdenom = n.m_invmass_angdenom;
    m_barrier1 = n.m_barrier1;
    m_barrier2 = n.m_barrier2;
    return *this;
}

EvtResonance2::EvtResonance2( const EvtVector4R& p4_p, const EvtVector4R& p4_d1,
                              const EvtVector4R& p4_d2, double ampl,
                              double theta, double gamma, double bwm, int spin,
                              bool invmass_angdenom, double barrier1,
                              double barrier2 ) :
    m_p4_p( p4_p ),
    m_p4_d1( p4_d1 ),
    m_p4_d2( p4_d2 ),
    m_ampl( ampl ),
    m_theta( theta ),
    m_gamma( gamma ),
    m_bwm( bwm ),
    m_barrier1( barrier1 ),
    m_barrier2( barrier2 ),
    m_spin( spin ),
    m_invmass_angdenom( invmass_angdenom )
{
}

EvtComplex EvtResonance2::resAmpl() const
{
    double pi180inv = 1.0 / EvtConst::radToDegrees;

    EvtComplex ampl;
    EvtVector4R p4_d3 = m_p4_p - m_p4_d1 - m_p4_d2;

    //get cos of the angle between the daughters from their 4-momenta
    //and the 4-momentum of the parent

    //in general, EvtDecayAngle(parent, part1+part2, part1) gives the angle
    //the missing particle (not listed in the arguments) makes
    //with part2 in the rest frame of both
    //listed particles (12)

    //angle 3 makes with 2 in rest frame of 12 (CS3)
    //double cos_phi_0 = EvtDecayAngle(m_p4_p, m_p4_d1+m_p4_d2, m_p4_d1);
    //angle 3 makes with 1 in 12 is, of course, -cos_phi_0

    //first compute several quantities...follow CLEO preprint 00-23

    double mAB = ( m_p4_d1 + m_p4_d2 ).mass();
    double mBC = ( m_p4_d2 + p4_d3 ).mass();
    double mAC = ( m_p4_d1 + p4_d3 ).mass();
    double mA = m_p4_d1.mass();
    double mB = m_p4_d2.mass();
    double mD = m_p4_p.mass();
    double mC = p4_d3.mass();

    double mR = m_bwm;
    double gammaR = m_gamma;
    double mdenom = m_invmass_angdenom ? mAB : mR;
    double pAB = sqrt( ( ( ( mAB * mAB - mA * mA - mB * mB ) *
                           ( mAB * mAB - mA * mA - mB * mB ) / 4.0 ) -
                         mA * mA * mB * mB ) /
                       ( mAB * mAB ) );
    double pR = sqrt( ( ( ( mR * mR - mA * mA - mB * mB ) *
                          ( mR * mR - mA * mA - mB * mB ) / 4.0 ) -
                        mA * mA * mB * mB ) /
                      ( mR * mR ) );

    double pD = ( ( ( mD * mD - mR * mR - mC * mC ) *
                    ( mD * mD - mR * mR - mC * mC ) / 4.0 ) -
                  mR * mR * mC * mC ) /
                ( mD * mD );
    if ( pD > 0 ) {
        pD = sqrt( pD );
    } else {
        pD = 0;
    }
    double pDAB = sqrt( ( ( ( mD * mD - mAB * mAB - mC * mC ) *
                            ( mD * mD - mAB * mAB - mC * mC ) / 4.0 ) -
                          mAB * mAB * mC * mC ) /
                        ( mD * mD ) );

    double fR = 1;
    double fD = 1;
    int power = 0;
    switch ( m_spin ) {
        case 0:
            fR = 1.0;
            fD = 1.0;
            power = 1;
            break;
        case 1:
            fR = sqrt( 1.0 + m_barrier1 * m_barrier1 * pR * pR ) /
                 sqrt( 1.0 + m_barrier1 * m_barrier1 * pAB * pAB );
            fD = sqrt( 1.0 + m_barrier2 * m_barrier2 * pD * pD ) /
                 sqrt( 1.0 + m_barrier2 * m_barrier2 * pDAB * pDAB );
            power = 3;
            break;
        case 2:
            fR = sqrt( ( 9 + 3 * pow( ( m_barrier1 * pR ), 2 ) +
                         pow( ( m_barrier1 * pR ), 4 ) ) /
                       ( 9 + 3 * pow( ( m_barrier1 * pAB ), 2 ) +
                         pow( ( m_barrier1 * pAB ), 4 ) ) );
            fD = sqrt( ( 9 + 3 * pow( ( m_barrier2 * pD ), 2 ) +
                         pow( ( m_barrier2 * pD ), 4 ) ) /
                       ( 9 + 3 * pow( ( m_barrier2 * pDAB ), 2 ) +
                         pow( ( m_barrier2 * pDAB ), 4 ) ) );
            power = 5;
            break;
        default:
            EvtGenReport( EVTGEN_INFO, "EvtGen" )
                << "Incorrect spin in EvtResonance2.cc\n";
    }

    double gammaAB = gammaR * pow( pAB / pR, power ) * ( mR / mAB ) * fR * fR;
    switch ( m_spin ) {
        case 0:
            ampl = m_ampl *
                   EvtComplex( cos( m_theta * pi180inv ),
                               sin( m_theta * pi180inv ) ) *
                   fR * fD /
                   ( mR * mR - mAB * mAB - EvtComplex( 0.0, mR * gammaAB ) );
            break;
        case 1:
            ampl = m_ampl *
                   EvtComplex( cos( m_theta * pi180inv ),
                               sin( m_theta * pi180inv ) ) *
                   ( fR * fD *
                     ( mAC * mAC - mBC * mBC +
                       ( ( mD * mD - mC * mC ) * ( mB * mB - mA * mA ) /
                         ( mdenom * mdenom ) ) ) /
                     ( mR * mR - mAB * mAB - EvtComplex( 0.0, mR * gammaAB ) ) );
            break;
        case 2:
            ampl = m_ampl *
                   EvtComplex( cos( m_theta * pi180inv ),
                               sin( m_theta * pi180inv ) ) *
                   fR * fD /
                   ( mR * mR - mAB * mAB - EvtComplex( 0.0, mR * gammaAB ) ) *
                   ( pow( ( mBC * mBC - mAC * mAC +
                            ( mD * mD - mC * mC ) * ( mA * mA - mB * mB ) /
                                ( mdenom * mdenom ) ),
                          2 ) -
                     ( 1.0 / 3.0 ) *
                         ( mAB * mAB - 2 * mD * mD - 2 * mC * mC +
                           pow( ( mD * mD - mC * mC ) / mdenom, 2 ) ) *
                         ( mAB * mAB - 2 * mA * mA - 2 * mB * mB +
                           pow( ( mA * mA - mB * mB ) / mdenom, 2 ) ) );
            break;

        default:
            EvtGenReport( EVTGEN_INFO, "EvtGen" )
                << "Incorrect spin in EvtResonance2.cc\n";
    }

    return ampl;
}
