
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

#include "EvtGenBase/EvtLASSAmp.hh"

#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtCyclic3.hh"
#include "EvtGenBase/EvtDalitzCoord.hh"
#include "EvtGenBase/EvtdFunction.hh"

#include <cassert>
#include <cmath>
#include <iostream>

using EvtCyclic3::Index;
using EvtCyclic3::Pair;
using std::endl;

EvtLASSAmp::EvtLASSAmp( EvtDalitzPlot* dp, EvtCyclic3::Pair pair, double m0,
                        double g0, double a, double r, double cutoff,
                        std::string subtype ) :
    EvtAmplitude<EvtDalitzPoint>(),
    m_pair( pair ),
    m_m0( m0 ),
    m_g0( g0 ),
    m_r( r ),
    m_a( a ),
    m_cutoff( cutoff ),
    m_subtype( subtype )
{
    m_dalitzSpace = dp;
    double ma = dp->m( first( pair ) );
    double mb = dp->m( second( pair ) );
    double E0a = 0.5 * ( m_m0 * m_m0 + ma * ma - mb * mb ) / m_m0;
    m_q0 = E0a * E0a - ma * ma;
    assert( m_q0 > 0 );
    m_q0 = sqrt( m_q0 );
}

EvtComplex EvtLASSAmp::amplitude( const EvtDalitzPoint& dalitzPoint ) const
{
    /*

    Parameterization of Kpi S-wave using LASS scattering data.
    - Nucl.Phys.B296, 493 (1988)
    - W.Dunwoodie,http://www.slac.stanford.edu/~wmd/kpi_swave/kpi_swave_fit.note

            m                                     m0^2*Gamma0/q0
    ----------------- + exp(2*i*delta) * --------------------------------
    q*cot(delta)-i*q                     m0^2-m^2 - i*m0*Gamma0*q/m*m0/q0


    where q = momentum of K or pi in Kpi system

          q*cot(delta) = 1/ a   + 1/2 * [ r * q**2 ]

	  a = scattering length

	  r = effective range

  */

    double s = dalitzPoint.q( m_pair );
    double m = sqrt( s );
    double q = dalitzPoint.p( first( m_pair ), m_pair );

    // elastic scattering
    double qcotd = 1. / m_a + 0.5 * m_r * q * q;
    EvtComplex lass_elastic = m < m_cutoff ? m / ( qcotd - EvtComplex( 0, q ) )
                                           : 0;

    // relative phase
    double cosd = 1;
    double sind = 0;
    if ( q > 0 ) {
        cosd = qcotd * qcotd / ( q * q );
        cosd = sqrt( cosd / ( 1 + cosd ) );
        sind = sqrt( 1 - cosd * cosd );
    }
    EvtComplex lass_phase( cosd, sind );
    lass_phase *= lass_phase;

    // K*(1430)
    double gamma = m_g0 * q / m * m_m0 / m_q0;
    EvtComplex lass_Kstar = ( m_m0 * m_m0 ) * ( m_g0 / m_q0 ) /
                            ( m_m0 * m_m0 - m * m -
                              EvtComplex( 0., m_m0 * gamma ) );

    EvtComplex theAmplitude( 0.0, 0.0 );

    if ( m_subtype == "LASS_ELASTIC" ) {
        theAmplitude = lass_elastic;

    } else if ( m_subtype == "LASS_RESONANT" ) {
        theAmplitude = lass_phase * lass_Kstar;

    } else {
        theAmplitude = lass_phase * lass_Kstar + lass_elastic;
    }

    return theAmplitude;
}
