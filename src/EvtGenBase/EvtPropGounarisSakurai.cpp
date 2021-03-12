
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

#include "EvtGenBase/EvtPropGounarisSakurai.hh"

#include <cmath>

EvtPropGounarisSakurai::EvtPropGounarisSakurai( EvtDalitzPlot* dp,
                                                EvtCyclic3::Pair pair,
                                                double m0, double g0 ) :
    EvtPropagator( m0, g0 ), m_pair( pair ), m_gbase( g0 )
{
    m_dalitzSpace = dp;
    m_m1 = dp->m( EvtCyclic3::first( m_pair ) );
    m_m2 = dp->m( EvtCyclic3::second( m_pair ) );
}

EvtAmplitude<EvtPoint1D>* EvtPropGounarisSakurai::clone() const
{
    return new EvtPropGounarisSakurai( *this );
}

EvtComplex EvtPropGounarisSakurai::amplitude( const EvtPoint1D& x ) const
{
    double m = x.value();
    double s = m * m;
    double m2 = m_m0 * m_m0;
    double f_width = m_gbase;
    double f_mass = m_m0;

    double A = ( 1 + dFun( m2 ) * f_width / f_mass );
    double B = s - m2 - fsFun( s );
    //  double C = sqrt(s)*_g0;//wrong!
    double C = sqrt( m2 ) * m_g0;    //correct!
    double D = B * B + C * C;

    EvtComplex rpt( A * B / D, -A * C / D );
    return rpt;
}

//  adapted from RhoPiTools
double EvtPropGounarisSakurai::fsFun( double s ) const
{
    double m2 = m_m0 * m_m0;

    EvtTwoBodyKine vd( m_m1, m_m2, sqrt( s ) );
    EvtTwoBodyKine vR( m_m1, m_m2, m_m0 );
    double k_s = vd.p();
    double k_Am2 = vR.p();
    //
    double f = m_gbase * m2 / pow( k_Am2, 3 ) *
               ( pow( k_s, 2 ) * ( hFun( s ) - hFun( m2 ) ) +
                 ( m2 - s ) * pow( k_Am2, 2 ) * dh_dsFun( m2 ) );

    return f;
}

double EvtPropGounarisSakurai::hFun( double s ) const
{
    double sm = m_m1 + m_m2;
    double SQRTs = sqrt( s );
    EvtTwoBodyKine vd( m_m1, m_m2, sqrt( s ) );
    double k_s = vd.p();

    return 2 / EvtConst::pi * ( k_s / SQRTs ) *
           log( ( SQRTs + 2 * k_s ) / ( sm ) );
}

double EvtPropGounarisSakurai::dh_dsFun( double s ) const
{
    EvtTwoBodyKine vd( m_m1, m_m2, sqrt( s ) );
    double k_s = vd.p();

    return hFun( s ) * ( 1 / ( 8 * pow( k_s, 2 ) ) - 1 / ( 2 * s ) ) +
           1 / ( 2 * EvtConst::pi * s );
}

double EvtPropGounarisSakurai::dFun( double s ) const
{
    double sm = m_m1 + m_m2;
    double sm24 = sm * sm / 4;
    double m = sqrt( s );
    EvtTwoBodyKine vd( m_m1, m_m2, sqrt( s ) );
    double k_m2 = vd.p();
    double pi = EvtConst::pi;

    return 3.0 / pi * sm24 / pow( k_m2, 2 ) * log( ( m + 2 * k_m2 ) / sm ) +
           m / ( 2 * pi * k_m2 ) - sm24 * m / ( pi * pow( k_m2, 3 ) );
}
