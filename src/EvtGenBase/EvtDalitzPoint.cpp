
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

#include "EvtGenBase/EvtDalitzPoint.hh"

#include "EvtGenBase/EvtPatches.hh"

#include <assert.h>
#include <math.h>
#include <stdio.h>
using namespace EvtCyclic3;

EvtDalitzPoint::EvtDalitzPoint() :
    m_mA( -1. ), m_mB( -1. ), m_mC( -1. ), m_qAB( -1. ), m_qBC( -1. ), m_qCA( -1. )
{
}

EvtDalitzPoint::EvtDalitzPoint( double mA, double mB, double mC, double qAB,
                                double qBC, double qCA ) :
    m_mA( mA ), m_mB( mB ), m_mC( mC ), m_qAB( qAB ), m_qBC( qBC ), m_qCA( qCA )
{
}

// Constructor from Zemach coordinates

EvtDalitzPoint::EvtDalitzPoint( double mA, double mB, double mC,
                                EvtCyclic3::Pair i, double qres, double qhel,
                                double qsum ) :
    m_mA( mA ), m_mB( mB ), m_mC( mC )
{
    double qi = qres + qsum / 3.;
    double qj = -qres / 2. + qhel + qsum / 3.;
    double qk = -qres / 2. - qhel + qsum / 3.;

    if ( i == AB ) {
        m_qAB = qi;
        m_qBC = qj;
        m_qCA = qk;
    } else if ( i == BC ) {
        m_qAB = qk;
        m_qBC = qi;
        m_qCA = qj;
    } else if ( i == CA ) {
        m_qAB = qj;
        m_qBC = qk;
        m_qCA = qi;
    }
}

EvtDalitzPoint::EvtDalitzPoint( const EvtDalitzPlot& dp,
                                const EvtDalitzCoord& x ) :
    m_mA( dp.m( A ) ), m_mB( dp.m( B ) ), m_mC( dp.m( C ) )
{
    if ( x.pair1() == AB )
        m_qAB = x.q1();
    else if ( x.pair2() == AB )
        m_qAB = x.q2();
    else
        m_qAB = dp.sum() - x.q1() - x.q2();

    if ( x.pair1() == BC )
        m_qBC = x.q1();
    else if ( x.pair2() == BC )
        m_qBC = x.q2();
    else
        m_qBC = dp.sum() - x.q1() - x.q2();

    if ( x.pair1() == CA )
        m_qCA = x.q1();
    else if ( x.pair2() == CA )
        m_qCA = x.q2();
    else
        m_qCA = dp.sum() - x.q1() - x.q2();
}

double EvtDalitzPoint::q( EvtCyclic3::Pair i ) const
{
    double ret = m_qAB;
    if ( BC == i )
        ret = m_qBC;
    else if ( CA == i )
        ret = m_qCA;

    return ret;
}

double EvtDalitzPoint::m( EvtCyclic3::Index i ) const
{
    double ret = m_mA;
    if ( B == i )
        ret = m_mB;
    else if ( C == i )
        ret = m_mC;

    return ret;
}

// Zemach variables

double EvtDalitzPoint::qres( EvtCyclic3::Pair i ) const
{
    return ( 2. * q( i ) - q( EvtCyclic3::prev( i ) ) -
             q( EvtCyclic3::next( i ) ) ) /
           3.;
}
double EvtDalitzPoint::qhel( EvtCyclic3::Pair i ) const
{
    Pair j = next( i );
    Pair k = prev( i );
    return ( q( j ) - q( k ) ) / 2.;
}
double EvtDalitzPoint::qsum() const
{
    return m_qAB + m_qBC + m_qCA;
}

double EvtDalitzPoint::qMin( EvtCyclic3::Pair i, EvtCyclic3::Pair j ) const
{
    EvtDalitzPlot dp = getDalitzPlot();
    return dp.qMin( i, j, q( j ) );
}

double EvtDalitzPoint::qMax( EvtCyclic3::Pair i, EvtCyclic3::Pair j ) const
{
    EvtDalitzPlot dp = getDalitzPlot();
    return dp.qMax( i, j, q( j ) );
}

double EvtDalitzPoint::pp( EvtCyclic3::Index i, EvtCyclic3::Index j ) const
{
    if ( i == j )
        return m( i ) * m( i );
    else
        return ( q( combine( i, j ) ) - m( i ) * m( i ) - m( j ) * m( j ) ) / 2.;
}

double EvtDalitzPoint::e( EvtCyclic3::Index i, EvtCyclic3::Pair j ) const
{
    EvtDalitzPlot dp = getDalitzPlot();
    return dp.e( i, j, q( j ) );
}

double EvtDalitzPoint::p( EvtCyclic3::Index i, EvtCyclic3::Pair j ) const
{
    EvtDalitzPlot dp = getDalitzPlot();
    return dp.p( i, j, q( j ) );
}

double EvtDalitzPoint::cosTh( EvtCyclic3::Pair pairAng,
                              EvtCyclic3::Pair pairRes ) const
{
    EvtDalitzPlot dp = getDalitzPlot();
    return dp.cosTh( pairAng, q( pairAng ), pairRes, q( pairRes ) );
}

EvtDalitzCoord EvtDalitzPoint::getDalitzPoint( EvtCyclic3::Pair i,
                                               EvtCyclic3::Pair j ) const
{
    return EvtDalitzCoord( i, q( i ), j, q( j ) );
}

EvtDalitzPlot EvtDalitzPoint::getDalitzPlot() const
{
    return EvtDalitzPlot( m_mA, m_mB, m_mC, bigM() );
}

bool EvtDalitzPoint::isValid() const
{
    // Check masses

    double M = bigM();
    if ( m_mA < 0 || m_mB < 0 || m_mC < 0 || M <= 0 )
        return false;
    if ( M < m_mA + m_mB + m_mC )
        return false;

    // Check that first coordinate is within absolute limits

    bool inside = false;
    EvtDalitzPlot dp = getDalitzPlot();

    if ( dp.qAbsMin( AB ) <= m_qAB && m_qAB <= dp.qAbsMax( AB ) )
        if ( qMin( BC, AB ) <= m_qBC && m_qBC <= qMax( BC, AB ) )
            inside = true;

    return inside;
}

double EvtDalitzPoint::bigM() const
{
    return sqrt( m_qAB + m_qBC + m_qCA - m_mA * m_mA - m_mB * m_mB - m_mC * m_mC );
}

void EvtDalitzPoint::print() const
{
    getDalitzPlot().print();
    printf( "%f %f %f\n", m_qAB, m_qBC, m_qCA );
}
