
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

#include "EvtGenBase/EvtdFunctionSingle.hh"

#include "EvtGenBase/EvtPatches.hh"

#include <assert.h>
#include <iostream>
#include <math.h>
#include <stdlib.h>

EvtdFunctionSingle::EvtdFunctionSingle()
{
    m_j = 0;
    m_m1 = 0;
    m_m2 = 0;
    m_coef = nullptr;
    m_kmin = 0;
    m_kmax = 0;
}

EvtdFunctionSingle::~EvtdFunctionSingle()
{
    if ( m_coef != nullptr )
        delete[] m_coef;
}

void EvtdFunctionSingle::init( int j, int m1, int m2 )
{
    assert( abs( m2 ) >= abs( m1 ) );
    assert( m2 >= 0 );

    m_j = j;
    m_m1 = m1;
    m_m2 = m2;

    m_kmin = m_m2 - m_m1;
    m_kmax = m_j - m_m1;

    assert( m_kmin <= m_kmax );

    m_coef = new double[( m_kmax - m_kmin ) / 2 + 1];

    int k;

    for ( k = m_kmin; k <= m_kmax; k += 2 ) {
        int sign = 1;
        if ( ( k - m_m2 + m_m1 ) % 4 != 0 )
            sign = -sign;
        double tmp = fact( ( m_j + m_m2 ) / 2 ) * fact( ( m_j - m_m2 ) / 2 ) *
                     fact( ( m_j + m_m1 ) / 2 ) * fact( ( m_j - m_m1 ) / 2 );
        m_coef[( k - m_kmin ) / 2] =
            sign * sqrt( tmp ) /
            ( fact( ( m_j + m_m2 - k ) / 2 ) * fact( k / 2 ) *
              fact( ( m_j - m_m1 - k ) / 2 ) * fact( ( k - m_m2 + m_m1 ) / 2 ) );
    }
}

double EvtdFunctionSingle::d( int j, int m1, int m2, double theta )
{
    assert( j == m_j );
    UNUSED( j );
    assert( m1 == m_m1 );
    assert( m2 == m_m2 );

    double c2 = cos( 0.5 * theta );
    double s2 = sin( 0.5 * theta );

    double d = 0.0;

    int k;
    for ( k = m_kmin; k <= m_kmax; k += 2 ) {
        d += m_coef[( k - m_kmin ) / 2] *
             pow( c2, ( 2 * m_j - 2 * k + m2 - m1 ) / 2 ) *
             pow( s2, ( 2 * k - m2 + m1 ) / 2 );
    }

    return d;
}

int EvtdFunctionSingle::fact( int n )
{
    assert( n >= 0 );

    int f = 1;

    int k;
    for ( k = 2; k <= n; k++ )
        f *= k;

    return f;
}
