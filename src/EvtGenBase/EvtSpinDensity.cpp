
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

#include "EvtGenBase/EvtSpinDensity.hh"

#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtReport.hh"

#include <cassert>
#include <cstdlib>
#include <iostream>

using std::endl;
using std::ostream;

EvtSpinDensity::EvtSpinDensity( const EvtSpinDensity& density )
{
    m_dim = 0;
    m_rho = nullptr;

    int i, j;
    setDim( density.m_dim );

    for ( i = 0; i < m_dim; i++ ) {
        for ( j = 0; j < m_dim; j++ ) {
            m_rho[i][j] = density.m_rho[i][j];
        }
    }
}

EvtSpinDensity& EvtSpinDensity::operator=( const EvtSpinDensity& density )
{
    int i, j;
    setDim( density.m_dim );

    for ( i = 0; i < m_dim; i++ ) {
        for ( j = 0; j < m_dim; j++ ) {
            m_rho[i][j] = density.m_rho[i][j];
        }
    }

    return *this;
}

EvtSpinDensity::~EvtSpinDensity()
{
    if ( m_dim != 0 ) {
        int i;
        for ( i = 0; i < m_dim; i++ )
            delete[] m_rho[i];
    }

    delete[] m_rho;
}

EvtSpinDensity::EvtSpinDensity()
{
    m_dim = 0;
    m_rho = nullptr;
}

void EvtSpinDensity::setDim( int n )
{
    if ( m_dim == n )
        return;
    if ( m_dim != 0 ) {
        int i;
        for ( i = 0; i < m_dim; i++ )
            delete[] m_rho[i];
        delete[] m_rho;
        m_rho = nullptr;
        m_dim = 0;
    }
    if ( n == 0 )
        return;
    m_dim = n;
    m_rho = new EvtComplexPtr[n];
    int i;
    for ( i = 0; i < n; i++ ) {
        m_rho[i] = new EvtComplex[n];
    }
}

int EvtSpinDensity::getDim() const
{
    return m_dim;
}

void EvtSpinDensity::set( int i, int j, const EvtComplex& rhoij )
{
    assert( i < m_dim && j < m_dim );
    m_rho[i][j] = rhoij;
}

const EvtComplex& EvtSpinDensity::get( int i, int j ) const
{
    assert( i < m_dim && j < m_dim );
    return m_rho[i][j];
}

void EvtSpinDensity::setDiag( int n )
{
    setDim( n );
    int i, j;

    for ( i = 0; i < n; i++ ) {
        for ( j = 0; j < n; j++ ) {
            m_rho[i][j] = EvtComplex( 0.0 );
        }
        m_rho[i][i] = EvtComplex( 1.0 );
    }
}

double EvtSpinDensity::normalizedProb( const EvtSpinDensity& d )
{
    int i, j;
    EvtComplex prob( 0.0, 0.0 );
    double norm = 0.0;

    if ( m_dim != d.m_dim ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Not matching dimensions in NormalizedProb" << endl;
        ::abort();
    }

    for ( i = 0; i < m_dim; i++ ) {
        norm += real( m_rho[i][i] );
        for ( j = 0; j < m_dim; j++ ) {
            prob += m_rho[i][j] * d.m_rho[i][j];
        }
    }

    if ( imag( prob ) > 0.00000001 * real( prob ) ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Imaginary probability:" << prob << " " << norm << endl;
    }
    if ( real( prob ) < 0.0 ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Negative probability:" << prob << " " << norm << endl;
    }

    return real( prob ) / norm;
}

int EvtSpinDensity::check()
{
    if ( m_dim < 1 ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "dim=" << m_dim << "in SpinDensity::Check" << endl;
    }

    int i, j;

    double trace( 0.0 );

    for ( i = 0; i < m_dim; i++ ) {
        trace += abs( m_rho[i][i] );
    }

    for ( i = 0; i < m_dim; i++ ) {
        if ( real( m_rho[i][i] ) < 0.0 )
            return 0;
        if ( imag( m_rho[i][i] ) * 1000000.0 > trace ) {
            EvtGenReport( EVTGEN_INFO, "EvtGen" ) << *this << endl;
            EvtGenReport( EVTGEN_INFO, "EvtGen" ) << trace << endl;
            EvtGenReport( EVTGEN_INFO, "EvtGen" ) << "Failing 1" << endl;
            return 0;
        }
    }

    for ( i = 0; i < m_dim; i++ ) {
        for ( j = i + 1; j < m_dim; j++ ) {
            if ( fabs( real( m_rho[i][j] - m_rho[j][i] ) ) >
                 0.00000001 * ( abs( m_rho[i][i] ) + abs( m_rho[j][j] ) ) ) {
                EvtGenReport( EVTGEN_INFO, "EvtGen" ) << "Failing 2" << endl;
                return 0;
            }
            if ( fabs( imag( m_rho[i][j] + m_rho[j][i] ) ) >
                 0.00000001 * ( abs( m_rho[i][i] ) + abs( m_rho[j][j] ) ) ) {
                EvtGenReport( EVTGEN_INFO, "EvtGen" ) << "Failing 3" << endl;
                return 0;
            }
        }
    }

    return 1;
}

ostream& operator<<( ostream& s, const EvtSpinDensity& d )
{
    int i, j;

    s << endl;
    s << "Dimension:" << d.m_dim << endl;

    for ( i = 0; i < d.m_dim; i++ ) {
        for ( j = 0; j < d.m_dim; j++ ) {
            s << d.m_rho[i][j] << " ";
        }
        s << endl;
    }

    return s;
}
