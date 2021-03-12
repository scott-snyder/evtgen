
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

#include "EvtGenBase/EvtValError.hh"

#include <cassert>
#include <cmath>
#include <iostream>

using std::endl;
using std::ostream;

EvtValError::EvtValError() :
    m_valKnown( 0 ), m_val( 0. ), m_errKnown( 0 ), m_err( 0. )
{
}

EvtValError::EvtValError( double val ) :
    m_valKnown( 1 ), m_val( val ), m_errKnown( 0 ), m_err( 0. )
{
}

EvtValError::EvtValError( double val, double err ) :
    m_valKnown( 1 ), m_val( val ), m_errKnown( 1 ), m_err( err )
{
}

EvtValError::EvtValError( const EvtValError& other ) :
    m_valKnown( other.m_valKnown ),
    m_val( other.m_val ),
    m_errKnown( other.m_errKnown ),
    m_err( other.m_err )
{
}

double EvtValError::prec() const
{
    assert( m_valKnown && m_errKnown );
    return ( m_val != 0 ) ? m_err / m_val : 0;
}

void EvtValError::operator=( const EvtValError& other )
{
    m_valKnown = other.m_valKnown;
    m_val = other.m_val;
    m_errKnown = other.m_errKnown;
    m_err = other.m_err;
}

void EvtValError::operator*=( const EvtValError& other )
{
    assert( m_valKnown && other.m_valKnown );

    // Relative errors add in quadrature
    if ( m_errKnown && other.m_errKnown )
        m_err = m_val * other.m_val *
                sqrt( prec() * prec() + other.prec() * other.prec() );
    else
        m_errKnown = 0;

    // Modify the value
    m_val *= other.m_val;
}

void EvtValError::operator/=( const EvtValError& other )
{
    assert( m_valKnown && other.m_valKnown && other.m_val != 0. );

    // Relative errors add in quadrature
    if ( m_errKnown && other.m_errKnown )
        m_err = m_val / other.m_val *
                sqrt( prec() * prec() + other.prec() * other.prec() );
    else
        m_errKnown = 0;

    // Modify the value
    m_val /= other.m_val;
}

void EvtValError::print( ostream& os ) const
{
    if ( m_valKnown )
        os << m_val;
    else
        os << "Undef";
    os << " +/- ";
    if ( m_errKnown )
        os << m_err;
    else
        os << "Undef";
    os << endl;
}

void EvtValError::operator+=( const EvtValError& other )
{
    assert( m_valKnown );
    assert( other.m_valKnown );
    m_val += other.m_val;

    // add errors in quadrature

    if ( m_errKnown && other.m_errKnown ) {
        m_err = sqrt( m_err * m_err + other.m_err * other.m_err );
    } else {
        m_errKnown = 0;
    }
}

void EvtValError::operator*=( double c )
{
    assert( m_valKnown );
    m_val *= c;
    if ( m_errKnown )
        m_err *= c;
}

EvtValError operator*( const EvtValError& x1, const EvtValError& x2 )
{
    EvtValError ret( x1 );
    ret *= x2;
    return ret;
}

EvtValError operator/( const EvtValError& x1, const EvtValError& x2 )
{
    EvtValError ret( x1 );
    ret /= x2;
    return ret;
}

EvtValError operator+( const EvtValError& x1, const EvtValError& x2 )
{
    EvtValError ret( x1 );
    ret += x2;
    return ret;
}

EvtValError operator*( const EvtValError& x, double c )
{
    EvtValError ret( x );
    ret *= c;
    return ret;
}

EvtValError operator*( double c, const EvtValError& x )
{
    EvtValError ret( x );
    ret *= c;
    return ret;
}

ostream& operator<<( ostream& os, const EvtValError& other )
{
    other.print( os );
    return os;
}
