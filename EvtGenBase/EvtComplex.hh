
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

#ifndef EVTCOMPLEX_HH
#define EVTCOMPLEX_HH

#include "EvtGenBase/EvtConst.hh"

#include <iostream>
#include <math.h>

class EvtComplex {
    inline friend EvtComplex operator*( double d, const EvtComplex& c );
    inline friend EvtComplex operator*( const EvtComplex& c, double d );
    inline friend EvtComplex operator/( const EvtComplex& c, double d );
    inline friend EvtComplex operator/( double d, const EvtComplex& c );
    inline friend EvtComplex operator*( const EvtComplex& c1,
                                        const EvtComplex& c2 );
    inline friend EvtComplex operator/( const EvtComplex& c1,
                                        const EvtComplex& c2 );
    inline friend EvtComplex operator+( const EvtComplex& c1,
                                        const EvtComplex& c2 );
    inline friend EvtComplex operator-( const EvtComplex& c1,
                                        const EvtComplex& c2 );
    inline friend EvtComplex operator-( const EvtComplex& c );
    inline friend EvtComplex conj( const EvtComplex& c );
    inline friend double abs( const EvtComplex& c );
    inline friend double abs2( const EvtComplex& c );
    inline friend double arg( const EvtComplex& c );
    inline friend double real( const EvtComplex& c );
    inline friend double imag( const EvtComplex& c );
    inline friend EvtComplex exp( const EvtComplex& c );
    friend std::ostream& operator<<( std::ostream& s, const EvtComplex& c );

  public:
    EvtComplex() : m_rpart( 0.0 ), m_ipart( 0.0 ) {}
    EvtComplex( double rpart, double ipart = 0.0 ) :
        m_rpart( rpart ), m_ipart( ipart )
    {
    }
    EvtComplex( const EvtComplex& c ) :
        m_rpart( c.m_rpart ), m_ipart( c.m_ipart )
    {
    }
    inline EvtComplex& operator*=( double d );
    inline EvtComplex& operator/=( double d );
    EvtComplex& operator*=( EvtComplex c );
    EvtComplex& operator/=( EvtComplex c );
    inline EvtComplex& operator=( const EvtComplex& c );
    inline EvtComplex& operator+=( const EvtComplex& c );
    inline EvtComplex& operator-=( const EvtComplex& c );
    inline EvtComplex& operator+=( double d );
    inline EvtComplex& operator-=( double d );
    inline int operator==( const EvtComplex c );
    inline int operator!=( const EvtComplex c );

  private:
    double m_rpart, m_ipart;
};

typedef EvtComplex* EvtComplexPtr;
typedef EvtComplexPtr* EvtComplexPtrPtr;
typedef EvtComplexPtrPtr* EvtComplexPtrPtrPtr;

EvtComplex& EvtComplex::operator=( const EvtComplex& c )
{
    m_rpart = c.m_rpart;
    m_ipart = c.m_ipart;

    return *this;
}

EvtComplex& EvtComplex::operator+=( const EvtComplex& c )
{
    m_rpart += c.m_rpart;
    m_ipart += c.m_ipart;

    return *this;
}

EvtComplex& EvtComplex::operator-=( const EvtComplex& c )
{
    m_rpart -= c.m_rpart;
    m_ipart -= c.m_ipart;

    return *this;
}

EvtComplex& EvtComplex::operator+=( double d )
{
    m_rpart += d;

    return *this;
}

EvtComplex& EvtComplex::operator-=( double d )
{
    m_rpart -= d;

    return *this;
}

EvtComplex operator*( double d, const EvtComplex& c )
{
    return EvtComplex( c.m_rpart * d, c.m_ipart * d );
}

EvtComplex operator*( const EvtComplex& c, double d )
{
    return EvtComplex( c.m_rpart * d, c.m_ipart * d );
}

EvtComplex operator/( const EvtComplex& c, double d )
{
    return EvtComplex( c.m_rpart / d, c.m_ipart / d );
}

EvtComplex& EvtComplex::operator*=( double d )
{
    m_rpart *= d;
    m_ipart *= d;

    return *this;
}

EvtComplex& EvtComplex::operator/=( double d )
{
    m_rpart /= d;
    m_ipart /= d;

    return *this;
}

EvtComplex operator/( double d, const EvtComplex& c )
{
    double Num = d / ( c.m_rpart * c.m_rpart + c.m_ipart * c.m_ipart );

    return EvtComplex( Num * c.m_rpart, -Num * c.m_ipart );
}

EvtComplex operator/( const EvtComplex& c1, const EvtComplex& c2 )
{
    double inv = 1.0 / ( c2.m_rpart * c2.m_rpart + c2.m_ipart * c2.m_ipart );

    return EvtComplex(
        inv * ( c1.m_rpart * c2.m_rpart + c1.m_ipart * c2.m_ipart ),
        inv * ( c1.m_ipart * c2.m_rpart - c1.m_rpart * c2.m_ipart ) );
}

EvtComplex operator*( const EvtComplex& c1, const EvtComplex& c2 )
{
    return EvtComplex( c1.m_rpart * c2.m_rpart - c1.m_ipart * c2.m_ipart,
                       c1.m_rpart * c2.m_ipart + c1.m_ipart * c2.m_rpart );
}

EvtComplex operator-( const EvtComplex& c1, const EvtComplex& c2 )
{
    return EvtComplex( c1.m_rpart - c2.m_rpart, c1.m_ipart - c2.m_ipart );
}

EvtComplex operator+( const EvtComplex& c1, const EvtComplex& c2 )
{
    return EvtComplex( c1.m_rpart + c2.m_rpart, c1.m_ipart + c2.m_ipart );
}

int EvtComplex::operator==( const EvtComplex c )
{
    return m_rpart == c.m_rpart && m_ipart == c.m_ipart;
}

int EvtComplex::operator!=( const EvtComplex c )
{
    return m_rpart != c.m_rpart || m_ipart != c.m_ipart;
}

EvtComplex operator-( const EvtComplex& c )
{
    return EvtComplex( -c.m_rpart, -c.m_ipart );
}

EvtComplex conj( const EvtComplex& c )
{
    return EvtComplex( c.m_rpart, -c.m_ipart );
}

double abs( const EvtComplex& c )
{
    double c2 = c.m_rpart * c.m_rpart + c.m_ipart * c.m_ipart;
    if ( c2 <= 0.0 )
        return 0.0;
    return sqrt( c2 );
}

double abs2( const EvtComplex& c )
{
    return c.m_rpart * c.m_rpart + c.m_ipart * c.m_ipart;
}

double arg( const EvtComplex& c )
{
    if ( ( c.m_rpart == 0 ) && ( c.m_ipart == 0 ) ) {
        return 0.0;
    }
    if ( c.m_rpart == 0 ) {
        if ( c.m_ipart > 0 ) {
            return EvtConst::pi / 2;
        } else {
            return -EvtConst::pi / 2;
        }
    } else {
        return atan2( c.m_ipart, c.m_rpart );
    }
}

double real( const EvtComplex& c )
{
    return c.m_rpart;
}

double imag( const EvtComplex& c )
{
    return c.m_ipart;
}

EvtComplex exp( const EvtComplex& c )
{
    return exp( c.m_rpart ) * EvtComplex( cos( c.m_ipart ), sin( c.m_ipart ) );
}

#endif
