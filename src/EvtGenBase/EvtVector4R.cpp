
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

#include "EvtGenBase/EvtVector4R.hh"

#include "EvtGenBase/EvtTensor4C.hh"
#include "EvtGenBase/EvtVector3R.hh"
#include "EvtGenBase/EvtVector4C.hh"

#include <cmath>
#include <iostream>

using std::ostream;

EvtVector4R::EvtVector4R()
{
    m_v[0] = 0.0;
    m_v[1] = 0.0;
    m_v[2] = 0.0;
    m_v[3] = 0.0;
}

EvtVector4R::EvtVector4R( double e, double p1, double p2, double p3 )
{
    m_v[0] = e;
    m_v[1] = p1;
    m_v[2] = p2;
    m_v[3] = p3;
}

double EvtVector4R::mass() const
{
    double m2 = m_v[0] * m_v[0] - m_v[1] * m_v[1] - m_v[2] * m_v[2] -
                m_v[3] * m_v[3];

    if ( m2 > 0.0 ) {
        return sqrt( m2 );
    } else {
        return 0.0;
    }
}

EvtVector4R rotateEuler( const EvtVector4R& rs, double alpha, double beta,
                         double gamma )
{
    EvtVector4R tmp( rs );
    tmp.applyRotateEuler( alpha, beta, gamma );
    return tmp;
}

EvtVector4R boostTo( const EvtVector4R& rs, const EvtVector4R& p4, bool inverse )
{
    EvtVector4R tmp( rs );
    tmp.applyBoostTo( p4, inverse );
    return tmp;
}

EvtVector4R boostTo( const EvtVector4R& rs, const EvtVector3R& boost,
                     bool inverse )
{
    EvtVector4R tmp( rs );
    tmp.applyBoostTo( boost, inverse );
    return tmp;
}

void EvtVector4R::applyRotateEuler( double phi, double theta, double ksi )
{
    double sp = sin( phi );
    double st = sin( theta );
    double sk = sin( ksi );
    double cp = cos( phi );
    double ct = cos( theta );
    double ck = cos( ksi );

    double x = ( ck * ct * cp - sk * sp ) * m_v[1] +
               ( -sk * ct * cp - ck * sp ) * m_v[2] + st * cp * m_v[3];
    double y = ( ck * ct * sp + sk * cp ) * m_v[1] +
               ( -sk * ct * sp + ck * cp ) * m_v[2] + st * sp * m_v[3];
    double z = -ck * st * m_v[1] + sk * st * m_v[2] + ct * m_v[3];

    m_v[1] = x;
    m_v[2] = y;
    m_v[3] = z;
}

ostream& operator<<( ostream& s, const EvtVector4R& v )
{
    s << "(" << v.m_v[0] << "," << v.m_v[1] << "," << v.m_v[2] << ","
      << v.m_v[3] << ")";

    return s;
}

void EvtVector4R::applyBoostTo( const EvtVector4R& p4, bool inverse )
{
    double e = p4.get( 0 );

    EvtVector3R boost( p4.get( 1 ) / e, p4.get( 2 ) / e, p4.get( 3 ) / e );

    applyBoostTo( boost, inverse );

    return;
}

void EvtVector4R::applyBoostTo( const EvtVector3R& boost, bool inverse )
{
    double bx, by, bz, gamma, b2;

    bx = boost.get( 0 );
    by = boost.get( 1 );
    bz = boost.get( 2 );

    double bxx = bx * bx;
    double byy = by * by;
    double bzz = bz * bz;

    b2 = bxx + byy + bzz;

    if ( b2 > 0.0 && b2 < 1.0 ) {
        gamma = 1.0 / sqrt( 1.0 - b2 );

        double gb2 = ( gamma - 1.0 ) / b2;

        double gb2xy = gb2 * bx * by;
        double gb2xz = gb2 * bx * bz;
        double gb2yz = gb2 * by * bz;

        double gbx = gamma * bx;
        double gby = gamma * by;
        double gbz = gamma * bz;

        double e2 = m_v[0];
        double px2 = m_v[1];
        double py2 = m_v[2];
        double pz2 = m_v[3];

        if ( inverse ) {
            m_v[0] = gamma * e2 - gbx * px2 - gby * py2 - gbz * pz2;

            m_v[1] = -gbx * e2 + gb2 * bxx * px2 + px2 + gb2xy * py2 +
                     gb2xz * pz2;

            m_v[2] = -gby * e2 + gb2 * byy * py2 + py2 + gb2xy * px2 +
                     gb2yz * pz2;

            m_v[3] = -gbz * e2 + gb2 * bzz * pz2 + pz2 + gb2yz * py2 +
                     gb2xz * px2;
        } else {
            m_v[0] = gamma * e2 + gbx * px2 + gby * py2 + gbz * pz2;

            m_v[1] = gbx * e2 + gb2 * bxx * px2 + px2 + gb2xy * py2 + gb2xz * pz2;

            m_v[2] = gby * e2 + gb2 * byy * py2 + py2 + gb2xy * px2 + gb2yz * pz2;

            m_v[3] = gbz * e2 + gb2 * bzz * pz2 + pz2 + gb2yz * py2 + gb2xz * px2;
        }
    }
}

EvtVector4R EvtVector4R::cross( const EvtVector4R& p2 ) const
{
    //Calcs the cross product.  Added by djl on July 27, 1995.
    //Modified for real vectros by ryd Aug 28-96

    EvtVector4R temp;

    temp.m_v[0] = 0.0;
    temp.m_v[1] = m_v[2] * p2.m_v[3] - m_v[3] * p2.m_v[2];
    temp.m_v[2] = m_v[3] * p2.m_v[1] - m_v[1] * p2.m_v[3];
    temp.m_v[3] = m_v[1] * p2.m_v[2] - m_v[2] * p2.m_v[1];

    return temp;
}

double EvtVector4R::d3mag() const

// returns the 3 momentum mag.
{
    double temp;

    temp = m_v[1] * m_v[1] + m_v[2] * m_v[2] + m_v[3] * m_v[3];

    temp = sqrt( temp );

    return temp;
}    // r3mag

double EvtVector4R::dot( const EvtVector4R& p2 ) const
{
    //Returns the dot product of the 3 momentum.  Added by
    //djl on July 27, 1995.  for real!!!

    double temp;

    temp = m_v[1] * p2.m_v[1];
    temp += m_v[2] * p2.m_v[2];
    temp += m_v[3] * p2.m_v[3];

    return temp;

}    //dot

// Functions below added by AJB

// Calculate ( \vec{p1} cross \vec{p2} ) \cdot \vec{p3} in rest frame of object
double EvtVector4R::scalartripler3( const EvtVector4R& p1, const EvtVector4R& p2,
                                    const EvtVector4R& p3 ) const
{
    EvtVector4C lc = dual( EvtGenFunctions::directProd( *this, p1 ) ).cont2( p2 );
    EvtVector4R l( real( lc.get( 0 ) ), real( lc.get( 1 ) ),
                   real( lc.get( 2 ) ), real( lc.get( 3 ) ) );

    return -1.0 / mass() * ( l * p3 );
}

// Calculate the 3-d dot product of 4-vectors p1 and p2 in the rest frame of
// 4-vector p0
double EvtVector4R::dotr3( const EvtVector4R& p1, const EvtVector4R& p2 ) const
{
    return 1 / mass2() * ( ( *this ) * p1 ) * ( ( *this ) * p2 ) - p1 * p2;
}

// Calculate the 3-d magnitude squared of 4-vector p1 in the rest frame of
// 4-vector p0
double EvtVector4R::mag2r3( const EvtVector4R& p1 ) const
{
    return Square( ( *this ) * p1 ) / mass2() - p1.mass2();
}

// Calculate the 3-d magnitude 4-vector p1 in the rest frame of 4-vector p0.
double EvtVector4R::magr3( const EvtVector4R& p1 ) const
{
    return sqrt( mag2r3( p1 ) );
}
