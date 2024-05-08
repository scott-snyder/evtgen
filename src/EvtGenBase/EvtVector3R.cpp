
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

#include "EvtGenBase/EvtVector3R.hh"

#include "EvtGenBase/EvtPatches.hh"

#include <iostream>
#include <math.h>
using std::ostream;

EvtVector3R::EvtVector3R()
{
    m_v[0] = m_v[1] = m_v[2] = 0.0;
}

EvtVector3R::EvtVector3R( double x, double y, double z )
{
    m_v[0] = x;
    m_v[1] = y;
    m_v[2] = z;
}

EvtVector3R rotateEuler( const EvtVector3R& v, double alpha, double beta,
                         double gamma )
{
    EvtVector3R tmp( v );
    tmp.applyRotateEuler( alpha, beta, gamma );
    return tmp;
}

void EvtVector3R::applyRotateEuler( double phi, double theta, double ksi )
{
    double temp[3];
    double sp, st, sk, cp, ct, ck;

    sp = sin( phi );
    st = sin( theta );
    sk = sin( ksi );
    cp = cos( phi );
    ct = cos( theta );
    ck = cos( ksi );

    temp[0] = ( ck * ct * cp - sk * sp ) * m_v[0] +
              ( -sk * ct * cp - ck * sp ) * m_v[1] + st * cp * m_v[2];
    temp[1] = ( ck * ct * sp + sk * cp ) * m_v[0] +
              ( -sk * ct * sp + ck * cp ) * m_v[1] + st * sp * m_v[2];
    temp[2] = -ck * st * m_v[0] + sk * st * m_v[1] + ct * m_v[2];

    m_v[0] = temp[0];
    m_v[1] = temp[1];
    m_v[2] = temp[2];
}

ostream& operator<<( ostream& s, const EvtVector3R& v )
{
    s << "(" << v.m_v[0] << "," << v.m_v[1] << "," << v.m_v[2] << ")";

    return s;
}

EvtVector3R cross( const EvtVector3R& p1, const EvtVector3R& p2 )
{
    //Calcs the cross product.  Added by djl on July 27, 1995.
    //Modified for real vectros by ryd Aug 28-96

    return EvtVector3R( p1.m_v[1] * p2.m_v[2] - p1.m_v[2] * p2.m_v[1],
                        p1.m_v[2] * p2.m_v[0] - p1.m_v[0] * p2.m_v[2],
                        p1.m_v[0] * p2.m_v[1] - p1.m_v[1] * p2.m_v[0] );
}

double EvtVector3R::d3mag() const
{
    // Returns the 3 momentum mag
    double temp;

    temp = m_v[0] * m_v[0] + m_v[1] * m_v[1] + m_v[2] * m_v[2];
    temp = sqrt( temp );

    return temp;
}

double EvtVector3R::dot( const EvtVector3R& p2 )
{
    double temp;

    temp = m_v[0] * p2.m_v[0];
    temp += m_v[1] * p2.m_v[1];
    temp += m_v[2] * p2.m_v[2];

    return temp;
}
