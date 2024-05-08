
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

#include "EvtGenBase/EvtFlatte.hh"

#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtConst.hh"
#include "EvtGenBase/EvtKine.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtVector4R.hh"

#include <math.h>

//operator

EvtFlatte& EvtFlatte::operator=( const EvtFlatte& n )
{
    if ( &n == this )
        return *this;
    m_p4_p = n.m_p4_p;
    m_p4_d1 = n.m_p4_d1;
    m_p4_d2 = n.m_p4_d2;
    m_ampl = n.m_ampl;
    m_theta = n.m_theta;
    m_mass = n.m_mass;
    m_params = n.m_params;
    return *this;
}

//constructor

EvtFlatte::EvtFlatte( const EvtVector4R& p4_p, const EvtVector4R& p4_d1,
                      const EvtVector4R& p4_d2, double ampl, double theta,
                      double mass, vector<EvtFlatteParam>& params ) :
    m_p4_p( p4_p ),
    m_p4_d1( p4_d1 ),
    m_p4_d2( p4_d2 ),
    m_ampl( ampl ),
    m_theta( theta ),
    m_mass( mass ),
    m_params( params )
{
}

//amplitude function

EvtComplex EvtFlatte::resAmpl()
{
    double pi180inv = 1.0 / EvtConst::radToDegrees;

    // SCALARS ONLY
    double mR = ( m_p4_d1 + m_p4_d2 ).mass();

    EvtComplex w;

    for ( vector<EvtFlatteParam>::const_iterator param = m_params.begin();
          param != m_params.end(); ++param ) {
        double m1 = ( *param ).m1();
        double m2 = ( *param ).m2();
        double g = ( *param ).g();
        w += ( g * g *
               sqrtCplx( ( 1 - ( ( m1 - m2 ) * ( m1 - m2 ) ) / ( mR * mR ) ) *
                         ( 1 - ( ( m1 + m2 ) * ( m1 + m2 ) ) / ( mR * mR ) ) ) );
    }

    EvtComplex denom = m_mass * m_mass - mR * mR - EvtComplex( 0, 1 ) * w;
    EvtComplex ampl = m_ampl *
                      EvtComplex( cos( m_theta * pi180inv ),
                                  sin( m_theta * pi180inv ) ) /
                      denom;
    return ampl;
}
