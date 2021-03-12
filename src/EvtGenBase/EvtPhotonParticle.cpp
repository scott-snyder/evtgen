
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

#include "EvtGenBase/EvtPhotonParticle.hh"

#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtVector4C.hh"

#include <cmath>
#include <iostream>

using std::endl;

void EvtPhotonParticle::init( EvtId part_n, const EvtVector4R& p4 )
{
    init( part_n, p4.get( 0 ), p4.get( 1 ), p4.get( 2 ), p4.get( 3 ) );
}

void EvtPhotonParticle::init( EvtId part_n, double e, double px, double py,
                              double pz )
{
    m_validP4 = true;
    setp( e, px, py, pz );
    setpart_num( part_n );

    setLifetime();
    setAttribute( "FSR", 0 );
}

EvtVector4C EvtPhotonParticle::epsParentPhoton( int i ) const
{
    EvtVector4C eps;

    switch ( i ) {
        case 0:
            eps.set( EvtComplex( 0.0, 0.0 ),
                     EvtComplex( -1.0 / sqrt( 2.0 ), 0.0 ),
                     EvtComplex( 0.0, -1.0 / sqrt( 2.0 ) ),
                     EvtComplex( 0.0, 0.0 ) );
            break;
        case 1:
            eps.set( EvtComplex( 0.0, 0.0 ), EvtComplex( 1.0 / sqrt( 2.0 ), 0.0 ),
                     EvtComplex( 0.0, -1.0 / sqrt( 2.0 ) ),
                     EvtComplex( 0.0, 0.0 ) );
            break;
        default:
            EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                << "EvtPhotonParticle.cc: Asked "
                << "for state:" << i << endl;
            ::abort();
            break;
    }
    // These are for photon along z axis.  Rotate to get
    // correct direction...

    const double px = this->getP4().get( 1 );
    const double py = this->getP4().get( 2 );
    const double pz = this->getP4().get( 3 );

    const double phi = atan2( py, px );
    const double theta = acos( pz / sqrt( px * px + py * py + pz * pz ) );
    eps.applyRotateEuler( phi, theta, -phi );

    return eps;
}

EvtVector4C EvtPhotonParticle::epsPhoton( int ) const
{
    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "EvtPhotonParticle.cc: Can not get "
        << "state in photons restframe." << endl;
    ;
    ::abort();
    return EvtVector4C();
}

EvtSpinDensity EvtPhotonParticle::rotateToHelicityBasis() const
{
    const EvtVector4C eplus( 0.0, -1.0 / sqrt( 2.0 ),
                             EvtComplex( 0.0, -1.0 / sqrt( 2.0 ) ), 0.0 );
    const EvtVector4C eminus( 0.0, 1.0 / sqrt( 2.0 ),
                              EvtComplex( 0.0, -1.0 / sqrt( 2.0 ) ), 0.0 );

    const EvtVector4C e1 = this->epsParentPhoton( 0 );
    const EvtVector4C e2 = this->epsParentPhoton( 1 );

    EvtSpinDensity R;
    R.setDim( 2 );

    R.set( 0, 0, ( eplus.conj() ) * e1 );
    R.set( 0, 1, ( eplus.conj() ) * e2 );

    R.set( 1, 0, ( eminus.conj() ) * e1 );
    R.set( 1, 1, ( eminus.conj() ) * e2 );

    return R;
}

EvtSpinDensity EvtPhotonParticle::rotateToHelicityBasis( double alpha,
                                                         double beta,
                                                         double gamma ) const
{
    EvtVector4C eplus( 0.0, -1.0 / sqrt( 2.0 ),
                       EvtComplex( 0.0, -1.0 / sqrt( 2.0 ) ), 0.0 );
    EvtVector4C eminus( 0.0, 1.0 / sqrt( 2.0 ),
                        EvtComplex( 0.0, -1.0 / sqrt( 2.0 ) ), 0.0 );

    eplus.applyRotateEuler( alpha, beta, gamma );
    eminus.applyRotateEuler( alpha, beta, gamma );

    const EvtVector4C e1 = this->epsParentPhoton( 0 );
    const EvtVector4C e2 = this->epsParentPhoton( 1 );

    EvtSpinDensity R;
    R.setDim( 2 );

    R.set( 0, 0, ( eplus.conj() ) * e1 );
    R.set( 0, 1, ( eplus.conj() ) * e2 );

    R.set( 1, 0, ( eminus.conj() ) * e1 );
    R.set( 1, 1, ( eminus.conj() ) * e2 );

    return R;
}
