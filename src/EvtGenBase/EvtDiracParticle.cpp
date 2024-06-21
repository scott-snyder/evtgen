
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

#include "EvtGenBase/EvtDiracParticle.hh"

#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtGammaMatrix.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtSpinDensity.hh"
#include "EvtGenBase/EvtVector4R.hh"

#include <iostream>
#include <math.h>
#include <stdlib.h>
using std::endl;

void EvtDiracParticle::init( EvtId part_n, const EvtVector4R& p4 )
{
    m_validP4 = true;
    setp( p4 );
    setpart_num( part_n );

    if ( EvtPDL::getStdHep( part_n ) == 0 ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Error in EvtDiracParticle::init, part_n=" << part_n.getId()
            << endl;
        ::abort();
    }

    if ( EvtPDL::getStdHep( part_n ) > 0 ) {
        m_spinorRest[0].set( EvtComplex( sqrt( 2.0 * mass() ), 0.0 ),
                             EvtComplex( 0.0, 0.0 ), EvtComplex( 0.0, 0.0 ),
                             EvtComplex( 0.0, 0.0 ) );
        m_spinorRest[1].set( EvtComplex( 0.0, 0.0 ),
                             EvtComplex( sqrt( 2.0 * mass() ), 0.0 ),
                             EvtComplex( 0.0, 0.0 ), EvtComplex( 0.0, 0.0 ) );

        m_spinorParent[0] = boostTo( m_spinorRest[0], p4 );
        m_spinorParent[1] = boostTo( m_spinorRest[1], p4 );

    } else {
        m_spinorRest[0].set( EvtComplex( 0.0, 0.0 ), EvtComplex( 0.0, 0.0 ),
                             EvtComplex( sqrt( 2.0 * mass() ), 0.0 ),
                             EvtComplex( 0.0, 0.0 ) );
        m_spinorRest[1].set( EvtComplex( 0.0, 0.0 ), EvtComplex( 0.0, 0.0 ),
                             EvtComplex( 0.0, 0.0 ),
                             EvtComplex( sqrt( 2.0 * mass() ), 0.0 ) );

        m_spinorParent[0] = boostTo( m_spinorRest[0], p4 );
        m_spinorParent[1] = boostTo( m_spinorRest[1], p4 );
    }

    setLifetime();
}

void EvtDiracParticle::init( EvtId part_n, const EvtVector4R& p4,
                             const EvtDiracSpinor& prod1,
                             const EvtDiracSpinor& prod2,
                             const EvtDiracSpinor& rest1,
                             const EvtDiracSpinor& rest2 )
{
    m_validP4 = true;
    setp( p4 );
    setpart_num( part_n );

    if ( EvtPDL::getStdHep( part_n ) == 0 ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Error in EvtDiracParticle::init, part_n=" << part_n.getId()
            << std::endl;
        ::abort();
    }
    m_spinorRest[0] = rest1;
    m_spinorRest[1] = rest2;
    m_spinorParent[0] = prod1;
    m_spinorParent[1] = prod2;

    setLifetime();
}

EvtSpinDensity EvtDiracParticle::rotateToHelicityBasis() const
{
    EvtDiracSpinor spplus;
    EvtDiracSpinor spminus;

    double sqmt2 = sqrt( 2. * ( getP4().mass() ) );

    if ( EvtPDL::getStdHep( getId() ) > 0 ) {
        spplus.set( 1.0, 0.0, 0.0, 0.0 );
        spminus.set( 0.0, 1.0, 0.0, 0.0 );
    } else {
        spplus.set( 0.0, 0.0, 1.0, 0.0 );
        spminus.set( 0.0, 0.0, 0.0, 1.0 );
    }

    EvtSpinDensity R;
    R.setDim( 2 );

    for ( int i = 0; i < 2; i++ ) {
        R.set( 0, i, ( spplus * m_spinorRest[i] ) / sqmt2 );
        R.set( 1, i, ( spminus * m_spinorRest[i] ) / sqmt2 );
    }

    return R;
}

EvtSpinDensity EvtDiracParticle::rotateToHelicityBasis( double alpha, double beta,
                                                        double gamma ) const
{
    EvtDiracSpinor spplus;
    EvtDiracSpinor spminus;

    double sqmt2 = sqrt( 2. * ( getP4().mass() ) );

    if ( EvtPDL::getStdHep( getId() ) > 0 ) {
        spplus.set( 1.0, 0.0, 0.0, 0.0 );
        spminus.set( 0.0, 1.0, 0.0, 0.0 );
    } else {
        spplus.set( 0.0, 0.0, 1.0, 0.0 );
        spminus.set( 0.0, 0.0, 0.0, 1.0 );
    }

    spplus.applyRotateEuler( alpha, beta, gamma );
    spminus.applyRotateEuler( alpha, beta, gamma );

    EvtSpinDensity R;
    R.setDim( 2 );

    for ( int i = 0; i < 2; i++ ) {
        R.set( 0, i, ( spplus * m_spinorRest[i] ) / sqmt2 );
        R.set( 1, i, ( spminus * m_spinorRest[i] ) / sqmt2 );
    }

    return R;
}
