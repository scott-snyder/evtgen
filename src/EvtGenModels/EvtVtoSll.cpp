
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

#include "EvtGenModels/EvtVtoSll.hh"

#include "EvtGenBase/EvtDiracSpinor.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtIdSet.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtTensor4C.hh"
#include "EvtGenBase/EvtVector4C.hh"
#include "EvtGenBase/EvtVectorParticle.hh"

#include <iostream>
#include <stdlib.h>
#include <string>

std::string EvtVtoSll::getName()
{
    return "VTOSLL";
}

EvtDecayBase* EvtVtoSll::clone()
{
    return new EvtVtoSll;
}

void EvtVtoSll::init()
{
    // check that there are 0 arguments
    checkNArg( 0 );
    checkNDaug( 3 );

    checkSpinParent( EvtSpinType::VECTOR );

    checkSpinDaughter( 0, EvtSpinType::SCALAR );
    checkSpinDaughter( 1, EvtSpinType::DIRAC );
    checkSpinDaughter( 2, EvtSpinType::DIRAC );

    // Work out whether we have electron mode
    const EvtIdSet leptons{ "e-", "e+" };
    if ( leptons.contains( getDaug( 1 ) ) || leptons.contains( getDaug( 2 ) ) ) {
        m_electronMode = true;
        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << " EvtVtoSll has dielectron final state" << std::endl;
    }
}

void EvtVtoSll::initProbMax()
{
    EvtGenReport( EVTGEN_INFO, "EvtGen" )
        << " EvtVtoSll is finding maximum probability ... " << std::endl;

    double theProbMax = 0;
    double theProbMax_q2 = 0;
    double theProbMax_theta = 0;

    EvtVectorParticle parent{};
    parent.noLifeTime();
    parent.init( getParentId(),
                 EvtVector4R( EvtPDL::getMass( getParentId() ), 0, 0, 0 ) );
    parent.setDiagonalSpinDensity();

    EvtAmp amp;
    EvtId daughters[3] = { getDaug( 0 ), getDaug( 1 ), getDaug( 2 ) };
    amp.init( getParentId(), 3, daughters );
    parent.makeDaughters( 3, daughters );
    EvtParticle* scalar = parent.getDaug( 0 );
    EvtParticle* lep1 = parent.getDaug( 1 );
    EvtParticle* lep2 = parent.getDaug( 2 );
    scalar->noLifeTime();
    lep1->noLifeTime();
    lep2->noLifeTime();

    EvtSpinDensity rho;
    rho.setDiag( parent.getSpinStates() );

    const double M0 = EvtPDL::getMass( getParentId() );
    const double mL = EvtPDL::getMass( getDaug( 0 ) );
    const double m1 = EvtPDL::getMass( getDaug( 1 ) );
    const double m2 = EvtPDL::getMass( getDaug( 2 ) );

    const double m12Sum = m1 + m2;
    const double m12Del = m1 - m2;

    const double q2min = ( m1 + m2 ) * ( m1 + m2 );
    const double q2max = ( M0 - mL ) * ( M0 - mL );
    const double mSqSum = M0 * M0 + mL * mL;

    EvtVector4R p4scalar, p4lep1, p4lep2, boost;

    EvtGenReport( EVTGEN_INFO, "EvtGen" )
        << " EvtVtoSll is probing whole phase space ..." << std::endl;

    double prob = 0;
    const int nsteps = 5000;
    for ( int i = 0; i <= nsteps; i++ ) {
        const double q2 = q2min + i * ( q2max - q2min ) / nsteps;
        const double eScalar = ( mSqSum - q2 ) / ( 2 * M0 );
        const double pstar = i == 0 ? 0
                                    : sqrt( q2 - m12Sum * m12Sum ) *
                                          sqrt( q2 - m12Del * m12Del ) /
                                          ( 2 * sqrt( q2 ) );

        boost.set( M0 - eScalar, 0, 0, +sqrt( eScalar * eScalar - mL * mL ) );
        if ( i != nsteps ) {
            p4scalar.set( eScalar, 0, 0, -sqrt( eScalar * eScalar - mL * mL ) );
        } else {
            p4scalar.set( mL, 0, 0, 0 );
        }

        const double pstarSq = pstar * pstar;
        const double E1star = sqrt( pstarSq + m1 * m1 );
        const double E2star = sqrt( pstarSq + m2 * m2 );

        for ( int j = 0; j <= 45; j++ ) {
            const double theta = j * EvtConst::pi / 45;

            const double pstarT = pstar * sin( theta );
            const double pstarZ = pstar * cos( theta );

            p4lep1.set( E1star, 0, pstarT, pstarZ );
            p4lep2.set( E2star, 0, -pstarT, -pstarZ );

            if ( i != nsteps )    // At maximal q2 we are already in correct frame as scalar child and W/Zvirtual are at rest
            {
                p4lep1 = boostTo( p4lep1, boost );
                p4lep2 = boostTo( p4lep2, boost );
            }
            scalar->init( getDaug( 0 ), p4scalar );
            lep1->init( getDaug( 1 ), p4lep1 );
            lep2->init( getDaug( 2 ), p4lep2 );
            calcAmp( parent, amp );
            prob = rho.normalizedProb( amp.getSpinDensity() );
            // In case of electron mode add pole
            if ( m_electronMode ) {
                prob /= 1.0 + m_poleSize / ( q2 * q2 );
            }

            if ( prob > theProbMax ) {
                theProbMax = prob;
                theProbMax_q2 = q2;
                theProbMax_theta = theta;
            }
        }
    }

    theProbMax *= 1.01;

    setProbMax( theProbMax );
    EvtGenReport( EVTGEN_INFO, "EvtGen" )
        << " EvtVtoSll set up maximum probability to " << theProbMax
        << " found at q2 = " << theProbMax_q2 << " ("
        << nsteps * ( theProbMax_q2 - q2min ) / ( q2max - q2min )
        << " %) and theta = " << theProbMax_theta * 180 / EvtConst::pi
        << std::endl;
}

void EvtVtoSll::decay( EvtParticle* parent )
{
    // Phase space initialization depends on what leptons are
    if ( m_electronMode ) {
        setWeight( parent->initializePhaseSpace( getNDaug(), getDaugs(), false,
                                                 m_poleSize, 1, 2 ) );
    } else {
        parent->initializePhaseSpace( getNDaug(), getDaugs() );
    }

    calcAmp( *parent, m_amp2 );
}

void EvtVtoSll::calcAmp( const EvtParticle& parent, EvtAmp& amp ) const
{
    const EvtParticle& l1 = *parent.getDaug( 1 );
    const EvtParticle& l2 = *parent.getDaug( 2 );

    const EvtVector4C l11 = EvtLeptonVCurrent( l1.spParent( 0 ),
                                               l2.spParent( 0 ) );
    const EvtVector4C l12 = EvtLeptonVCurrent( l1.spParent( 0 ),
                                               l2.spParent( 1 ) );
    const EvtVector4C l21 = EvtLeptonVCurrent( l1.spParent( 1 ),
                                               l2.spParent( 0 ) );
    const EvtVector4C l22 = EvtLeptonVCurrent( l1.spParent( 1 ),
                                               l2.spParent( 1 ) );

    const EvtVector4C eps0 = parent.eps( 0 );
    const EvtVector4C eps1 = parent.eps( 1 );
    const EvtVector4C eps2 = parent.eps( 2 );

    const EvtVector4R P = parent.getP4Restframe();
    const EvtVector4R k = l1.getP4() + l2.getP4();
    const double k2 = k * k;

    const EvtTensor4C T(
        dual( EvtGenFunctions::directProd( P, ( 1.0 / k2 ) * k ) ) );

    double M2 = parent.mass();
    M2 *= M2;
    double m2 = l1.mass();
    m2 *= m2;

    const double norm = 1.0 / sqrt( 2 * M2 + 4 * m2 - 4 * m2 * m2 / M2 );

    amp.vertex( 0, 0, 0, norm * ( eps0 * T.cont2( l11 ) ) );
    amp.vertex( 0, 0, 1, norm * ( eps0 * T.cont2( l12 ) ) );
    amp.vertex( 0, 1, 0, norm * ( eps0 * T.cont2( l21 ) ) );
    amp.vertex( 0, 1, 1, norm * ( eps0 * T.cont2( l22 ) ) );

    amp.vertex( 1, 0, 0, norm * ( eps1 * T.cont2( l11 ) ) );
    amp.vertex( 1, 0, 1, norm * ( eps1 * T.cont2( l12 ) ) );
    amp.vertex( 1, 1, 0, norm * ( eps1 * T.cont2( l21 ) ) );
    amp.vertex( 1, 1, 1, norm * ( eps1 * T.cont2( l22 ) ) );

    amp.vertex( 2, 0, 0, norm * ( eps2 * T.cont2( l11 ) ) );
    amp.vertex( 2, 0, 1, norm * ( eps2 * T.cont2( l12 ) ) );
    amp.vertex( 2, 1, 0, norm * ( eps2 * T.cont2( l21 ) ) );
    amp.vertex( 2, 1, 1, norm * ( eps2 * T.cont2( l22 ) ) );

    return;
}
