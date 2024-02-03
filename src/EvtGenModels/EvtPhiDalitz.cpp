
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

#include "EvtGenModels/EvtPhiDalitz.hh"

#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtVector3C.hh"
#include "EvtGenBase/EvtVector3R.hh"
#include "EvtGenBase/EvtVector4C.hh"
#include "EvtGenBase/EvtVector4R.hh"

#include <math.h>
#include <stdlib.h>
#include <string>

// Implementation of KLOE measurement
// PL B561: 55-60 (2003) + Erratum B609:449-450 (2005)
// or hep-ex/0303016v2

std::string EvtPhiDalitz::getName()
{
    return "PHI_DALITZ";
}

EvtDecayBase* EvtPhiDalitz::clone()
{
    return new EvtPhiDalitz;
}

void EvtPhiDalitz::init()
{
    // check that there are 0 arguments
    checkNArg( 0 );
    checkNDaug( 3 );

    checkSpinParent( EvtSpinType::VECTOR );

    checkSpinDaughter( 0, EvtSpinType::SCALAR );
    checkSpinDaughter( 1, EvtSpinType::SCALAR );
    checkSpinDaughter( 2, EvtSpinType::SCALAR );

    // results taken from KLOE results: arxiv.org/abs/hep-ex/0303016
    m_mRho = 0.7758;
    m_gRho = 0.1439;
    m_aD = 0.78;
    m_phiD = -2.47;
    m_aOmega = 0.0071;
    m_phiOmega = -0.22;

    m_locPip = -1;
    m_locPim = -1;
    m_locPi0 = -1;

    for ( int i = 0; i < 3; i++ ) {
        if ( getDaug( i ) == EvtPDL::getId( "pi+" ) )
            m_locPip = i;
        if ( getDaug( i ) == EvtPDL::getId( "pi-" ) )
            m_locPim = i;
        if ( getDaug( i ) == EvtPDL::getId( "pi0" ) )
            m_locPi0 = i;
    }
    if ( m_locPip == -1 || m_locPim == -1 || m_locPi0 == -1 ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << getModelName()
            << "generator expects daughters to be pi+ pi- pi0\n";
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Found " << EvtPDL::name( getDaug( 0 ) ) << " "
            << EvtPDL::name( getDaug( 1 ) ) << " "
            << EvtPDL::name( getDaug( 2 ) ) << std::endl;
    }
}

void EvtPhiDalitz::initProbMax()
{
    setProbMax( 300.0 );
}

void EvtPhiDalitz::decay( EvtParticle* p )
{
    const EvtId PIP = EvtPDL::getId( "pi+" );
    const EvtId PIM = EvtPDL::getId( "pi-" );
    const EvtId PIZ = EvtPDL::getId( "pi0" );
    const EvtId OMEGA = EvtPDL::getId( "omega" );

    p->initializePhaseSpace( getNDaug(), getDaugs() );

    const EvtVector4R Ppip = p->getDaug( m_locPip )->getP4();
    const EvtVector4R Ppim = p->getDaug( m_locPim )->getP4();
    const EvtVector4R Ppi0 = p->getDaug( m_locPi0 )->getP4();
    const EvtVector4R Qm = ( Ppim + Ppi0 );
    const EvtVector4R Qp = ( Ppip + Ppi0 );
    const EvtVector4R Q0 = ( Ppip + Ppim );
    const double m_pip = EvtPDL::getMeanMass( PIP );
    const double m_pim = EvtPDL::getMeanMass( PIM );
    const double m_pi0 = EvtPDL::getMeanMass( PIZ );
    const double Mrhop = m_mRho;
    const double Mrhom = m_mRho;
    const double Mrho0 = m_mRho;
    const double M2rhop = pow( Mrhop, 2 );
    const double M2rhom = pow( Mrhom, 2 );
    const double M2rho0 = pow( Mrho0, 2 );
    const double M2omega = pow( EvtPDL::getMeanMass( OMEGA ), 2 );

    const double Wrhop = m_gRho;
    const double Wrhom = m_gRho;
    const double Wrho0 = m_gRho;
    const double Womega = EvtPDL::getWidth( OMEGA );

    const double QmM = Qm.mass();
    const double QmM2 = Qm.mass2();
    const double QpM = Qp.mass();
    const double QpM2 = Qp.mass2();
    const double Q0M = Q0.mass();
    const double Q0M2 = Q0.mass2();

    //Rho- Resonance Amplitude
    const double qm = calc_q( QmM, m_pim, m_pi0 );
    const double qm_0 = calc_q( Mrhom, m_pim, m_pi0 );
    const double Gm = Wrhom * pow( qm / qm_0, 3 ) * ( M2rhom / QmM2 );
    const EvtComplex Drhom( ( QmM2 - M2rhom ), QmM * Gm );
    const EvtComplex A1( M2rhom / Drhom );

    //Rho+ Resonance Amplitude
    const double qp = calc_q( QpM, m_pip, m_pi0 );
    const double qp_0 = calc_q( Mrhop, m_pip, m_pi0 );
    const double Gp = Wrhop * pow( qp / qp_0, 3 ) * ( M2rhop / QpM2 );
    const EvtComplex Drhop( ( QpM2 - M2rhop ), QpM * Gp );
    const EvtComplex A2( M2rhop / Drhop );

    //Rho0 Resonance Amplitude
    const double q0 = calc_q( Q0M, m_pip, m_pim );
    const double q0_0 = calc_q( Mrho0, m_pip, m_pim );
    const double G0 = Wrho0 * pow( q0 / q0_0, 3 ) * ( M2rho0 / Q0M2 );
    const EvtComplex Drho0( ( Q0M2 - M2rho0 ), Q0M * G0 );
    const EvtComplex A3( M2rho0 / Drho0 );

    //Omega Resonance Amplitude
    const EvtComplex OmegaA( m_aOmega * cos( m_phiOmega ),
                             m_aOmega * sin( m_phiOmega ) );
    const EvtComplex DOmega( ( Q0M2 - M2omega ), Q0M * Womega );
    const EvtComplex A4( OmegaA * M2omega / DOmega );

    //Direct Decay Amplitude
    const EvtComplex A5( m_aD * cos( m_phiD ), m_aD * sin( m_phiD ) );

    const EvtComplex Atot = A1 + A2 + A3 + A4 + A5;

    // Polarization

    const EvtVector4C ep0 = p->eps( 0 );
    const EvtVector4C ep1 = p->eps( 1 );
    const EvtVector4C ep2 = p->eps( 2 );

    const EvtVector3R p1( Ppip.get( 1 ), Ppip.get( 2 ), Ppip.get( 3 ) );
    const EvtVector3R p2( Ppim.get( 1 ), Ppim.get( 2 ), Ppim.get( 3 ) );
    const EvtVector3R q = cross( p1, p2 );

    const EvtVector3C e1( ep0.get( 1 ), ep0.get( 2 ), ep0.get( 3 ) );
    const EvtVector3C e2( ep1.get( 1 ), ep1.get( 2 ), ep1.get( 3 ) );
    const EvtVector3C e3( ep2.get( 1 ), ep2.get( 2 ), ep2.get( 3 ) );

    //This is an approximate formula of the maximum value that
    //|q| can have.
    const double pM = p->mass();
    const double mSum = Ppip.mass() + Ppim.mass() + Ppi0.mass();
    const double norm = 10.26 / ( pM * pM - mSum * mSum );

    vertex( 0, norm * e1 * q * Atot );
    vertex( 1, norm * e2 * q * Atot );
    vertex( 2, norm * e3 * q * Atot );

    return;
}

double EvtPhiDalitz::calc_q( double M, double m1, double m2 ) const
{
    const double m12Sum = m1 + m2;

    if ( M > m12Sum ) {
        const double MSq = M * M;
        const double m12Diff = m1 - m2;

        return sqrt( ( MSq - ( m12Sum ) * ( m12Sum ) ) *
                     ( MSq - ( m12Diff ) * ( m12Diff ) ) ) /
               ( 2.0 * M );
    } else {
        return 0.0;
    }
}
