/***********************************************************************
* Copyright 1998-2023 CERN for the benefit of the EvtGen authors       *
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

#include "EvtGenModels/EvtBcVPPHad.hh"

#include "EvtGenBase/EvtDiracSpinor.hh"
#include "EvtGenBase/EvtIdSet.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtSpinType.hh"
#include "EvtGenBase/EvtTensor4C.hh"
#include "EvtGenBase/EvtVector4R.hh"

#include <iostream>

std::string EvtBcVPPHad::getName()
{
    return "BC_VPPHAD";
}

EvtDecayBase* EvtBcVPPHad::clone()
{
    return new EvtBcVPPHad;
}

//======================================================

void EvtBcVPPHad::init()
{
    checkNArg( 1 );
    checkSpinParent( EvtSpinType::SCALAR );
    checkSpinDaughter( 0, EvtSpinType::VECTOR );
    checkSpinDaughter( 1, EvtSpinType::DIRAC );
    checkSpinDaughter( 2, EvtSpinType::DIRAC );
    for ( int i = 3; i <= ( getNDaug() - 1 ); i++ ) {
        checkSpinDaughter( i, EvtSpinType::SCALAR );
    }

    m_idVector = getDaug( 0 ).getId();
    m_whichFit = int( getArg( 0 ) + 0.1 );
    m_FFModel = std::make_unique<EvtBCVFF2>( m_idVector, m_whichFit );

    m_WCurr = std::make_unique<EvtWHad>();

    // determine the code of final hadronic state
    const EvtIdSet thePis{ "pi+", "pi-", "pi0" };
    const EvtIdSet theK{ "K+", "K-", "K_S0" };
    const EvtIdSet thePs{ "p+", "anti-p-" };

    if ( getNDaug() == 4 && thePs.contains( getDaug( 1 ) ) &&
         thePs.contains( getDaug( 2 ) ) && thePis.contains( getDaug( 3 ) ) ) {
        m_outCode = 1;    // p+ p- pi+
    } else {
        EvtGenReport( EVTGEN_ERROR, "EvtBcPPHad::init has unknown decay mode" );
    }
    EvtGenReport( EVTGEN_INFO, "EvtBcVPPHad" )
        << ": m_outCode = " << m_outCode << ", m_whichFit = " << m_whichFit
        << std::endl;
}

//======================================================

void EvtBcVPPHad::initProbMax()
{
    if ( m_outCode == 1 ) {
        // p+ p- pi+
        if ( m_idVector == EvtPDL::getId( "J/psi" ).getId() ) {
            if ( m_whichFit == 2 ) {
                setProbMax( 550.0 );
            } else {
                setProbMax( 1100.0 );
            }
        } else if ( m_idVector == EvtPDL::getId( "psi(2S)" ).getId() ) {
            if ( m_whichFit == 2 ) {
                setProbMax( 7.0 );
            } else {
                setProbMax( 50.0 );
            }
        }
    } else {
        EvtGenReport(
            EVTGEN_ERROR,
            "probmax: Have not yet implemented this final state in BC_VPPHAD model, m_outCode = " )
            << m_outCode << std::endl;
    }
}

//======================================================
EvtVector4C EvtBcVPPHad::hardCurrPP( EvtParticle* parent, int i1, int i2 ) const
{
    EvtVector4C hardCur;
    const EvtVector4R p1 = parent->getDaug( 1 )->getP4();
    const EvtDiracSpinor sp1 = parent->getDaug( 1 )->spParent( i1 );

    const EvtVector4R p2 = parent->getDaug( 2 )->getP4();
    const EvtDiracSpinor sp2 = parent->getDaug( 2 )->spParent( i2 );

    if ( m_outCode == 1 ) {
        const EvtVector4R k = parent->getDaug( 3 )->getP4();
        hardCur = m_WCurr->WCurrent_ppPi( p1, sp1, p2, sp2, k );
    }
    return hardCur;
}

//======================================================
void EvtBcVPPHad::decay( EvtParticle* parent )
{
    parent->initializePhaseSpace( getNDaug(), getDaugs() );

    EvtParticle* Jpsi = parent->getDaug( 0 );

    const EvtVector4R p4b( parent->mass(), 0., 0., 0. ),    // Bc momentum
        p4meson = Jpsi->getP4(),                            // J/psi momenta
        Q = p4b - p4meson, p4Sum = p4meson + p4b;
    const double Q2 = Q.mass2();

    // Calculate Bc -> V W form-factors
    double a1f( 0.0 ), a2f( 0.0 ), vf( 0.0 ), a0f( 0.0 );

    const double mMeson = Jpsi->mass();
    const double mB = parent->mass();
    const double mVar = mB + mMeson;

    m_FFModel->getvectorff( parent->getId(), Jpsi->getId(), Q2, mMeson, &a1f,
                            &a2f, &vf, &a0f );

    const double a3f = ( mVar / ( 2.0 * mMeson ) ) * a1f -
                       ( ( mB - mMeson ) / ( 2.0 * mMeson ) ) * a2f;

    // Calculate Bc -> V W current
    EvtTensor4C H = a1f * mVar * EvtTensor4C::g();
    H.addDirProd( ( -a2f / mVar ) * p4b, p4Sum );
    H += EvtComplex( 0.0, vf / mVar ) *
         dual( EvtGenFunctions::directProd( p4Sum, Q ) );
    H.addDirProd( ( a0f - a3f ) * 2.0 * ( mMeson / Q2 ) * p4b, Q );

    for ( int i1 = 0; i1 < 2; ++i1 ) {
        for ( int i2 = 0; i2 < 2; ++i2 ) {
            const EvtVector4C hardCur = hardCurrPP( parent, i1, i2 );
            const EvtVector4C Heps = H.cont2( hardCur );
            for ( int i = 0; i < 3; i++ ) {
                // psi-meson polarization vector
                const EvtVector4C eps = Jpsi->epsParent( i ).conj();
                const EvtComplex amp = eps * Heps;
                vertex( i, i1, i2, amp );
            }
        }
    }
}
