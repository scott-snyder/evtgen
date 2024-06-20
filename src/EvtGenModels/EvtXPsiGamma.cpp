
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

#include "EvtGenModels/EvtXPsiGamma.hh"

#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtTensor4C.hh"
#include "EvtGenBase/EvtTensorParticle.hh"
#include "EvtGenBase/EvtVector4C.hh"

#include <iostream>
#include <stdlib.h>
#include <string>

using namespace std;

std::string EvtXPsiGamma::getName()
{
    return "X38722-+_PSI_GAMMA";
}

EvtDecayBase* EvtXPsiGamma::clone()
{
    return new EvtXPsiGamma;
}

void EvtXPsiGamma::init()
{
    checkNArg( 0, 6 );

    if ( getNArg() == 0 ) {
        // X -> omega psi, rho0 psi couplings from table II in F. Brazzi et al, arXiv:1103.3155
        m_gOmega = 1.58;
        m_gPOmega = -0.74;
        m_gRho = -0.29;
        m_gPRho = 0.28;

        // Decay constants used in F. Brazzi et al, arXiv:1103.3155, taken from J. Sakurai, Currents and Mesons (1969)
        m_fOmega = 0.036;
        m_fRho = 0.121;

    } else {
        m_gOmega = getArg( 0 );
        m_gPOmega = getArg( 1 );
        m_gRho = getArg( 2 );
        m_gPRho = getArg( 3 );
        m_fOmega = getArg( 4 );
        m_fRho = getArg( 5 );
    }

    checkNDaug( 2 );

    checkSpinParent( EvtSpinType::TENSOR );

    checkSpinDaughter( 1, EvtSpinType::VECTOR );

    m_ID0 = getDaug( 0 );

    if ( m_ID0 != EvtPDL::getId( "gamma" ) &&
         m_ID0 != EvtPDL::getId( "omega" ) && m_ID0 != EvtPDL::getId( "rho0" ) ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << " " << getName() << " - Decays with '" << getDaug( 0 ).getName()
            << "' as first daughter are not supported. Choose among: 'gamma', 'omega', or 'rho0'."
            << std::endl;
        ::abort();
    };
}

void EvtXPsiGamma::initProbMax()
{
    double theProbMax = 1.;

    // Create a tensor parent at rest and initialize it
    // Use noLifeTime() cludge to avoid generating random numbers

    EvtTensorParticle parent{};
    parent.noLifeTime();
    parent.init( getParentId(),
                 EvtVector4R( EvtPDL::getMass( getParentId() ), 0, 0, 0 ) );
    parent.setDiagonalSpinDensity();

    // Create daughters and initialize amplitude
    EvtAmp amp;
    EvtId daughters[2] = { getDaug( 0 ), getDaug( 1 ) };
    amp.init( getParentId(), 2, daughters );
    parent.makeDaughters( 2, daughters );

    EvtParticle* child1 = parent.getDaug( 0 );
    EvtParticle* child2 = parent.getDaug( 1 );

    child1->noLifeTime();
    child2->noLifeTime();

    EvtSpinDensity rho;
    rho.setDiag( parent.getSpinStates() );

    // Momentum of daughters in parent's frame
    const double parentMass = EvtPDL::getMass( getParentId() );

    // The daughter CMS momentum pstar (and thus the phase space) is larger if the mass of the daughters is lower.
    // Thus the probability is maximal for the minimal resonance mass for rho0 and omega resonances.
    // For photons the minimal mass is always zero.
    const double d1Mass = EvtPDL::getMinMass( getDaug( 0 ) );

    const double d2Mass = EvtPDL::getMass( getDaug( 1 ) );

    const double pstar = calcPstar( parentMass, d1Mass, d2Mass );

    EvtVector4R p4_1, p4_2;

    const int nsteps = 180;

    double prob_max = 0;
    double theta_max = 0;

    for ( int i = 0; i <= nsteps; i++ ) {
        const double theta = i * EvtConst::pi / nsteps;

        p4_1.set( sqrt( pow( pstar, 2 ) + pow( d1Mass, 2 ) ), 0,
                  +pstar * sin( theta ), +pstar * cos( theta ) );

        p4_2.set( sqrt( pow( pstar, 2 ) + pow( d2Mass, 2 ) ), 0,
                  -pstar * sin( theta ), -pstar * cos( theta ) );

        child1->init( getDaug( 0 ), p4_1 );
        child2->init( getDaug( 1 ), p4_2 );

        calcAmp( parent, amp );

        const double i_prob = rho.normalizedProb( amp.getSpinDensity() );

        if ( i_prob > prob_max ) {
            prob_max = i_prob;
            theta_max = theta;
        }
    }

    EvtGenReport( EVTGEN_INFO, "EvtGen" )
        << " " << getName() << " - probability " << prob_max
        << " found for p* (child momenta in parent's frame) = " << pstar
        << ", at theta* = " << theta_max << std::endl;

    theProbMax *= 1.01 * prob_max;

    // For wide resonances we have to account for the phase space increasing with pstar
    if ( m_ID0 != EvtPDL::getId( "gamma" ) )
        theProbMax /= pstar;

    setProbMax( theProbMax );

    EvtGenReport( EVTGEN_INFO, "EvtGen" )
        << " " << getName() << " - set up maximum probability to " << theProbMax
        << std::endl;
}

void EvtXPsiGamma::decay( EvtParticle* parent )
{
    parent->initializePhaseSpace( getNDaug(), getDaugs() );

    if ( m_ID0 != EvtPDL::getId( "gamma" ) ) {
        // This weight compensates for the phase space becoming small at resonance masses
        // close to kinematic boundary (only for omega and rho which have large widths)
        setWeight( calcPstar( parent->mass(), parent->getDaug( 0 )->mass(),
                              parent->getDaug( 1 )->mass() ) );
    }

    calcAmp( *parent, m_amp2 );
}

void EvtXPsiGamma::calcAmp( EvtParticle& parent, EvtAmp& amp )
{
    static const double mOmega2 =
        pow( EvtPDL::getMeanMass( EvtPDL::getId( "omega" ) ), 2 );
    static const double mRho2 =
        pow( EvtPDL::getMeanMass( EvtPDL::getId( "rho0" ) ), 2 );

    if ( m_ID0 == EvtPDL::getId( "gamma" ) ) {
        for ( int iPsi = 0; iPsi < 3; iPsi++ ) {
            for ( int iGamma = 0; iGamma < 2; iGamma++ ) {
                for ( int iChi = 0; iChi < 4; iChi++ ) {
                    const EvtComplex T2 = fT2(
                        parent.getDaug( 1 )->getP4(),
                        parent.getDaug( 0 )->getP4(), parent.epsTensor( iChi ),
                        parent.getDaug( 1 )->epsParent( iPsi ).conj(),
                        parent.getDaug( 0 )->epsParentPhoton( iGamma ).conj() );
                    const EvtComplex T3 = fT3(
                        parent.getDaug( 1 )->getP4(),
                        parent.getDaug( 0 )->getP4(), parent.epsTensor( iChi ),
                        parent.getDaug( 1 )->epsParent( iPsi ).conj(),
                        parent.getDaug( 0 )->epsParentPhoton( iGamma ).conj() );
                    amp.vertex( iChi, iGamma, iPsi,
                                ( m_fOmega / mOmega2 * m_gOmega +
                                  m_fRho / mRho2 * m_gRho ) *
                                        T2 +
                                    ( m_fOmega / mOmega2 * m_gPOmega +
                                      m_fRho / mRho2 * m_gPRho ) *
                                        T3 );
                }
            }
        }
    } else if ( m_ID0 == EvtPDL::getId( "omega" ) ) {
        for ( int iPsi = 0; iPsi < 3; iPsi++ ) {
            for ( int iVect = 0; iVect < 3; iVect++ ) {
                for ( int iChi = 0; iChi < 4; iChi++ ) {
                    const EvtComplex T2 = fT2(
                        parent.getDaug( 1 )->getP4(),
                        parent.getDaug( 0 )->getP4(), parent.epsTensor( iChi ),
                        parent.getDaug( 1 )->epsParent( iPsi ).conj(),
                        parent.getDaug( 0 )->epsParent( iVect ).conj() );
                    const EvtComplex T3 = fT3(
                        parent.getDaug( 1 )->getP4(),
                        parent.getDaug( 0 )->getP4(), parent.epsTensor( iChi ),
                        parent.getDaug( 1 )->epsParent( iPsi ).conj(),
                        parent.getDaug( 0 )->epsParent( iVect ).conj() );
                    amp.vertex( iChi, iVect, iPsi,
                                m_gOmega * T2 + m_gPOmega * T3 );
                }
            }
        }
        // This is for the m_ID0 == EvtPDL::getId( "rho0" ) case
    } else {
        for ( int iPsi = 0; iPsi < 3; iPsi++ ) {
            for ( int iVect = 0; iVect < 3; iVect++ ) {
                for ( int iChi = 0; iChi < 4; iChi++ ) {
                    const EvtComplex T2 = fT2(
                        parent.getDaug( 1 )->getP4(),
                        parent.getDaug( 0 )->getP4(), parent.epsTensor( iChi ),
                        parent.getDaug( 1 )->epsParent( iPsi ).conj(),
                        parent.getDaug( 0 )->epsParent( iVect ).conj() );
                    const EvtComplex T3 = fT3(
                        parent.getDaug( 1 )->getP4(),
                        parent.getDaug( 0 )->getP4(), parent.epsTensor( iChi ),
                        parent.getDaug( 1 )->epsParent( iPsi ).conj(),
                        parent.getDaug( 0 )->epsParent( iVect ).conj() );
                    amp.vertex( iChi, iVect, iPsi, m_gRho * T2 + m_gPRho * T3 );
                }
            }
        }
    }
}

double EvtXPsiGamma::calcPstar( double parentMass, double d1Mass,
                                double d2Mass ) const
{
    const double pstar =
        sqrt( pow( parentMass, 2 ) - pow( ( d1Mass + d2Mass ), 2 ) ) *
        sqrt( pow( parentMass, 2 ) - pow( ( d2Mass - d1Mass ), 2 ) ) /
        ( 2 * parentMass );

    return pstar;
}

EvtComplex EvtXPsiGamma::fT2( EvtVector4R p, EvtVector4R q, EvtTensor4C epsPI,
                              EvtVector4C epsEps, EvtVector4C epsEta ) const
{
    // T2 term from F. Brazzi et al, arXiv:1103.3155, eq. (10)
    const EvtTensor4C epsPQ = dual(
        EvtGenFunctions::directProd( q, p ) );    // e_{mu nu a b} p^a q^b;

    const EvtVector4C tmp1 = epsPI.cont1( epsEps );
    const EvtVector4C tmp2 = epsPQ.cont1( tmp1 );
    const EvtComplex T2temp =
        tmp2 * epsEta;    // eps^a pi_{a mu} e_{mu nu rho si} p_nu q_rho eta_si

    const EvtVector4C tmp3 = epsPI.cont1( epsEta );
    const EvtVector4C tmp4 = epsPQ.cont1( tmp3 );
    const EvtComplex T2 =
        T2temp +
        tmp4 * epsEps;    // T2 - eta^a pi_{a mu} e_{mu nu rho si} q_nu p_rho eps_si

    return T2;
}

EvtComplex EvtXPsiGamma::fT3( EvtVector4R p, EvtVector4R q, EvtTensor4C epsPI,
                              EvtVector4C epsEps, EvtVector4C epsEta ) const
{
    // T3 term from F. Brazzi et al, arXiv:1103.3155, eq. (11)
    const EvtVector4R Q = p - q;
    const EvtVector4R P = p + q;
    const EvtVector4C tmp1 = epsPI.cont1( Q );    // Q_a pi_{a mu}
    const EvtTensor4C tmp3 = dual( EvtGenFunctions::directProd(
        P, epsEps ) );    // e_{mu nu rho si} P^rho eps^si
    const EvtVector4C tmp4 = tmp3.cont1( tmp1 );
    const EvtComplex T3 =
        tmp4 * epsEta;    // Q_a pi_{a mu} e_{mu nu rho si} P^rho eps_si eta_nu
    return T3;
}
