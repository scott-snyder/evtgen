
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

#include "EvtGenModels/EvtLambdaP_BarGamma.hh"

#include "EvtGenBase/EvtDiracParticle.hh"
#include "EvtGenBase/EvtDiracSpinor.hh"
#include "EvtGenBase/EvtGammaMatrix.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtPhotonParticle.hh"
#include "EvtGenBase/EvtSpinType.hh"

#include <stdlib.h>
using std::cout;
using std::endl;

EvtLambdaP_BarGamma::EvtLambdaP_BarGamma() :
    m_mLambdab( 5.624 ),       // Lambda_b mass
    m_mLambda0( 1.115684 ),    // Lambda0 mass
    m_c7Eff( -0.31 ),          // Wilson coefficient
    m_mb( 4.4 ),               // running b mass
    m_mV( 5.42 ),              // pole mass vector current
    m_mA( 5.86 ),              // pole mass axial current
    m_GF( 1.166E-5 ),          // Fermi constant
    m_gLambdab( 16 ),          // coupling constant Lambda_b -> B- p
    m_e0( 1 ),                 // electromagnetic coupling (+1)
    m_g1( 0.64 ),              // heavy-light form factors at q_mSqare
    m_g2( -0.10 ),
    m_f1( 0.64 ),
    m_f2( -0.31 ),
    m_VtbVtsStar( 0.038 )    // |V_tb V_ts^*|
{
}

std::string EvtLambdaP_BarGamma::getName()
{
    return "B_TO_LAMBDA_PBAR_GAMMA";
}

EvtDecayBase* EvtLambdaP_BarGamma::clone()
{
    return new EvtLambdaP_BarGamma;
}

void EvtLambdaP_BarGamma::init()
{
    // no arguments, daughter lambda p_bar gamma
    checkNArg( 0 );
    checkNDaug( 3 );

    checkSpinParent( EvtSpinType::SCALAR );
    checkSpinDaughter( 0, EvtSpinType::DIRAC );
    checkSpinDaughter( 1, EvtSpinType::DIRAC );
    checkSpinDaughter( 2, EvtSpinType::PHOTON );
}

// initialize phasespace and calculate the amplitude
void EvtLambdaP_BarGamma::decay( EvtParticle* p )
{
    EvtComplex I( 0, 1 );

    p->initializePhaseSpace( getNDaug(), getDaugs() );

    EvtDiracParticle* theLambda = static_cast<EvtDiracParticle*>(
        p->getDaug( 0 ) );
    EvtVector4R lambdaMomentum = theLambda->getP4Lab();

    EvtDiracParticle* theAntiP = static_cast<EvtDiracParticle*>( p->getDaug( 1 ) );

    EvtPhotonParticle* thePhoton = static_cast<EvtPhotonParticle*>(
        p->getDaug( 2 ) );
    EvtVector4R photonMomentum =
        thePhoton->getP4Lab();    // get momentum in the same frame

    // loop over all possible spin states
    for ( int i = 0; i < 2; ++i ) {
        EvtDiracSpinor lambdaPol = theLambda->spParent( i );
        for ( int j = 0; j < 2; ++j ) {
            EvtDiracSpinor antiP_Pol = theAntiP->spParent( j );
            for ( int k = 0; k < 2; ++k ) {
                EvtVector4C photonPol = thePhoton->epsParentPhoton(
                    k );    // one of two possible polarization states
                EvtGammaMatrix photonGamma;    // sigma[mu][nu] * epsilon[mu] * k[nu] (watch lower indices)
                for ( int mu = 0; mu < 4; ++mu )
                    for ( int nu = 0; nu < 4; ++nu )
                        photonGamma += EvtGammaMatrix::sigmaLower( mu, nu ) *
                                       photonPol.get( mu ) *
                                       photonMomentum.get( nu );

                EvtComplex amp = -I * m_gLambdab * lambdaPol.adjoint() *
                                 ( ( constA() * EvtGammaMatrix::id() +
                                     constB() * EvtGammaMatrix::g5() ) *
                                   photonGamma *
                                   ( EvtGenFunctions::slash( lambdaMomentum ) +
                                     EvtGenFunctions::slash( photonMomentum ) +
                                     m_mLambdab * EvtGammaMatrix::id() ) /
                                   ( ( lambdaMomentum + photonMomentum ) *
                                         ( lambdaMomentum + photonMomentum ) -
                                     m_mLambdab * m_mLambdab ) *
                                   EvtGammaMatrix::g5() * antiP_Pol );
                // use of parentheses so I do not have to define EvtDiracSpinor*EvtGammaMatrix, which shouldn't be defined to prevent errors in indexing

                vertex( i, j, k, amp );
            }
        }
    }
}

void EvtLambdaP_BarGamma::initProbMax()
{
    // setProbMax(1);
    setProbMax( 9.0000E-13 );    // found by trial and error
}

// form factors at 0
double EvtLambdaP_BarGamma::f0( double fqm, int n ) const
{
    return fqm *
           pow( 1 - pow( m_mLambdab - m_mLambda0, 2 ) / ( m_mV * m_mV ), n );
}

double EvtLambdaP_BarGamma::g0( double gqm, int n ) const
{
    return gqm *
           pow( 1 - pow( m_mLambdab - m_mLambda0, 2 ) / ( m_mA * m_mA ), n );
}

double EvtLambdaP_BarGamma::constA() const
{
    return m_GF / sqrt( 2. ) * m_e0 / ( 8 * EvtConst::pi * EvtConst::pi ) * 2 *
           m_c7Eff * m_mb * m_VtbVtsStar * ( f0( m_f1 ) - f0( m_f2 ) );
}

double EvtLambdaP_BarGamma::constB() const
{
    return m_GF / sqrt( 2. ) * m_e0 / ( 8 * EvtConst::pi * EvtConst::pi ) * 2 *
           m_c7Eff * m_mb * m_VtbVtsStar *
           ( g0( m_g1 ) - ( m_mLambdab - m_mLambda0 ) /
                              ( m_mLambdab + m_mLambda0 ) * g0( m_g2 ) );
}
