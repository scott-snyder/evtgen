
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

#include "EvtGenModels/EvtD0gammaDalitz.hh"

#include "EvtGenBase/EvtConst.hh"
#include "EvtGenBase/EvtDecayTable.hh"
#include "EvtGenBase/EvtFlatte.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtResonance.hh"
#include "EvtGenBase/EvtResonance2.hh"

#include <cmath>
#include <cstdlib>
#include <string>

// Initialize the static variables.
const EvtSpinType::spintype& EvtD0gammaDalitz::m_SCALAR = EvtSpinType::SCALAR;
const EvtSpinType::spintype& EvtD0gammaDalitz::m_VECTOR = EvtSpinType::VECTOR;
const EvtSpinType::spintype& EvtD0gammaDalitz::m_TENSOR = EvtSpinType::TENSOR;

const EvtDalitzReso::CouplingType& EvtD0gammaDalitz::m_EtaPic =
    EvtDalitzReso::EtaPic;
const EvtDalitzReso::CouplingType& EvtD0gammaDalitz::m_PicPicKK =
    EvtDalitzReso::PicPicKK;

const EvtDalitzReso::NumType& EvtD0gammaDalitz::m_RBW =
    EvtDalitzReso::RBW_CLEO_ZEMACH;
const EvtDalitzReso::NumType& EvtD0gammaDalitz::m_GS =
    EvtDalitzReso::GS_CLEO_ZEMACH;
const EvtDalitzReso::NumType& EvtD0gammaDalitz::m_KMAT = EvtDalitzReso::K_MATRIX;

const EvtCyclic3::Pair& EvtD0gammaDalitz::m_AB = EvtCyclic3::AB;
const EvtCyclic3::Pair& EvtD0gammaDalitz::m_AC = EvtCyclic3::AC;
const EvtCyclic3::Pair& EvtD0gammaDalitz::m_BC = EvtCyclic3::BC;

std::string EvtD0gammaDalitz::getName() const
{
    return "D0GAMMADALITZ";
}

EvtDecayBase* EvtD0gammaDalitz::clone() const
{
    return new EvtD0gammaDalitz;
}

void EvtD0gammaDalitz::init()
{
    // check that there are 0 arguments
    checkNArg( 0 );

    // Check that this model is valid for the specified decay.
    checkNDaug( 3 );
    checkSpinParent( m_SCALAR );
    checkSpinDaughter( 0, m_SCALAR );
    checkSpinDaughter( 1, m_SCALAR );
    checkSpinDaughter( 2, m_SCALAR );

    // Get the values of the EvtId objects from the data files.
    readPDGValues();

    // Get the EvtId of the D0 and its 3 daughters.
    getParentId();

    EvtId dau[3];
    for ( int index = 0; index < 3; index++ ) {
        dau[index] = getDaug( index );
    }

    // Look for K0bar h+ h-. The order will be K[0SL] h+ h-
    for ( int index = 0; index < 3; index++ ) {
        if ( ( dau[index] == m_K0B ) || ( dau[index] == m_KS ) ||
             ( dau[index] == m_KL ) ) {
            m_d1 = index;
        } else if ( ( dau[index] == m_PIP ) || ( dau[index] == m_KP ) ) {
            m_d2 = index;
        } else if ( ( dau[index] == m_PIM ) || ( dau[index] == m_KM ) ) {
            m_d3 = index;
        } else {
            reportInvalidAndExit();
        }
    }

    // Check if we're dealing with Ks pi pi or with Ks K K.
    m_isKsPiPi = false;
    if ( dau[m_d2] == m_PIP || dau[m_d2] == m_PIM ) {
        m_isKsPiPi = true;
    }
}

void EvtD0gammaDalitz::initProbMax()
{
    setProbMax( 5200. );
}

void EvtD0gammaDalitz::decay( EvtParticle* part )
{
    // Check if the D is from a B+- -> D0 K+- decay with the appropriate model.
    EvtParticle* parent =
        part->getParent();    // If there are no mistakes, should be B+ or B-.
    if ( parent != nullptr &&
         EvtDecayTable::getInstance()->getDecayFunc( parent )->getName() ==
             "BTODDALITZCPK" ) {
        EvtId parId = parent->getId();
        if ( ( parId == m_BP ) || ( parId == m_BM ) || ( parId == m_B0 ) ||
             ( parId == m_B0B ) ) {
            m_bFlavor = parId;
        } else {
            reportInvalidAndExit();
        }
    } else {
        reportInvalidAndExit();
    }

    // Read the D decay parameters from the B decay model.
    // Gamma angle in rad.
    double gamma = EvtDecayTable::getInstance()->getDecayFunc( parent )->getArg(
        0 );
    // Strong phase in rad.
    double delta = EvtDecayTable::getInstance()->getDecayFunc( parent )->getArg(
        1 );
    // Ratio between B->D0K and B->D0barK
    double rB = EvtDecayTable::getInstance()->getDecayFunc( parent )->getArg( 2 );

    // Same structure for all of these decays.
    part->initializePhaseSpace( getNDaug(), getDaugs() );
    EvtVector4R pA = part->getDaug( m_d1 )->getP4();
    EvtVector4R pB = part->getDaug( m_d2 )->getP4();
    EvtVector4R pC = part->getDaug( m_d3 )->getP4();

    // Squared invariant masses.
    double mSqAB = ( pA + pB ).mass2();
    double mSqAC = ( pA + pC ).mass2();
    double mSqBC = ( pB + pC ).mass2();

    EvtComplex amp( 1.0, 0.0 );

    // Direct and conjugated amplitudes.
    EvtComplex ampDir;
    EvtComplex ampCnj;

    if ( m_isKsPiPi ) {
        // Direct and conjugated Dalitz points.
        EvtDalitzPoint pointDir( m_mKs, m_mPi, m_mPi, mSqAB, mSqBC, mSqAC );
        EvtDalitzPoint pointCnj( m_mKs, m_mPi, m_mPi, mSqAC, mSqBC, mSqAB );

        // Direct and conjugated amplitudes.
        ampDir = dalitzKsPiPi( pointDir );
        ampCnj = dalitzKsPiPi( pointCnj );
    } else {
        // Direct and conjugated Dalitz points.
        EvtDalitzPoint pointDir( m_mKs, m_mK, m_mK, mSqAB, mSqBC, mSqAC );
        EvtDalitzPoint pointCnj( m_mKs, m_mK, m_mK, mSqAC, mSqBC, mSqAB );

        // Direct and conjugated amplitudes.
        ampDir = dalitzKsKK( pointDir );
        ampCnj = dalitzKsKK( pointCnj );
    }

    if ( m_bFlavor == m_BP || m_bFlavor == m_B0 ) {
        amp = ampCnj + rB * exp( EvtComplex( 0., delta + gamma ) ) * ampDir;
    } else {
        amp = ampDir + rB * exp( EvtComplex( 0., delta - gamma ) ) * ampCnj;
    }

    vertex( amp );

    return;
}

EvtComplex EvtD0gammaDalitz::dalitzKsPiPi( const EvtDalitzPoint& point ) const
{
    static const EvtDalitzPlot plot( m_mKs, m_mPi, m_mPi, m_mD0 );

    EvtComplex amp = 0.;

    // This corresponds to relativistic Breit-Wigner distributions. Not K-matrix.
    // Defining resonances.
    static const EvtDalitzReso KStarm( plot, m_BC, m_AC, m_VECTOR, 0.893606,
                                       0.0463407, m_RBW );
    static const EvtDalitzReso KStarp( plot, m_BC, m_AB, m_VECTOR, 0.893606,
                                       0.0463407, m_RBW );
    static const EvtDalitzReso rho0( plot, m_AC, m_BC, m_VECTOR, 0.7758, 0.1464,
                                     m_GS );
    static const EvtDalitzReso omega( plot, m_AC, m_BC, m_VECTOR, 0.78259,
                                      0.00849, m_RBW );
    static const EvtDalitzReso f0_980( plot, m_AC, m_BC, m_SCALAR, 0.975, 0.044,
                                       m_RBW );
    static const EvtDalitzReso f0_1370( plot, m_AC, m_BC, m_SCALAR, 1.434,
                                        0.173, m_RBW );
    static const EvtDalitzReso f2_1270( plot, m_AC, m_BC, m_TENSOR, 1.2754,
                                        0.1851, m_RBW );
    static const EvtDalitzReso K0Starm_1430( plot, m_BC, m_AC, m_SCALAR, 1.459,
                                             0.175, m_RBW );
    static const EvtDalitzReso K0Starp_1430( plot, m_BC, m_AB, m_SCALAR, 1.459,
                                             0.175, m_RBW );
    static const EvtDalitzReso K2Starm_1430( plot, m_BC, m_AC, m_TENSOR, 1.4256,
                                             0.0985, m_RBW );
    static const EvtDalitzReso K2Starp_1430( plot, m_BC, m_AB, m_TENSOR, 1.4256,
                                             0.0985, m_RBW );
    static const EvtDalitzReso sigma( plot, m_AC, m_BC, m_SCALAR, 0.527699,
                                      0.511861, m_RBW );
    static const EvtDalitzReso sigma2( plot, m_AC, m_BC, m_SCALAR, 1.03327,
                                       0.0987890, m_RBW );
    static const EvtDalitzReso KStarm_1680( plot, m_BC, m_AC, m_VECTOR, 1.677,
                                            0.205, m_RBW );

    // Adding terms to the amplitude with their corresponding amplitude and phase terms.
    amp += EvtComplex( .848984, .893618 );
    amp += EvtComplex( -1.16356, 1.19933 ) * KStarm.evaluate( point );
    amp += EvtComplex( .106051, -.118513 ) * KStarp.evaluate( point );
    amp += EvtComplex( 1.0, 0.0 ) * rho0.evaluate( point );
    amp += EvtComplex( -.0249569, .0388072 ) * omega.evaluate( point );
    amp += EvtComplex( -.423586, -.236099 ) * f0_980.evaluate( point );
    amp += EvtComplex( -2.16486, 3.62385 ) * f0_1370.evaluate( point );
    amp += EvtComplex( .217748, -.133327 ) * f2_1270.evaluate( point );
    amp += EvtComplex( 1.62128, 1.06816 ) * K0Starm_1430.evaluate( point );
    amp += EvtComplex( .148802, .0897144 ) * K0Starp_1430.evaluate( point );
    amp += EvtComplex( 1.15489, -.773363 ) * K2Starm_1430.evaluate( point );
    amp += EvtComplex( .140865, -.165378 ) * K2Starp_1430.evaluate( point );
    amp += EvtComplex( -1.55556, -.931685 ) * sigma.evaluate( point );
    amp += EvtComplex( -.273791, -.0535596 ) * sigma2.evaluate( point );
    amp += EvtComplex( -1.69720, .128038 ) * KStarm_1680.evaluate( point );

    return amp;
}

EvtComplex EvtD0gammaDalitz::dalitzKsKK( const EvtDalitzPoint& point ) const
{
    static const EvtDalitzPlot plot( m_mKs, m_mK, m_mK, m_mD0 );

    // Defining resonances.
    static const EvtDalitzReso a00_980( plot, m_AC, m_BC, m_SCALAR, 0.999,
                                        m_RBW, .550173, .324, m_EtaPic );
    static const EvtDalitzReso phi( plot, m_AC, m_BC, m_VECTOR, 1.01943,
                                    .00459319, m_RBW );
    static const EvtDalitzReso a0p_980( plot, m_AC, m_AB, m_SCALAR, 0.999,
                                        m_RBW, .550173, .324, m_EtaPic );
    static const EvtDalitzReso f0_1370( plot, m_AC, m_BC, m_SCALAR, 1.350, .265,
                                        m_RBW );
    static const EvtDalitzReso a0m_980( plot, m_AB, m_AC, m_SCALAR, 0.999,
                                        m_RBW, .550173, .324, m_EtaPic );
    static const EvtDalitzReso f0_980( plot, m_AC, m_BC, m_SCALAR, 0.965, m_RBW,
                                       .695, .165, m_PicPicKK );
    static const EvtDalitzReso f2_1270( plot, m_AC, m_BC, m_TENSOR, 1.2754,
                                        .1851, m_RBW );
    static const EvtDalitzReso a00_1450( plot, m_AC, m_BC, m_SCALAR, 1.474,
                                         .265, m_RBW );
    static const EvtDalitzReso a0p_1450( plot, m_AC, m_AB, m_SCALAR, 1.474,
                                         .265, m_RBW );
    static const EvtDalitzReso a0m_1450( plot, m_AB, m_AC, m_SCALAR, 1.474,
                                         .265, m_RBW );

    // Adding terms to the amplitude with their corresponding amplitude and phase terms.
    EvtComplex amp( 0., 0. );    // Phase space amplitude.
    amp += EvtComplex( 1.0, 0.0 ) * a00_980.evaluate( point );
    amp += EvtComplex( -.126314, .188701 ) * phi.evaluate( point );
    amp += EvtComplex( -.561428, .0135338 ) * a0p_980.evaluate( point );
    amp += EvtComplex( .035, -.00110488 ) * f0_1370.evaluate( point );
    amp += EvtComplex( -.0872735, .0791190 ) * a0m_980.evaluate( point );
    amp += EvtComplex( 0., 0. ) * f0_980.evaluate( point );
    amp += EvtComplex( .257341, -.0408343 ) * f2_1270.evaluate( point );
    amp += EvtComplex( -.0614342, -.649930 ) * a00_1450.evaluate( point );
    amp += EvtComplex( -.104629, .830120 ) * a0p_1450.evaluate( point );
    amp += EvtComplex( 0., 0. ) * a0m_1450.evaluate( point );

    return 2.8 *
           amp;    // Multiply by 2.8 in order to reuse the same probmax as Ks pi pi.
}

void EvtD0gammaDalitz::readPDGValues()
{
    // Define the EvtIds.
    m_BP = EvtPDL::getId( "B+" );
    m_BM = EvtPDL::getId( "B-" );
    m_B0 = EvtPDL::getId( "B0" );
    m_B0B = EvtPDL::getId( "anti-B0" );
    m_D0 = EvtPDL::getId( "D0" );
    m_D0B = EvtPDL::getId( "anti-D0" );
    m_KM = EvtPDL::getId( "K-" );
    m_KP = EvtPDL::getId( "K+" );
    m_K0 = EvtPDL::getId( "K0" );
    m_K0B = EvtPDL::getId( "anti-K0" );
    m_KL = EvtPDL::getId( "K_L0" );
    m_KS = EvtPDL::getId( "K_S0" );
    m_PIM = EvtPDL::getId( "pi-" );
    m_PIP = EvtPDL::getId( "pi+" );

    // Read the relevant masses.
    m_mD0 = EvtPDL::getMass( m_D0 );
    m_mKs = EvtPDL::getMass( m_KS );
    m_mPi = EvtPDL::getMass( m_PIP );
    m_mK = EvtPDL::getMass( m_KP );
}

void EvtD0gammaDalitz::reportInvalidAndExit() const
{
    EvtGenReport( EVTGEN_ERROR, "EvtD0gammaDalitz" )
        << "EvtD0gammaDalitz: Invalid mode." << std::endl;
    exit( 1 );
}
