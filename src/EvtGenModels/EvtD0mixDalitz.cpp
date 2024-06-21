
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

#include "EvtGenModels/EvtD0mixDalitz.hh"

#include "EvtGenBase/EvtDalitzPlot.hh"
#include "EvtGenBase/EvtDalitzReso.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtResonance.hh"

#include <cmath>    // for std::fabs

// Initialize the static variables.
const EvtSpinType::spintype& EvtD0mixDalitz::m_SCALAR = EvtSpinType::SCALAR;
const EvtSpinType::spintype& EvtD0mixDalitz::m_VECTOR = EvtSpinType::VECTOR;
const EvtSpinType::spintype& EvtD0mixDalitz::m_TENSOR = EvtSpinType::TENSOR;

const EvtDalitzReso::CouplingType& EvtD0mixDalitz::m_EtaPic = EvtDalitzReso::EtaPic;
const EvtDalitzReso::CouplingType& EvtD0mixDalitz::m_PicPicKK =
    EvtDalitzReso::PicPicKK;

const EvtDalitzReso::NumType& EvtD0mixDalitz::m_RBW =
    EvtDalitzReso::RBW_CLEO_ZEMACH;
const EvtDalitzReso::NumType& EvtD0mixDalitz::m_GS = EvtDalitzReso::GS_CLEO_ZEMACH;
const EvtDalitzReso::NumType& EvtD0mixDalitz::m_KMAT = EvtDalitzReso::K_MATRIX;

const EvtCyclic3::Pair& EvtD0mixDalitz::m_AB = EvtCyclic3::AB;
const EvtCyclic3::Pair& EvtD0mixDalitz::m_AC = EvtCyclic3::AC;
const EvtCyclic3::Pair& EvtD0mixDalitz::m_BC = EvtCyclic3::BC;

void EvtD0mixDalitz::init()
{
    // check that there are 0 arguments
    checkNDaug( 3 );

    if ( getNArg() ) {
        if ( getNArg() == 2 ) {
            m_x = getArg( 0 );
            m_y = getArg( 1 );
        } else if ( getNArg() == 4 ) {
            m_x = getArg( 0 );
            m_y = getArg( 1 );
            m_qp = EvtComplex( getArg( 2 ), getArg( 3 ) );
        } else if ( getNArg() == 5 ) {
            m_x = getArg( 0 );
            m_y = getArg( 1 );
            m_qp = EvtComplex( getArg( 2 ), getArg( 3 ) );
            m_isRBWmodel = !getArg(
                4 );    // RBW by default. If arg4 is set, do K-matrix.
        } else {
            EvtGenReport( EVTGEN_ERROR, "EvtD0mixDalitz" )
                << "Number of arguments for this model must be 0, 2, 4 or 5:"
                << std::endl
                << "[ x y ][ qp.re qp.im ][ doK-matrix ]" << std::endl
                << "Check your dec file." << std::endl;
            exit( 1 );
        }
    }

    checkSpinParent( m_SCALAR );
    checkSpinDaughter( 0, m_SCALAR );
    checkSpinDaughter( 1, m_SCALAR );
    checkSpinDaughter( 2, m_SCALAR );

    readPDGValues();

    // Get the EvtId of the D0 and its (3) daughters.
    EvtId parId = getParentId();

    EvtId dau[3];
    for ( int index = 0; index < 3; index++ )
        dau[index] = getDaug( index );

    if ( parId == m_D0 )    // Look for K0bar h+ h-. The order must be K[0SL] h+ h-
        for ( int index = 0; index < 3; index++ )
            if ( ( dau[index] == m_K0B ) || ( dau[index] == m_KS ) ||
                 ( dau[index] == m_KL ) )
                m_d1 = index;
            else if ( ( dau[index] == m_PIP ) || ( dau[index] == m_KP ) )
                m_d2 = index;
            else if ( ( dau[index] == m_PIM ) || ( dau[index] == m_KM ) )
                m_d3 = index;
            else
                reportInvalidAndExit();
    else if ( parId == m_D0B )    // Look for K0 h+ h-. The order must be K[0SL] h- h+
        for ( int index = 0; index < 3; index++ )
            if ( ( dau[index] == m_K0 ) || ( dau[index] == m_KS ) ||
                 ( dau[index] == m_KL ) )
                m_d1 = index;
            else if ( ( dau[index] == m_PIM ) || ( dau[index] == m_KM ) )
                m_d2 = index;
            else if ( ( dau[index] == m_PIP ) || ( dau[index] == m_KP ) )
                m_d3 = index;
            else
                reportInvalidAndExit();
    else
        reportInvalidAndExit();

    // If the D meson is a D0bar, the expressions should use p/q instead of q/p.
    if ( parId == m_D0B )
        m_qp = 1.0 / m_qp;

    // At this point, if parId is D0bar, the amplitude is the D0bar amplitude, the conjugated amplitude
    //    is the amplitude of the D0 decay, and m_qp means p/q, so it is like changing the meaning of
    //    A <-> Abar, and p <-> q. It is just a trick so after this point the code for D0bar can be the
    //    same as the code for D0.

    // Check if we're dealing with Ks pi pi or with Ks K K.
    m_isKsPiPi = false;
    if ( dau[m_d2] == m_PIP || dau[m_d2] == m_PIM )
        m_isKsPiPi = true;
}

void EvtD0mixDalitz::decay( EvtParticle* part )
{
    // Same structure for all of these decays.
    part->initializePhaseSpace( getNDaug(), getDaugs() );
    EvtVector4R pA = part->getDaug( m_d1 )->getP4();
    EvtVector4R pB = part->getDaug( m_d2 )->getP4();
    EvtVector4R pC = part->getDaug( m_d3 )->getP4();

    // Squared invariant masses.
    double m2AB = ( pA + pB ).mass2();
    double m2AC = ( pA + pC ).mass2();
    double m2BC = ( pB + pC ).mass2();

    // Dalitz amplitudes of the decay of the particle and that of the antiparticle.
    EvtComplex ampDalitz;
    EvtComplex ampAntiDalitz;

    if ( m_isKsPiPi ) {    // For Ks pi pi
        EvtDalitzPoint point( m_mKs, m_mPi, m_mPi, m2AB, m2BC, m2AC );
        EvtDalitzPoint antiPoint( m_mKs, m_mPi, m_mPi, m2AC, m2BC, m2AB );

        ampDalitz = dalitzKsPiPi( point );
        ampAntiDalitz = dalitzKsPiPi( antiPoint );
    } else {    // For Ks K K
        EvtDalitzPoint point( m_mKs, m_mK, m_mK, m2AB, m2BC, m2AC );
        EvtDalitzPoint antiPoint( m_mKs, m_mK, m_mK, m2AC, m2BC, m2AB );

        ampDalitz = dalitzKsKK( point );
        ampAntiDalitz = dalitzKsKK( antiPoint );
    }

    // Assume there's no direct CP violation.
    EvtComplex barAOverA = ampAntiDalitz / ampDalitz;

    // CP violation in the interference. m_qp implements CP violation in the mixing.
    EvtComplex chi = m_qp * barAOverA;

    // Generate a negative exponential life time. p( gt ) = ( 1 - y ) * e^{ - ( 1 - y ) gt }
    double gt = -log( EvtRandom::Flat() ) / ( 1.0 - std::fabs( m_y ) );
    part->setLifetime( gt / m_gamma );

    // Compute time dependent amplitude.
    EvtComplex amp = 0.5 * ampDalitz * exp( -std::fabs( m_y ) * gt / 2.0 ) *
                     ( ( 1.0 + chi ) * h1( gt ) + ( 1.0 - chi ) * h2( gt ) );

    vertex( amp );

    return;
}

void EvtD0mixDalitz::readPDGValues()
{
    // Define the EvtIds.
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

    // Compute the decay rate from the parameter in the evt.pdl file.
    m_ctau = EvtPDL::getctau( EvtPDL::getId( "D0" ) );

    m_gamma = 1.0 / m_ctau;    // ALERT: Gamma is not 1 / tau.
}

EvtComplex EvtD0mixDalitz::dalitzKsPiPi( const EvtDalitzPoint& point )
{
    static const EvtDalitzPlot plot( m_mKs, m_mPi, m_mPi, m_mD0 );

    EvtComplex amp = 0.;

    if ( m_isRBWmodel ) {
        // This corresponds to relativistic Breit-Wigner distributions. Not K-matrix.
        // Defining resonances.
        static EvtDalitzReso KStarm( plot, m_BC, m_AC, m_VECTOR, 0.893606,
                                     0.0463407, m_RBW );
        static EvtDalitzReso KStarp( plot, m_BC, m_AB, m_VECTOR, 0.893606,
                                     0.0463407, m_RBW );
        static EvtDalitzReso rho0( plot, m_AC, m_BC, m_VECTOR, 0.7758, 0.1464,
                                   m_GS );
        static EvtDalitzReso omega( plot, m_AC, m_BC, m_VECTOR, 0.78259,
                                    0.00849, m_RBW );
        static EvtDalitzReso f0_980( plot, m_AC, m_BC, m_SCALAR, 0.975, 0.044,
                                     m_RBW );
        static EvtDalitzReso f0_1370( plot, m_AC, m_BC, m_SCALAR, 1.434, 0.173,
                                      m_RBW );
        static EvtDalitzReso f2_1270( plot, m_AC, m_BC, m_TENSOR, 1.2754,
                                      0.1851, m_RBW );
        static EvtDalitzReso K0Starm_1430( plot, m_BC, m_AC, m_SCALAR, 1.459,
                                           0.175, m_RBW );
        static EvtDalitzReso K0Starp_1430( plot, m_BC, m_AB, m_SCALAR, 1.459,
                                           0.175, m_RBW );
        static EvtDalitzReso K2Starm_1430( plot, m_BC, m_AC, m_TENSOR, 1.4256,
                                           0.0985, m_RBW );
        static EvtDalitzReso K2Starp_1430( plot, m_BC, m_AB, m_TENSOR, 1.4256,
                                           0.0985, m_RBW );
        static EvtDalitzReso sigma( plot, m_AC, m_BC, m_SCALAR, 0.527699,
                                    0.511861, m_RBW );
        static EvtDalitzReso sigma2( plot, m_AC, m_BC, m_SCALAR, 1.03327,
                                     0.0987890, m_RBW );
        static EvtDalitzReso KStarm_1680( plot, m_BC, m_AC, m_VECTOR, 1.677,
                                          0.205, m_RBW );

        // Adding terms to the amplitude with their corresponding amplitude and phase terms.
        amp += EvtComplex( 0.848984, 0.893618 );
        amp += EvtComplex( -1.16356, 1.19933 ) * KStarm.evaluate( point );
        amp += EvtComplex( 0.106051, -0.118513 ) * KStarp.evaluate( point );
        amp += EvtComplex( 1.0, 0.0 ) * rho0.evaluate( point );
        amp += EvtComplex( -0.0249569, 0.0388072 ) * omega.evaluate( point );
        amp += EvtComplex( -0.423586, -0.236099 ) * f0_980.evaluate( point );
        amp += EvtComplex( -2.16486, 3.62385 ) * f0_1370.evaluate( point );
        amp += EvtComplex( 0.217748, -0.133327 ) * f2_1270.evaluate( point );
        amp += EvtComplex( 1.62128, 1.06816 ) * K0Starm_1430.evaluate( point );
        amp += EvtComplex( 0.148802, 0.0897144 ) * K0Starp_1430.evaluate( point );
        amp += EvtComplex( 1.15489, -0.773363 ) * K2Starm_1430.evaluate( point );
        amp += EvtComplex( 0.140865, -0.165378 ) * K2Starp_1430.evaluate( point );
        amp += EvtComplex( -1.55556, -0.931685 ) * sigma.evaluate( point );
        amp += EvtComplex( -0.273791, -0.0535596 ) * sigma2.evaluate( point );
        amp += EvtComplex( -1.69720, 0.128038 ) * KStarm_1680.evaluate( point );
    } else {
        // This corresponds to the complete model (RBW, GS, LASS and K-matrix).
        // Defining resonances.
        static EvtDalitzReso KStarm( plot, m_BC, m_AC, m_VECTOR, 0.893619,
                                     0.0466508, m_RBW );
        static EvtDalitzReso KStarp( plot, m_BC, m_AB, m_VECTOR, 0.893619,
                                     0.0466508, m_RBW );
        static EvtDalitzReso rho0( plot, m_AC, m_BC, m_VECTOR, 0.7758, 0.1464,
                                   m_GS );
        static EvtDalitzReso omega( plot, m_AC, m_BC, m_VECTOR, 0.78259,
                                    0.00849, m_RBW );
        static EvtDalitzReso f2_1270( plot, m_AC, m_BC, m_TENSOR, 1.2754,
                                      0.1851, m_RBW );
        static EvtDalitzReso K0Starm_1430( plot, m_AC, 1.46312, 0.232393, 1.0746,
                                           -1.83214, .803516, 2.32788, 1.0,
                                           -5.31306 );    // LASS
        static EvtDalitzReso K0Starp_1430( plot, m_AB, 1.46312, 0.232393, 1.0746,
                                           -1.83214, .803516, 2.32788, 1.0,
                                           -5.31306 );    // LASS
        static EvtDalitzReso K2Starm_1430( plot, m_BC, m_AC, m_TENSOR, 1.4256,
                                           0.0985, m_RBW );
        static EvtDalitzReso K2Starp_1430( plot, m_BC, m_AB, m_TENSOR, 1.4256,
                                           0.0985, m_RBW );
        static EvtDalitzReso KStarm_1680( plot, m_BC, m_AC, m_VECTOR, 1.677,
                                          0.205, m_RBW );

        // Defining K-matrix.
        static EvtComplex fr12( 1.87981, -0.628378 );
        static EvtComplex fr13( 4.3242, 2.75019 );
        static EvtComplex fr14( 3.22336, 0.271048 );
        static EvtComplex fr15( 0.0, 0.0 );
        static EvtDalitzReso Pole1( plot, m_BC, "Pole1", m_KMAT, fr12, fr13,
                                    fr14, fr15, -0.0694725 );
        static EvtDalitzReso Pole2( plot, m_BC, "Pole2", m_KMAT, fr12, fr13,
                                    fr14, fr15, -0.0694725 );
        static EvtDalitzReso Pole3( plot, m_BC, "Pole3", m_KMAT, fr12, fr13,
                                    fr14, fr15, -0.0694725 );
        static EvtDalitzReso Pole4( plot, m_BC, "Pole4", m_KMAT, fr12, fr13,
                                    fr14, fr15, -0.0694725 );
        static EvtDalitzReso kmatrix( plot, m_BC, "f11prod", m_KMAT, fr12, fr13,
                                      fr14, fr15, -0.0694725 );

        // Adding terms to the amplitude with their corresponding amplitude and phase terms.
        amp += EvtComplex( -1.31394, 1.14072 ) * KStarm.evaluate( point );
        amp += EvtComplex( 0.116239, -0.107287 ) * KStarp.evaluate( point );
        amp += EvtComplex( 1.0, 0.0 ) * rho0.evaluate( point );
        amp += EvtComplex( -0.0313343, 0.0424013 ) * omega.evaluate( point );
        amp += EvtComplex( 0.559412, -0.232336 ) * f2_1270.evaluate( point );
        amp += EvtComplex( 7.35400, -3.67637 ) * K0Starm_1430.evaluate( point );
        amp += EvtComplex( 0.255913, -0.190459 ) * K0Starp_1430.evaluate( point );
        amp += EvtComplex( 1.05397, -0.936297 ) * K2Starm_1430.evaluate( point );
        amp += EvtComplex( -0.00760136, -0.0908624 ) *
               K2Starp_1430.evaluate( point );
        amp += EvtComplex( -1.45336, -0.164494 ) * KStarm_1680.evaluate( point );
        amp += EvtComplex( -1.81830, 9.10680 ) * Pole1.evaluate( point );
        amp += EvtComplex( 10.1751, 3.87961 ) * Pole2.evaluate( point );
        amp += EvtComplex( 23.6569, -4.94551 ) * Pole3.evaluate( point );
        amp += EvtComplex( 0.0725431, -9.16264 ) * Pole4.evaluate( point );
        amp += EvtComplex( -2.19449, -7.62666 ) * kmatrix.evaluate( point );

        amp *= 0.97;    // Multiply by a constant in order to use the same maximum as RBW model.
    }

    return amp;
}

EvtComplex EvtD0mixDalitz::dalitzKsKK( const EvtDalitzPoint& point )
{
    static const EvtDalitzPlot plot( m_mKs, m_mK, m_mK, m_mD0 );

    // Defining resonances.
    static EvtDalitzReso a00_980( plot, m_AC, m_BC, m_SCALAR, 0.999, m_RBW,
                                  0.550173, 0.324, m_EtaPic );
    static EvtDalitzReso phi( plot, m_AC, m_BC, m_VECTOR, 1.01943, 0.00459319,
                              m_RBW );
    static EvtDalitzReso a0p_980( plot, m_AC, m_AB, m_SCALAR, 0.999, m_RBW,
                                  0.550173, 0.324, m_EtaPic );
    static EvtDalitzReso f0_1370( plot, m_AC, m_BC, m_SCALAR, 1.350, 0.265,
                                  m_RBW );
    static EvtDalitzReso a0m_980( plot, m_AB, m_AC, m_SCALAR, 0.999, m_RBW,
                                  0.550173, 0.324, m_EtaPic );
    static EvtDalitzReso f0_980( plot, m_AC, m_BC, m_SCALAR, 0.965, m_RBW,
                                 0.695, 0.165, m_PicPicKK );
    static EvtDalitzReso f2_1270( plot, m_AC, m_BC, m_TENSOR, 1.2754, 0.1851,
                                  m_RBW );
    static EvtDalitzReso a00_1450( plot, m_AC, m_BC, m_SCALAR, 1.474, 0.265,
                                   m_RBW );
    static EvtDalitzReso a0p_1450( plot, m_AC, m_AB, m_SCALAR, 1.474, 0.265,
                                   m_RBW );
    static EvtDalitzReso a0m_1450( plot, m_AB, m_AC, m_SCALAR, 1.474, 0.265,
                                   m_RBW );

    // Adding terms to the amplitude with their corresponding amplitude and phase terms.
    EvtComplex amp( 0., 0. );    // Phase space amplitude.
    amp += EvtComplex( 1.0, 0.0 ) * a00_980.evaluate( point );
    amp += EvtComplex( -0.126314, 0.188701 ) * phi.evaluate( point );
    amp += EvtComplex( -0.561428, 0.0135338 ) * a0p_980.evaluate( point );
    amp += EvtComplex( 0.035, -0.00110488 ) * f0_1370.evaluate( point );
    amp += EvtComplex( -0.0872735, 0.0791190 ) * a0m_980.evaluate( point );
    amp += EvtComplex( 0.0, 0.0 ) * f0_980.evaluate( point );
    amp += EvtComplex( 0.257341, -0.0408343 ) * f2_1270.evaluate( point );
    amp += EvtComplex( -0.0614342, -0.649930 ) * a00_1450.evaluate( point );
    amp += EvtComplex( -0.104629, 0.830120 ) * a0p_1450.evaluate( point );
    amp += EvtComplex( 0.0, 0.0 ) * a0m_1450.evaluate( point );

    return 2.8 *
           amp;    // Multiply by 2.8 in order to reuse the same probmax as Ks pi pi.
}

// < f | H | D^0 (t) > = 1/2 * [ ( 1 + \chi_f ) * A_f * e_1(gt) + ( 1 - \chi_f ) * A_f * e_2(gt) ]
// < f | H | D^0 (t) > = 1/2 * exp( -gamma t / 2 ) * [ ( 1 + \chi_f ) * A_f * h_1(t) + ( 1 - \chi_f ) * A_f * h_2(t) ]
// e{1,2}( gt ) = exp( -gt / 2 ) * h{1,2}( gt ).
EvtComplex EvtD0mixDalitz::h1( const double& gt ) const
{
    return exp( -EvtComplex( m_y, m_x ) * gt / 2. );
}

EvtComplex EvtD0mixDalitz::h2( const double& gt ) const
{
    return exp( EvtComplex( m_y, m_x ) * gt / 2. );
}
