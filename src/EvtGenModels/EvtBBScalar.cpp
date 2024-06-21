
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

#include "EvtGenModels/EvtBBScalar.hh"

#include "EvtGenBase/EvtDiracSpinor.hh"
#include "EvtGenBase/EvtGammaMatrix.hh"
#include "EvtGenBase/EvtSpinType.hh"
#include "EvtGenBase/EvtTensor4C.hh"

#include <cmath>

using namespace std;

const float pi = 3.14159;
const EvtComplex EvtBBScalar::m_I = EvtComplex( 0, 1 );
const EvtComplex EvtBBScalar::m_V_ub =
    EvtComplex( 3.67e-3 * cos( 60 / 180 * pi ), 3.67e-3 * cos( 60 / 180 * pi ) );
const EvtComplex EvtBBScalar::m_V_us_star = EvtComplex( 0.22, 0 );
const EvtComplex EvtBBScalar::m_a1 = EvtComplex( 1.05, 0 );
const EvtComplex EvtBBScalar::m_V_tb = EvtComplex( 0.99915, 0 );
const EvtComplex EvtBBScalar::m_V_ts_star =
    EvtComplex( -0.04029 - 0.000813 * cos( 60 / 180 * pi ),
                -0.000813 * cos( 60 / 180 * pi ) );
const EvtComplex EvtBBScalar::m_a4 = EvtComplex( -387.3e-4, -121e-4 );
const EvtComplex EvtBBScalar::m_a6 = EvtComplex( -555.3e-4, -121e-4 );
const double EvtBBScalar::m_x[] = { 420.96, -10485.50, 100639.97, -433916.61,
                                    613780.15 };
const double EvtBBScalar::m_y[] = { 292.62, -735.73 };
const double EvtBBScalar::m_ms = 0.120;
const double EvtBBScalar::m_mu = 0.029 * 0.120;
const double EvtBBScalar::m_mb = 4.88;

EvtBBScalar::EvtBBScalar() :
    EvtDecayAmp(), m_massRatio( 0 ), m_baryonMassSum( 0 )
{
    FormFactor dummy;
    dummy.m_value = 0.36;
    dummy.m_sigma1 = 0.43;
    dummy.m_sigma2 = 0.0;
    dummy.m_mV = 5.42;
    m_f1Map.insert( make_pair( string( "K" ), dummy ) );
    dummy.m_sigma1 = 0.70;
    dummy.m_sigma2 = 0.27;
    m_f0Map.insert( make_pair( string( "K" ), dummy ) );
    dummy.m_value = 0.29;
    dummy.m_sigma1 = 0.48;
    dummy.m_sigma2 = 0.0;
    dummy.m_mV = 5.32;
    m_f1Map.insert( make_pair( string( "pi" ), dummy ) );
    dummy.m_sigma1 = 0.76;
    dummy.m_sigma2 = 0.28;
    m_f0Map.insert( make_pair( string( "pi" ), dummy ) );
}

std::string EvtBBScalar::getName()
{
    return "B_TO_2BARYON_SCALAR";
}

EvtBBScalar* EvtBBScalar::clone()
{
    return new EvtBBScalar;
}

void EvtBBScalar::setKnownBaryonTypes( const EvtId& baryon )
{
    int baryonId = EvtPDL::getStdHep( baryon );
    if ( EvtPDL::getStdHep( EvtPDL::getId( "Lambda0" ) ) == baryonId or
         EvtPDL::getStdHep( EvtPDL::getId( "anti-Lambda0" ) ) == baryonId ) {
        m_baryonCombination.set( Lambda );
    } else if ( EvtPDL::getStdHep( EvtPDL::getId( "p+" ) ) == baryonId or
                EvtPDL::getStdHep( EvtPDL::getId( "anti-p-" ) ) == baryonId ) {
        m_baryonCombination.set( Proton );
    } else if ( EvtPDL::getStdHep( EvtPDL::getId( "n0" ) ) == baryonId or
                EvtPDL::getStdHep( EvtPDL::getId( "anti-n0" ) ) == baryonId ) {
        m_baryonCombination.set( Neutron );
    } else if ( EvtPDL::getStdHep( EvtPDL::getId( "Sigma0" ) ) == baryonId or
                EvtPDL::getStdHep( EvtPDL::getId( "anti-Sigma0" ) ) == baryonId ) {
        m_baryonCombination.set( Sigma0 );
    } else if ( EvtPDL::getStdHep( EvtPDL::getId( "Sigma-" ) ) == baryonId or
                EvtPDL::getStdHep( EvtPDL::getId( "anti-Sigma+" ) ) == baryonId ) {
        m_baryonCombination.set( Sigma_minus );
    } else if ( EvtPDL::getStdHep( EvtPDL::getId( "Xi0" ) ) == baryonId or
                EvtPDL::getStdHep( EvtPDL::getId( "anti-Xi0" ) ) == baryonId ) {
        m_baryonCombination.set( Xi0 );
    } else if ( EvtPDL::getStdHep( EvtPDL::getId( "Xi-" ) ) == baryonId or
                EvtPDL::getStdHep( EvtPDL::getId( "anti-Xi+" ) ) == baryonId ) {
        m_baryonCombination.set( Xi_minus );
    } else {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "EvtBBScalar::init: Don't know what to do with this type as the first or second baryon\n";
        exit( 2 );
    }
}

double EvtBBScalar::baryonF1F2( double t ) const
{
    // check for known form factors for combination of baryons
    if ( m_baryonCombination.test( Lambda ) and
         m_baryonCombination.test( Proton ) ) {
        return -sqrt( 1.5 ) * G_p( t );
    } else if ( m_baryonCombination.test( Sigma0 ) and
                m_baryonCombination.test( Proton ) ) {
        return -sqrt( 0.5 ) * ( G_p( t ) + 2 * G_n( t ) );
    } else if ( m_baryonCombination.test( Sigma_minus ) and
                m_baryonCombination.test( Neutron ) ) {
        return -G_p( t ) - 2 * G_n( t );
    } else if ( m_baryonCombination.test( Xi0 ) and
                m_baryonCombination.test( Sigma_minus ) ) {
        return G_p( t ) - G_n( t );
    } else if ( m_baryonCombination.test( Xi_minus ) and
                m_baryonCombination.test( Sigma0 ) ) {
        return sqrt( 0.5 ) * ( G_p( t ) - G_n( t ) );
    } else if ( m_baryonCombination.test( Xi_minus ) and
                m_baryonCombination.test( Lambda ) ) {
        return sqrt( 1.5 ) * ( G_p( t ) + G_n( t ) );
    } else {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "EvtBBScalar::baryonF1F2: Don't know what to do with this type as the first or second baryon\n";
        exit( 2 );
    }
}

double EvtBBScalar::formFactorFit( double t, const vector<double>& params ) const
{
    static const double gamma = 2.148;
    static const double Lambda_0 = 0.3;
    double result = 0;
    for ( size_t i = 0; i < params.size(); ++i ) {
        result += params[i] / pow( t, static_cast<int>( i + 1 ) );
    }
    return result * pow( log( t / pow( Lambda_0, 2 ) ), -gamma );
}

double EvtBBScalar::G_p( double t ) const
{
    const vector<double> v_x( m_x, m_x + 5 );
    return formFactorFit( t, v_x );
}

double EvtBBScalar::G_n( double t ) const
{
    const vector<double> v_y( m_y, m_y + 2 );
    return -formFactorFit( t, v_y );
}

double EvtBBScalar::baryon_gA( double t ) const
{
    // check for known form factors for combination of baryons
    if ( m_baryonCombination.test( Lambda ) and
         m_baryonCombination.test( Proton ) ) {
        return -1 / sqrt( 6. ) * ( D_A( t ) + 3 * F_A( t ) );
    } else if ( m_baryonCombination.test( Sigma0 ) and
                m_baryonCombination.test( Proton ) ) {
        return 1 / sqrt( 2. ) * ( D_A( t ) - F_A( t ) );
    } else if ( m_baryonCombination.test( Sigma_minus ) and
                m_baryonCombination.test( Neutron ) ) {
        return D_A( t ) - F_A( t );
    } else if ( m_baryonCombination.test( Xi0 ) and
                m_baryonCombination.test( Sigma_minus ) ) {
        return D_A( t ) + F_A( t );
    } else if ( m_baryonCombination.test( Xi_minus ) and
                m_baryonCombination.test( Sigma0 ) ) {
        return 1 / sqrt( 2. ) * ( D_A( t ) + F_A( t ) );
    } else if ( m_baryonCombination.test( Xi_minus ) and
                m_baryonCombination.test( Lambda ) ) {
        return -1 / sqrt( 6. ) * ( D_A( t ) - 3 * F_A( t ) );
    } else {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "EvtBBScalar::baryon_gA: Don't know what to do with this type as the first or second baryon\n";
        exit( 2 );
    }
}

double EvtBBScalar::baryon_gP( double t ) const
{
    // check for known form factors for combination of baryons
    if ( m_baryonCombination.test( Lambda ) and
         m_baryonCombination.test( Proton ) ) {
        return -1 / sqrt( 6. ) * ( D_P( t ) + 3 * F_P( t ) );
    } else if ( m_baryonCombination.test( Sigma0 ) and
                m_baryonCombination.test( Proton ) ) {
        return 1 / sqrt( 2. ) * ( D_P( t ) - F_P( t ) );
    } else if ( m_baryonCombination.test( Sigma_minus ) and
                m_baryonCombination.test( Neutron ) ) {
        return D_P( t ) - F_P( t );
    } else if ( m_baryonCombination.test( Xi0 ) and
                m_baryonCombination.test( Sigma_minus ) ) {
        return D_P( t ) + F_P( t );
    } else if ( m_baryonCombination.test( Xi_minus ) and
                m_baryonCombination.test( Sigma0 ) ) {
        return 1 / sqrt( 2. ) * ( D_P( t ) + F_P( t ) );
    } else if ( m_baryonCombination.test( Xi_minus ) and
                m_baryonCombination.test( Lambda ) ) {
        return -1 / sqrt( 6. ) * ( D_P( t ) - 3 * F_P( t ) );
    } else {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "EvtBBScalar::baryon_gP: Don't know what to do with this type as the first or second baryon\n";
        exit( 2 );
    }
}

double EvtBBScalar::baryon_fS( double t ) const
{
    // check for known form factors for combination of baryons
    if ( m_baryonCombination.test( Lambda ) and
         m_baryonCombination.test( Proton ) ) {
        return -1 / sqrt( 6. ) * ( D_S( t ) + 3 * F_S( t ) );
    } else if ( m_baryonCombination.test( Sigma0 ) and
                m_baryonCombination.test( Proton ) ) {
        return 1 / sqrt( 2. ) * ( D_S( t ) - F_S( t ) );
    } else if ( m_baryonCombination.test( Sigma_minus ) and
                m_baryonCombination.test( Neutron ) ) {
        return D_S( t ) - F_S( t );
    } else if ( m_baryonCombination.test( Xi0 ) and
                m_baryonCombination.test( Sigma_minus ) ) {
        return D_S( t ) + F_S( t );
    } else if ( m_baryonCombination.test( Xi_minus ) and
                m_baryonCombination.test( Sigma0 ) ) {
        return 1 / sqrt( 2. ) * ( D_S( t ) + F_S( t ) );
    } else if ( m_baryonCombination.test( Xi_minus ) and
                m_baryonCombination.test( Lambda ) ) {
        return -1 / sqrt( 6. ) * ( D_S( t ) - 3 * F_S( t ) );
    } else {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "EvtBBScalar::baryon_fS: Don't know what to do with this type as the first or second baryon\n";
        exit( 2 );
    }
}

double EvtBBScalar::D_A( double t ) const
{
    const double d_tilde[] = { m_x[0] - 1.5 * m_y[0], -478 };
    const vector<double> v_d_tilde( d_tilde, d_tilde + 2 );
    return formFactorFit( t, v_d_tilde );
}

double EvtBBScalar::F_A( double t ) const
{
    const double f_tilde[] = { 2. / 3 * m_x[0] + 0.5 * m_y[0], -478 };
    const vector<double> v_f_tilde( f_tilde, f_tilde + 2 );
    return formFactorFit( t, v_f_tilde );
}

double EvtBBScalar::D_P( double t ) const
{
    const double d_bar[] = { 1.5 * m_y[0] * m_massRatio, /*-952*/ 0 };
    const vector<double> v_d_bar( d_bar, d_bar + 2 );
    return formFactorFit( t, v_d_bar );
}

double EvtBBScalar::F_P( double t ) const
{
    const double f_bar[] = { ( m_x[0] - 0.5 * m_y[0] ) * m_massRatio, /*-952*/ 0 };
    const vector<double> v_f_bar( f_bar, f_bar + 2 );
    return formFactorFit( t, v_f_bar );
}

double EvtBBScalar::D_S( double t ) const
{
    return -1.5 * m_massRatio * G_n( t );
}

double EvtBBScalar::F_S( double t ) const
{
    return ( G_p( t ) + 0.5 * G_n( t ) ) * m_massRatio;
}

double EvtBBScalar::baryon_hA( double t ) const
{
    return ( 1 / m_massRatio * baryon_gP( t ) - baryon_gA( t ) ) *
           pow( m_baryonMassSum, 2 ) / t;
}

void EvtBBScalar::init()
{
    // no arguments, daughter lambda p_bar pi
    // charge conservation is checked by base class
    checkNArg( 0 );
    checkNDaug( 3 );
    checkSpinParent( EvtSpinType::SCALAR );
    checkSpinDaughter( 0, EvtSpinType::DIRAC );
    checkSpinDaughter( 1, EvtSpinType::DIRAC );
    checkSpinDaughter( 2, EvtSpinType::SCALAR );
    EvtId baryon1 = getDaug( 0 );
    EvtId baryon2 = getDaug( 1 );
    EvtId scalar = getDaug( 2 );
    int scalarId = EvtPDL::getStdHep( scalar );

    // Different form factors for the B-pi or B-K transition.
    if ( scalarId == EvtPDL::getStdHep( EvtPDL::getId( "pi+" ) ) or
         scalarId == EvtPDL::getStdHep( EvtPDL::getId( "pi-" ) ) or
         scalarId == EvtPDL::getStdHep( EvtPDL::getId( "pi0" ) ) ) {
        m_scalarType = "pi";
    } else if ( scalarId == EvtPDL::getStdHep( EvtPDL::getId( "K+" ) ) or
                scalarId == EvtPDL::getStdHep( EvtPDL::getId( "K-" ) ) or
                scalarId == EvtPDL::getStdHep( EvtPDL::getId( "K0" ) ) or
                scalarId == EvtPDL::getStdHep( EvtPDL::getId( "anti-K0" ) ) ) {
        m_scalarType = "K";
    } else {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "EvtBBScalar::init: Can only deal with Kaons or pions as the third particle\n"
            << "\tFound: " << scalarId << endl;
        exit( 2 );
    }
    // check for known particles
    setKnownBaryonTypes( baryon1 );
    setKnownBaryonTypes( baryon2 );
    double mass1 = EvtPDL::getMass( baryon1 );
    double mass2 = EvtPDL::getMass( baryon2 );
    // This whole model deals only with baryons that differ in s-u
    if ( mass1 > mass2 )
        m_massRatio = ( mass1 - mass2 ) / ( m_ms - m_mu );
    else
        m_massRatio = ( mass2 - mass1 ) / ( m_ms - m_mu );
    m_baryonMassSum = mass1 + mass2;
}

// initialize phasespace and calculate the amplitude
void EvtBBScalar::decay( EvtParticle* p )
{
    p->initializePhaseSpace( getNDaug(), getDaugs() );
    EvtVector4R B_Momentum = p->getP4Lab();
    EvtDiracParticle* theLambda = dynamic_cast<EvtDiracParticle*>(
        p->getDaug( 0 ) );
    EvtDiracParticle* theAntiP = dynamic_cast<EvtDiracParticle*>(
        p->getDaug( 1 ) );
    EvtScalarParticle* theScalar = dynamic_cast<EvtScalarParticle*>(
        p->getDaug( 2 ) );
    EvtVector4R scalarMomentum = theScalar->getP4Lab();

    // The amplitude consists of three matrix elements. These will be calculated one by one here.

    // loop over all possible spin states
    for ( int i = 0; i < 2; ++i ) {
        EvtDiracSpinor lambdaPol = theLambda->spParent( i );
        for ( int j = 0; j < 2; ++j ) {
            EvtDiracSpinor antiP_Pol = theAntiP->spParent( j );
            EvtVector4C theAmplitudePartA = amp_A( B_Momentum, scalarMomentum );
            EvtComplex amplitude;
            for ( int index = 0; index < 4; ++index ) {
                amplitude += theAmplitudePartA.get( index ) *
                             ( m_const_B * amp_B( theLambda, lambdaPol,
                                                  theAntiP, antiP_Pol, index ) +
                               m_const_C * amp_C( theLambda, lambdaPol,
                                                  theAntiP, antiP_Pol, index ) );
            }
            vertex( i, j, amplitude );
        }
    }
}

void EvtBBScalar::initProbMax()
{
    // setProbMax(1);
    setProbMax( 0.2 );    // found by trial and error
}

// Form factor f1 for B-pi transition
double EvtBBScalar::B_pi_f1( double t ) const
{
    FormFactor f = m_f1Map[m_scalarType];
    double mv2 = f.m_mV * f.m_mV;
    return f.m_value / ( ( 1 - t / mv2 ) * ( 1 - f.m_sigma1 * t / mv2 +
                                             f.m_sigma2 * t * t / mv2 / mv2 ) );
}

// Form factor f0 for B-pi transition
double EvtBBScalar::B_pi_f0( double t ) const
{
    FormFactor f = m_f0Map[m_scalarType];
    double mv2 = f.m_mV * f.m_mV;
    return f.m_value /
           ( 1 - f.m_sigma1 * t / mv2 + f.m_sigma2 * t * t / mv2 / mv2 );
}

// constants of the B and C parts of the amplitude
const EvtComplex EvtBBScalar::m_const_B = m_V_ub * m_V_us_star * m_a1 -
                                          m_V_tb * m_V_ts_star * m_a4;
const EvtComplex EvtBBScalar::m_const_C = 2 * m_a6 * m_V_tb * m_V_ts_star;

// part A of the amplitude, see hep-ph/0204185
const EvtVector4C EvtBBScalar::amp_A( const EvtVector4R& p4B,
                                      const EvtVector4R& p4Scalar )
{
    double mB2 = p4B.mass2();
    double mScalar2 = p4Scalar.mass2();
    double t = ( p4B - p4Scalar ).mass2();
    return ( ( p4B + p4Scalar ) - ( mB2 - mScalar2 ) / t * ( p4B - p4Scalar ) ) *
               B_pi_f1( t ) +
           ( mB2 - mScalar2 ) / t * ( p4B - p4Scalar ) * B_pi_f0( t );
}

// part B of the amplitude, Vector and Axial Vector parts
const EvtComplex EvtBBScalar::amp_B( const EvtDiracParticle* baryon1,
                                     const EvtDiracSpinor& b1Pol,
                                     const EvtDiracParticle* baryon2,
                                     const EvtDiracSpinor& b2Pol, int index )
{
    return amp_B_vectorPart( baryon1, b1Pol, baryon2, b2Pol, index ) -
           amp_B_axialPart( baryon1, b1Pol, baryon2, b2Pol, index );
}

const EvtComplex EvtBBScalar::amp_B_vectorPart( const EvtDiracParticle* baryon1,
                                                const EvtDiracSpinor& b1Pol,
                                                const EvtDiracParticle* baryon2,
                                                const EvtDiracSpinor& b2Pol,
                                                int index )
{
    double t = ( baryon1->getP4Lab() + baryon2->getP4Lab() ).mass2();
    EvtGammaMatrix gamma;
    for ( int i = 0; i < 4; ++i ) {
        gamma += EvtTensor4C::g().get( index, i ) * EvtGammaMatrix::g( i );
    }
    // The F2 contribution that is written out in the paper is neglected here.
    // see hep-ph/0204185
    EvtDiracSpinor A = EvtComplex( baryonF1F2( t ) ) * b2Pol;
    EvtDiracSpinor Adjb1Pol = b1Pol.adjoint();
    EvtDiracSpinor gammaA = gamma * A;
    return Adjb1Pol * gammaA;
    //    return b1Pol.adjoint()*(gamma*(EvtComplex(baryonF1F2(t))*b2Pol));
}

const EvtComplex EvtBBScalar::amp_B_axialPart( const EvtDiracParticle* baryon1,
                                               const EvtDiracSpinor& b1Pol,
                                               const EvtDiracParticle* baryon2,
                                               const EvtDiracSpinor& b2Pol,
                                               int index )
{
    EvtGammaMatrix gamma;
    for ( int i = 0; i < 4; ++i ) {
        gamma += EvtTensor4C::g().get( index, i ) * EvtGammaMatrix::g( i );
    }
    double t = ( baryon1->getP4Lab() + baryon2->getP4Lab() ).mass2();
    double mSum = baryon1->mass() + baryon2->mass();
    EvtVector4C momentum_upper = ( baryon1->getP4Lab() + baryon2->getP4Lab() );
    EvtVector4C momentum;
    for ( int mu = 0; mu < 0; ++mu ) {
        EvtComplex dummy;
        for ( int i = 0; i < 4; ++i ) {
            dummy += EvtTensor4C::g().get( index, i ) * momentum_upper.get( i );
        }
        momentum.set( mu, dummy );
    }
    return b1Pol.adjoint() *
           ( ( ( baryon_gA( t ) * gamma + EvtGammaMatrix::id() * baryon_hA( t ) /
                                              mSum * momentum.get( index ) ) *
               EvtGammaMatrix::g5() ) *
             b2Pol );
}

// part C of the amplitude, Scalar and Pseudoscalar parts
const EvtComplex EvtBBScalar::amp_C( const EvtDiracParticle* baryon1,
                                     const EvtDiracSpinor& b1Pol,
                                     const EvtDiracParticle* baryon2,
                                     const EvtDiracSpinor& b2Pol, int index )
{
    EvtVector4C baryonSumP4_upper = baryon1->getP4Lab() + baryon2->getP4Lab();
    EvtVector4C baryonSumP4;
    for ( int mu = 0; mu < 4; ++mu ) {
        EvtComplex dummy;
        for ( int i = 0; i < 4; ++i ) {
            dummy += EvtTensor4C::g().get( mu, i ) * baryonSumP4_upper.get( i );
        }
        baryonSumP4.set( mu, dummy );
    }
    double t = ( baryon1->getP4Lab() + baryon2->getP4Lab() ).mass2();
    return baryonSumP4.get( index ) / ( m_mb - m_mu ) *
           ( amp_C_scalarPart( b1Pol, b2Pol, t ) +
             amp_C_pseudoscalarPart( b1Pol, b2Pol, t ) );
}

const EvtComplex EvtBBScalar::amp_C_scalarPart( const EvtDiracSpinor& b1Pol,
                                                const EvtDiracSpinor& b2Pol,
                                                double t )
{
    return baryon_fS( t ) * b1Pol.adjoint() * b2Pol;
}

const EvtComplex EvtBBScalar::amp_C_pseudoscalarPart(
    const EvtDiracSpinor& b1Pol, const EvtDiracSpinor& b2Pol, double t )
{
    return baryon_gP( t ) * b1Pol.adjoint() * ( EvtGammaMatrix::g5() * b2Pol );
}
