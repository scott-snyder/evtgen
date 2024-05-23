
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

#include "EvtGenModels/EvtBcToNPi.hh"

#include "EvtGenBase/EvtIdSet.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtScalarParticle.hh"
#include "EvtGenBase/EvtSpinType.hh"
#include "EvtGenBase/EvtTensor4C.hh"
#include "EvtGenBase/EvtVector4C.hh"

#include <iostream>
using std::endl;

EvtBcToNPi::EvtBcToNPi( bool printAuthorInfo )
{
    m_nCall = 0;
    m_maxAmp2 = 0;
    if ( printAuthorInfo == true ) {
        this->printAuthorInfo();
    }
}

std::string EvtBcToNPi::getName() const
{
    return "EvtBcToNPi";
}

EvtDecayBase* EvtBcToNPi::clone() const
{
    return new EvtBcToNPi;
}

void EvtBcToNPi::init()
{
    // check spins
    checkSpinParent( EvtSpinType::SCALAR );

    // the others are scalar
    for ( int i = 1; i <= ( getNDaug() - 1 ); i++ ) {
        checkSpinDaughter( i, EvtSpinType::SCALAR );
    };

    m_beta = -0.108;
    m_mRho = 0.775;
    m_gammaRho = 0.149;
    m_mRhopr = 1.364;
    m_gammaRhopr = 0.400;
    m_mA1 = 1.23;
    m_gammaA1 = 0.4;

    // read arguments
    if ( EvtPDL::getSpinType( getDaug( 0 ) ) == EvtSpinType::VECTOR ) {
        checkNArg( 10 );
        int n = 0;
        m_maxProb = getArg( n++ );
        m_FA0_N = getArg( n++ );
        m_FA0_c1 = getArg( n++ );
        m_FA0_c2 = getArg( n++ );
        m_FAp_N = getArg( n++ );
        m_FAp_c1 = getArg( n++ );
        m_FAp_c2 = getArg( n++ );
        m_FV_N = getArg( n++ );
        m_FV_c1 = getArg( n++ );
        m_FV_c2 = getArg( n++ );
        m_FAm_N = 0;
        m_FAm_c1 = 0;
        m_FAm_c2 = 0;
    } else if ( EvtPDL::getSpinType( getDaug( 0 ) ) == EvtSpinType::SCALAR ) {
        checkNArg( 4 );
        int n = 0;
        m_maxProb = getArg( n++ );
        m_Fp_N = getArg( n++ );
        m_Fp_c1 = getArg( n++ );
        m_Fp_c2 = getArg( n++ );
        m_Fm_N = 0;
        m_Fm_c1 = 0;
        m_Fm_c2 = 0;
    } else {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Have not yet implemented this final state in BCPSINPI model"
            << endl;
        EvtGenReport( EVTGEN_ERROR, "EvtGen" ) << "Ndaug=" << getNDaug() << endl;
        for ( int id = 0; id < ( getNDaug() - 1 ); id++ )
            EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                << "Daug " << id << " " << EvtPDL::name( getDaug( id ) ).c_str()
                << endl;
        return;
    };

    if ( getNDaug() < 2 || getNDaug() > 4 ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Have not yet implemented this final state in BCPSINPI model"
            << endl;
        EvtGenReport( EVTGEN_ERROR, "EvtGen" ) << "Ndaug=" << getNDaug() << endl;
        for ( int id = 0; id < ( getNDaug() - 1 ); id++ )
            EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                << "Daug " << id << " " << EvtPDL::name( getDaug( id ) ).c_str()
                << endl;
        return;
    }
}

double EvtBcToNPi::energy1( double M, double m1, double m2 )
{
    return ( M * M + m1 * m1 - m2 * m2 ) / ( 2 * M );
}

double EvtBcToNPi::mom1( double M, double m1, double m2 )
{
    double e1 = energy1( M, m1, m2 );
    return sqrt( e1 * e1 - m1 * m1 );
}

void EvtBcToNPi::initProbMax()
{
    if ( m_maxProb > 0. )
        setProbMax( m_maxProb );
    else {
        EvtId id = getParentId();
        EvtScalarParticle* p = new EvtScalarParticle();
        p->init( id, EvtPDL::getMass( id ), 0., 0., 0. );
        p->setDiagonalSpinDensity();
        // add daughters
        p->makeDaughters( getNDaug(), getDaugs() );

        // fill the momenta
        if ( getNDaug() == 2 ) {
            double M = EvtPDL::getMass( id ),
                   m1 = EvtPDL::getMass( getDaug( 0 ) ),
                   m2 = EvtPDL::getMass( getDaug( 1 ) );
            double p1 = mom1( M, m1, m2 );
            p->getDaug( 0 )->setP4(
                EvtVector4R( energy1( M, m1, m2 ), 0., 0., p1 ) );
            p->getDaug( 1 )->setP4(
                EvtVector4R( energy1( M, m2, m1 ), 0., 0., -p1 ) );
        } else if ( getNDaug() == 3 ) {
            double M = EvtPDL::getMass( id ),
                   m1 = EvtPDL::getMass( getDaug( 0 ) ),
                   m2 = EvtPDL::getMass( getDaug( 1 ) ),
                   m3 = EvtPDL::getMass( getDaug( 2 ) );
            double pRho = mom1( M, m1, m_mRho ), pPi = mom1( m_mRho, m2, m3 );
            p->getDaug( 0 )->setP4(
                EvtVector4R( energy1( M, m1, m_mRho ), 0., 0., pRho ) );
            EvtVector4R p4Rho( energy1( M, m_mRho, m1 ), 0., 0., -pRho );
            EvtVector4R p4_2( energy1( m_mRho, m2, m3 ), 0., 0., pPi );
            p4_2.applyBoostTo( p4Rho );
            EvtVector4R p4_3( energy1( m_mRho, m2, m3 ), 0., 0., -pPi );
            p4_3.applyBoostTo( p4Rho );
            p->getDaug( 1 )->setP4( p4_2 );
            p->getDaug( 2 )->setP4( p4_3 );

        } else if ( getNDaug() == 4 ) {
            double M = EvtPDL::getMass( id ),
                   m1 = EvtPDL::getMass( getDaug( 0 ) ),
                   m2 = EvtPDL::getMass( getDaug( 1 ) ),
                   m3 = EvtPDL::getMass( getDaug( 2 ) ),
                   m4 = EvtPDL::getMass( getDaug( 3 ) );
            if ( M < m1 + m_mA1 )
                return;
            double pA1 = mom1( M, m1, m_mA1 ), pRho = mom1( m_mA1, m_mRho, m4 ),
                   pPi = mom1( m_mRho, m2, m3 );
            p->getDaug( 0 )->setP4(
                EvtVector4R( energy1( M, m1, m_mRho ), 0., 0., pA1 ) );
            EvtVector4R p4A1( energy1( M, m_mA1, m1 ), 0., 0., -pA1 );
            EvtVector4R p4Rho( energy1( m_mA1, m_mRho, m4 ), 0, 0, pRho );
            p4Rho.applyBoostTo( p4A1 );
            EvtVector4R p4_4( energy1( m_mA1, m4, m_mRho ), 0, 0, -pRho );
            p4_4.applyBoostTo( p4A1 );
            p->getDaug( 3 )->setP4( p4_4 );
            EvtVector4R p4_2( energy1( m_mRho, m2, m3 ), 0, 0, pPi );
            p4_2.applyBoostTo( p4Rho );
            p->getDaug( 1 )->setP4( p4_2 );
            EvtVector4R p4_3( energy1( m_mRho, m2, m3 ), 0, 0, -pPi );
            p4_2.applyBoostTo( p4Rho );
            p->getDaug( 2 )->setP4( p4_3 );
        };

        m_amp2.init( p->getId(), getNDaug(), getDaugs() );

        decay( p );

        EvtSpinDensity rho = m_amp2.getSpinDensity();

        double prob = p->getSpinDensityForward().normalizedProb( rho );

        if ( prob > 0 )
            setProbMax( 0.9 * prob );
    };
}

void EvtBcToNPi::decay( EvtParticle* root_particle )
{
    ++m_nCall;

    EvtIdSet thePis{ "pi+", "pi-", "pi0" };
    EvtComplex I = EvtComplex( 0.0, 1.0 );

    root_particle->initializePhaseSpace( getNDaug(), getDaugs() );

    EvtVector4R p( root_particle->mass(), 0., 0., 0. ),    // Bc momentum
        k = root_particle->getDaug( 0 )->getP4(),          // J/psi momenta
        Q = p - k;

    double Q2 = Q.mass2();

    // check pi-mesons and calculate hadronic current
    EvtVector4C hardCur;
    bool foundHadCurr = false;

    if ( getNDaug() == 2 )    // Bc -> psi pi+
    {
        hardCur = Q;
        foundHadCurr = true;
    } else if ( getNDaug() == 3 )    // Bc -> psi pi+ pi0
    {
        EvtVector4R p1, p2;
        p1 = root_particle->getDaug( 1 )->getP4(),        // pi+ momenta
            p2 = root_particle->getDaug( 2 )->getP4(),    // pi0 momentum
            hardCur = Fpi( p1, p2 ) * ( p1 - p2 );
        foundHadCurr = true;
    } else if ( getNDaug() == 4 )    // Bc -> psi pi+ pi pi
    {
        int diffPi( 0 ), samePi1( 0 ), samePi2( 0 );
        if ( getDaug( 1 ) == getDaug( 2 ) ) {
            diffPi = 3;
            samePi1 = 1;
            samePi2 = 2;
        }
        if ( getDaug( 1 ) == getDaug( 3 ) ) {
            diffPi = 2;
            samePi1 = 1;
            samePi2 = 3;
        }
        if ( getDaug( 2 ) == getDaug( 3 ) ) {
            diffPi = 1;
            samePi1 = 2;
            samePi2 = 3;
        }

        EvtVector4R p1 = root_particle->getDaug( samePi1 )->getP4();
        EvtVector4R p2 = root_particle->getDaug( samePi2 )->getP4();
        EvtVector4R p3 = root_particle->getDaug( diffPi )->getP4();

        EvtComplex BA1;
        double GA1 = m_gammaA1 * pi3G( Q2, samePi1 ) /
                     pi3G( m_mA1 * m_mA1, samePi1 );
        EvtComplex denBA1( m_mA1 * m_mA1 - Q.mass2(), -1. * m_mA1 * GA1 );
        BA1 = m_mA1 * m_mA1 / denBA1;

        hardCur = BA1 * ( ( p1 - p3 ) -
                          ( Q * ( Q * ( p1 - p3 ) ) / Q2 ) * Fpi( p2, p3 ) +
                          ( p2 - p3 ) -
                          ( Q * ( Q * ( p2 - p3 ) ) / Q2 ) * Fpi( p1, p3 ) );
        foundHadCurr = true;
    }

    if ( !foundHadCurr ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Have not yet implemented this final state in BCNPI model"
            << endl;
        EvtGenReport( EVTGEN_ERROR, "EvtGen" ) << "Ndaug=" << getNDaug() << endl;
        int id;
        for ( id = 0; id < ( getNDaug() - 1 ); id++ )
            EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                << "Daug " << id << " " << EvtPDL::name( getDaug( id ) ).c_str()
                << endl;
        ::abort();
    };

    double amp2 = 0.;
    if ( root_particle->getDaug( 0 )->getSpinType() == EvtSpinType::VECTOR ) {
        EvtTensor4C H;
        double FA0 = m_FA0_N * exp( m_FA0_c1 * Q2 + m_FA0_c2 * Q2 * Q2 );
        double FAp = m_FAp_N * exp( m_FAp_c1 * Q2 + m_FAp_c2 * Q2 * Q2 );
        double FAm = m_FAm_N * exp( m_FAm_c1 * Q2 + m_FAm_c2 * Q2 * Q2 );
        double FV = m_FV_N * exp( m_FV_c1 * Q2 + m_FV_c2 * Q2 * Q2 );
        H = -FA0 * EvtTensor4C::g() -
            FAp * EvtGenFunctions::directProd( p, p + k ) +
            FAm * EvtGenFunctions::directProd( p, p - k ) +
            2 * I * FV * dual( EvtGenFunctions::directProd( p, k ) );
        EvtVector4C Heps = H.cont2( hardCur );

        for ( int i = 0; i < 3; i++ ) {
            EvtVector4C eps = root_particle->getDaug( 0 )
                                  ->epsParent( i )
                                  .conj();    // psi-meson polarization vector
            EvtComplex amp = eps * Heps;
            vertex( i, amp );
            amp2 += pow( abs( amp ), 2 );
        }
    } else if ( root_particle->getDaug( 0 )->getSpinType() ==
                EvtSpinType::SCALAR ) {
        double Fp = m_Fp_N * exp( m_Fp_c1 * Q2 + m_Fp_c2 * Q2 * Q2 );
        double Fm = m_Fm_N * exp( m_Fm_c1 * Q2 + m_Fm_c2 * Q2 * Q2 );
        EvtVector4C H = Fp * ( p + k ) + Fm * ( p - k );
        EvtComplex amp = H * hardCur;
        vertex( amp );
        amp2 += pow( abs( amp ), 2 );
    };
    if ( amp2 > m_maxAmp2 )
        m_maxAmp2 = amp2;

    return;
}

EvtComplex EvtBcToNPi::Fpi( EvtVector4R q1, EvtVector4R q2 )
{
    double m1 = q1.mass();
    double m2 = q2.mass();

    EvtVector4R Q = q1 + q2;
    double mQ2 = Q * Q;

    // momenta in the rho->pipi decay
    double dRho = m_mRho * m_mRho - m1 * m1 - m2 * m2;
    double pPiRho = ( 1.0 / m_mRho ) *
                    sqrt( ( dRho * dRho ) / 4.0 - m1 * m1 * m2 * m2 );

    double dRhopr = m_mRhopr * m_mRhopr - m1 * m1 - m2 * m2;
    double pPiRhopr = ( 1.0 / m_mRhopr ) *
                      sqrt( ( dRhopr * dRhopr ) / 4.0 - m1 * m1 * m2 * m2 );

    double dQ = mQ2 - m1 * m1 - m2 * m2;
    double pPiQ = ( 1.0 / sqrt( mQ2 ) ) *
                  sqrt( ( dQ * dQ ) / 4.0 - m1 * m1 * m2 * m2 );

    double gammaRho = m_gammaRho * m_mRho / sqrt( mQ2 ) *
                      pow( ( pPiQ / pPiRho ), 3 );
    EvtComplex BRhoDem( m_mRho * m_mRho - mQ2, -1.0 * m_mRho * gammaRho );
    EvtComplex BRho = m_mRho * m_mRho / BRhoDem;

    double gammaRhopr = m_gammaRhopr * m_mRhopr / sqrt( mQ2 ) *
                        pow( ( pPiQ / pPiRhopr ), 3 );
    EvtComplex BRhoprDem( m_mRhopr * m_mRhopr - mQ2, -1.0 * m_mRho * gammaRhopr );
    EvtComplex BRhopr = m_mRhopr * m_mRhopr / BRhoprDem;

    return ( BRho + m_beta * BRhopr ) / ( 1 + m_beta );
}

double EvtBcToNPi::pi3G( double m2, int dupD )
{
    double mPi = EvtPDL::getMeanMass( getDaug( dupD ) );
    if ( m2 > ( m_mRho + mPi ) ) {
        return m2 * ( 1.623 + 10.38 / m2 - 9.32 / ( m2 * m2 ) +
                      0.65 / ( m2 * m2 * m2 ) );
    } else {
        double t1 = m2 - 9.0 * mPi * mPi;
        return 4.1 * pow( t1, 3.0 ) * ( 1.0 - 3.3 * t1 + 5.8 * t1 * t1 );
    }
}

void EvtBcToNPi::printAuthorInfo()
{
    EvtGenReport( EVTGEN_INFO, "EvtGen" )
        << "Defining EvtBcToNPi model: Bc -> V + npi and Bc -> P + npi decays\n"
        << "from A.V. Berezhnoy, A.K. Likhoded, A.V. Luchinsky: "
        << "Phys.Rev.D 82, 014012 (2010) and arXiV:1104.0808." << endl;
}
