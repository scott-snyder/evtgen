
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

#include "EvtGenModels/EvtPsi2JpsiPiPi.hh"

#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtTensor4C.hh"
#include "EvtGenBase/EvtVector4C.hh"
#include "EvtGenBase/EvtVector4R.hh"

#include <cmath>

EvtPsi2JpsiPiPi::EvtPsi2JpsiPiPi() :
    m_tree( false ),
    m_phi( 0.0 ),
    m_cosPhi( 1.0 ),
    m_cos2Phi( 1.0 ),
    m_sinPhi( 0.0 ),
    m_sin2Phi( 0.0 )
{
    this->setNLOArrays();
}

void EvtPsi2JpsiPiPi::setNLOArrays()
{
    // Parameters for NLO corrections obtained by fitting distributions
    // shown in Fig 2 of the article
    m_c0[0] = 1.21214;
    m_c0[1] = -2.517;
    m_c0[2] = 4.66947;
    m_c0[3] = 15.0853;
    m_c0[4] = -49.7381;
    m_c0[5] = 35.5604;

    m_c1[0] = -6.74237;
    m_c1[1] = 84.2391;
    m_c1[2] = -389.74;
    m_c1[3] = 823.902;
    m_c1[4] = -808.538;
    m_c1[5] = 299.1;

    m_c2[0] = -1.25073;
    m_c2[1] = 16.2666;
    m_c2[2] = -74.6453;
    m_c2[3] = 156.789;
    m_c2[4] = -154.185;
    m_c2[5] = 57.5711;

    m_s1[0] = -8.01579;
    m_s1[1] = 93.9513;
    m_s1[2] = -451.713;
    m_s1[3] = 1049.67;
    m_s1[4] = -1162.9;
    m_s1[5] = 492.364;

    m_s2[0] = 3.04459;
    m_s2[1] = -26.0901;
    m_s2[2] = 81.1557;
    m_s2[3] = -112.875;
    m_s2[4] = 66.0432;
    m_s2[5] = -10.0446;
}

std::string EvtPsi2JpsiPiPi::getName() const
{
    return "PSI2JPSIPIPI";
}

EvtDecayBase* EvtPsi2JpsiPiPi::clone() const
{
    return new EvtPsi2JpsiPiPi;
}

void EvtPsi2JpsiPiPi::initProbMax()
{
    // Should be OK for all m_phi values
    setProbMax( 1.1 );
}

void EvtPsi2JpsiPiPi::init()
{
    checkNArg( 0, 1 );

    if ( getNArg() == 0 ) {
        m_tree = true;
        m_phi = 0.0;

    } else {
        m_tree = false;
        m_phi = getArg( 0 );    // LO vs NLO mixing angle in radians
    }

    double twoPhi = 2.0 * m_phi;
    m_cosPhi = cos( m_phi );
    m_cos2Phi = cos( twoPhi );
    m_sinPhi = sin( m_phi );
    m_sin2Phi = sin( twoPhi );
}

void EvtPsi2JpsiPiPi::decay( EvtParticle* root )
{
    root->initializePhaseSpace( getNDaug(), getDaugs() );

    EvtVector4R p4 =
        root->getDaug( 0 )->getP4();    // J-psi momentum in psi2 rest frame
    EvtVector4R k1 = root->getDaug( 1 )->getP4();    // pi+ momentum in psi2 rest frame
    double mPiSq = k1.mass2();                       // squared pion mass
    EvtVector4R k2 = root->getDaug( 2 )->getP4();    // pi- momentum in psi2 rest frame
    EvtVector4R tq = k1 - k2;
    EvtVector4R p3 = k1 + k2;
    double p3Sq = p3.mass2();
    double mpipi = p3.mass();
    double corr( 1.0 );

    if ( !m_tree ) {
        // Calculate NLO corrections
        corr = 0.0;
        for ( int iq = 0; iq < m_nQ; ++iq ) {
            corr += ( m_c0[iq] + m_c1[iq] * m_cosPhi + m_c2[iq] * m_cos2Phi +
                      m_s1[iq] * m_sinPhi + m_s2[iq] * m_sin2Phi ) *
                    std::pow( mpipi, iq );
        }
    }

    double mSqTerm = 2.0 * mPiSq / p3Sq;
    EvtTensor4C p3Prod = EvtGenFunctions::directProd( p3, p3 );

    // Eq 14 from the article
    EvtTensor4C L = EvtGenFunctions::directProd( tq, tq ) +
                    ( ( 1.0 - 2.0 * mSqTerm ) / 3.0 ) *
                        ( p3Sq * EvtTensor4C::g() - p3Prod );

    EvtTensor4C T = ( 2.0 / 3.0 ) * ( 1.0 + mSqTerm ) * p3Prod - L;

    for ( int iPsi2 = 0; iPsi2 < 5; ++iPsi2 ) {
        EvtTensor4C epsX = root->epsTensor(
            iPsi2 );    // psi2 polarization tensor in psi2 rest frame
        EvtTensor4C epsXT = cont22( epsX, T );

        for ( int iPsi = 0; iPsi < 3; ++iPsi ) {
            EvtVector4C epsPsi = root->getDaug( 0 )->epsParent(
                iPsi );    // Jpsi polarization vector in psi2 rest frame
            EvtTensor4C epeps = dual( EvtGenFunctions::directProd( epsPsi, p4 ) );
            EvtTensor4C ttt = cont22( epeps, epsXT );

            // Eq 13 from the article
            EvtComplex amp = ttt.trace();

            // NLO corrections
            amp *= corr;

            // Set vertex amplitude component
            vertex( iPsi2, iPsi, amp );
        }
    }
}
