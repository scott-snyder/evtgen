
/***********************************************************************
* Copyright 1998-2022 CERN for the benefit of the EvtGen authors       *
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

#include "EvtGenModels/EvtRareLbToLllFFlQCD.hh"

#include "EvtGenBase/EvtConst.hh"
#include "EvtGenBase/EvtIdSet.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtVector4R.hh"

#include "EvtGenModels/EvtWilsonCoefficients.hh"

#include <cmath>

//-----------------------------------------------------------------------------
// Implementation file for class : EvtRareLbToLllFFlQCD
//
// 2016-04-19 : Michal Kreps
// 2014-10-22 : Michal Kreps
//-----------------------------------------------------------------------------

void EvtRareLbToLllFFlQCD::init()
{
    const EvtId LbID = EvtPDL::getId( std::string( "Lambda_b0" ) );
    const EvtId LID = EvtPDL::getId( std::string( "Lambda0" ) );
    const EvtId BID = EvtPDL::getId( std::string( "B+" ) );
    const EvtId KID = EvtPDL::getId( std::string( "K-" ) );
    const double m1 = EvtPDL::getMass( LbID );
    const double m2 = EvtPDL::getMass( LID );
    const double mB = EvtPDL::getMass( BID );
    const double mK = EvtPDL::getMass( KID );
    m_t0 = ( m1 - m2 ) * ( m1 - m2 );
    m_tplus = ( mB + mK ) * ( mB + mK );

    m_fconsts[0][0] = 0.4221;
    m_fconsts[0][1] = -1.1386;
    m_fconsts[0][2] = 5.416;
    m_fconsts[1][0] = 0.5182;
    m_fconsts[1][1] = -1.3495;
    m_fconsts[1][2] = 5.416;
    m_fconsts[2][0] = 0.3725;
    m_fconsts[2][1] = -0.9389;
    m_fconsts[2][2] = 5.711;

    m_gconsts[0][0] = 0.3563;
    m_gconsts[0][1] = -1.0612;
    m_gconsts[0][2] = 5.750;
    m_gconsts[1][0] = 0.3563;
    m_gconsts[1][1] = -1.1357;
    m_gconsts[1][2] = 5.750;
    m_gconsts[2][0] = 0.4028;
    m_gconsts[2][1] = -1.0290;
    m_gconsts[2][2] = 5.367;

    m_hconsts[0][0] = 0.4960;
    m_hconsts[0][1] = -1.1275;
    m_hconsts[0][2] = 5.416;
    m_hconsts[1][0] = 0.3876;
    m_hconsts[1][1] = -0.9623;
    m_hconsts[1][2] = 5.416;
    m_hconsts[2][0] = 0;
    m_hconsts[2][1] = 0;
    m_hconsts[2][2] = 0;

    m_htildaconsts[0][0] = 0.3403;
    m_htildaconsts[0][1] = -0.7697;
    m_htildaconsts[0][2] = 5.750;
    m_htildaconsts[1][0] = 0.3403;
    m_htildaconsts[1][1] = -0.8008;
    m_htildaconsts[1][2] = 5.750;
    m_htildaconsts[2][0] = 0;
    m_htildaconsts[2][1] = 0;
    m_htildaconsts[2][2] = 0;

    EvtGenReport( EVTGEN_INFO, "EvtGen" )
        << " EvtRareLbToLll is using form factors from arXiv:1602.01399 "
        << std::endl;
}

//=============================================================================

void EvtRareLbToLllFFlQCD::getFF( const EvtParticle& parent,
                                  const EvtParticle& lambda,
                                  EvtRareLbToLllFFBase::FormFactors& FF ) const
{
    // Find the FF information for this particle, start by setting all to zero
    FF.areZero();

    const double m1 = parent.getP4().mass();
    const double m2 = lambda.getP4().mass();
    //  double m21=m2/m1;
    EvtVector4R p4parent;
    p4parent.set( parent.mass(), 0, 0, 0 );
    const double q2 = ( p4parent - lambda.getP4() ).mass2();

    const double massSum = m1 + m2;
    const double massDiff = m1 - m2;
    const double massSumSq = massSum * massSum;
    const double massDiffSq = massDiff * massDiff;
    const double q2Sum = q2 - massSumSq;
    const double q2Diff = q2 - massDiffSq;

    double f[3];
    double g[3];
    double h[2];
    double htilda[2];

    for ( int i = 0; i <= 2; ++i ) {
        f[i] = formFactorParametrization( q2, m_fconsts[i][0], m_fconsts[i][1],
                                          m_fconsts[i][2] );
        g[i] = formFactorParametrization( q2, m_gconsts[i][0], m_gconsts[i][1],
                                          m_gconsts[i][2] );
    }
    for ( int i = 0; i <= 1; ++i ) {
        h[i] = formFactorParametrization( q2, m_hconsts[i][0], m_hconsts[i][1],
                                          m_hconsts[i][2] );
        htilda[i] = formFactorParametrization( q2, m_htildaconsts[i][0],
                                               m_htildaconsts[i][1],
                                               m_htildaconsts[i][2] );
    }

    // Both v^2==v'^2==1 by definition
    FF.m_F[0] = f[1];
    FF.m_F[1] = m1 *
                ( ( f[1] - f[0] ) * massSum +
                  massDiff *
                      ( q2 * ( f[2] - f[1] ) - ( f[2] - f[0] ) * massSumSq ) /
                      q2 ) /
                q2Sum;
    FF.m_F[2] = -m2 *
                ( massSum * ( f[0] - f[1] ) +
                  massDiff *
                      ( q2 * ( f[2] - f[1] ) - massSumSq * ( f[2] - f[0] ) ) /
                      q2 ) /
                q2Sum;

    FF.m_G[0] = g[1];
    FF.m_G[1] = m1 / q2Diff *
                ( massDiff * ( g[0] - g[1] ) +
                  massSum *
                      ( q2 * ( g[1] - g[2] ) + massDiffSq * ( g[2] - g[0] ) ) /
                      q2 );
    FF.m_G[2] = -m2 / q2Diff *
                ( massDiff * ( g[1] - g[0] ) +
                  massSum *
                      ( q2 * ( g[1] - g[2] ) + massDiffSq * ( g[2] - g[0] ) ) /
                      q2 );

    FF.m_FT[0] = -massSum * h[1];

    FF.m_FT[1] = -m1 / q2Sum *
                 ( 2 * h[1] * m2 * massSum - h[0] * ( q2 - massSum * massDiff ) );
    FF.m_FT[2] = -m2 / q2Sum *
                 ( 2 * h[1] * m1 * massSum - h[0] * ( q2 + massSum * massDiff ) );

    FF.m_GT[0] = massDiff * htilda[1];

    FF.m_GT[1] = m1 / q2Diff *
                 ( 2 * htilda[1] * massDiff * m2 +
                   htilda[0] * ( q2 - massSum * massDiff ) );
    FF.m_GT[2] = m2 / q2Diff *
                 ( -2 * htilda[1] * massDiff * m1 +
                   htilda[0] * ( q2 + massSum * massDiff ) );

    return;
}

double EvtRareLbToLllFFlQCD::formFactorParametrization( const double q2,
                                                        const double a0,
                                                        const double a1,
                                                        const double pole ) const
{
    const double z = zvar( q2 );
    return 1. / ( 1. - q2 / ( pole * pole ) ) * ( a0 + a1 * z );
}

double EvtRareLbToLllFFlQCD::zvar( const double q2 ) const
{
    const double a = std::sqrt( m_tplus - q2 );
    const double b = std::sqrt( m_tplus - m_t0 );

    return ( a - b ) / ( a + b );
}
