
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

#include "EvtGenModels/EvtRareLbToLllFFGutsche.hh"

#include "EvtGenBase/EvtIdSet.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtVector4R.hh"

//-----------------------------------------------------------------------------
// Implementation file for class : EvtRareLbToLllFFGutsche
//
// 2014-10-21 : Michal Kreps
//-----------------------------------------------------------------------------

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================

void EvtRareLbToLllFFGutsche::init()
{
    m_Vconsts[0][0] = 0.107;
    m_Vconsts[0][1] = 2.27;
    m_Vconsts[0][2] = 1.367;
    m_Vconsts[1][0] = 0.043;
    m_Vconsts[1][1] = 2.411;
    m_Vconsts[1][2] = 1.531;
    m_Vconsts[2][0] = -0.003;
    m_Vconsts[2][1] = 2.815;
    m_Vconsts[2][2] = 2.041;

    m_Aconsts[0][0] = 0.104;
    m_Aconsts[0][1] = 2.232;
    m_Aconsts[0][2] = 1.328;
    m_Aconsts[1][0] = -0.003;
    m_Aconsts[1][1] = 2.955;
    m_Aconsts[1][2] = 3.620;
    m_Aconsts[2][0] = -0.052;
    m_Aconsts[2][1] = 2.437;
    m_Aconsts[2][2] = 1.559;

    m_TVconsts[0][0] = -0.043;
    m_TVconsts[0][1] = 2.411;
    m_TVconsts[0][2] = 1.531;
    m_TVconsts[1][0] = -0.105;
    m_TVconsts[1][1] = 2.27118;
    m_TVconsts[1][2] = 1.36776;
    m_TVconsts[2][0] = 0;    // Not used anywhere
    m_TVconsts[2][1] = 0;
    m_TVconsts[2][2] = 0;

    m_TAconsts[0][0] = 0.003;
    m_TAconsts[0][1] = 2.955;
    m_TAconsts[0][2] = 3.620;
    m_TAconsts[1][0] = -0.105;
    m_TAconsts[1][1] = 2.233;
    m_TAconsts[1][2] = 1.328;
    m_TAconsts[2][0] = 0;    // Not used anywhere
    m_TAconsts[2][1] = 0;
    m_TAconsts[2][2] = 0;

    EvtGenReport( EVTGEN_INFO, "EvtGen" )
        << " EvtRareLbToLll is using form factors from arXiv:1301.3737 "
        << std::endl;
}

//=============================================================================

void EvtRareLbToLllFFGutsche::getFF( const EvtParticle& parent,
                                     const EvtParticle& lambda,
                                     EvtRareLbToLllFFBase::FormFactors& FF ) const
{
    // Find the FF information for this particle, start by setting all to zero
    FF.areZero();

    /*
  if ( ! ( m_Parents.contains(parent->getId()) &&
           m_Daughters.contains(lambda->getId()) ) )
  {
    EvtGenReport(EVTGEN_ERROR,"EvtGen") << " EvtRareLbToLllFFGutsche: Unknown mother and/or daughter. " << std::endl;
    return;
  }
*/

    const double m1 = parent.getP4().mass();
    const double m2 = lambda.getP4().mass();
    EvtVector4R p4parent;
    p4parent.set( parent.mass(), 0, 0, 0 );
    const double q2 = ( p4parent - lambda.getP4() ).mass2();
    const double m21 = m2 / m1;
    const double shat = q2 / m1 / m1;

    double fV[3];
    double fA[3];
    for ( int i = 0; i <= 2; ++i ) {
        fV[i] = formFactorParametrization( shat, m_Vconsts[i][0],
                                           m_Vconsts[i][1], m_Vconsts[i][2] );
        fA[i] = formFactorParametrization( shat, m_Aconsts[i][0],
                                           m_Aconsts[i][1], m_Aconsts[i][2] );
    }
    double fTV[2];
    double fTA[2];
    for ( int i = 0; i <= 1; ++i ) {
        fTV[i] = formFactorParametrization( shat, m_TVconsts[i][0],
                                            m_TVconsts[i][1], m_TVconsts[i][2] );
        fTA[i] = formFactorParametrization( shat, m_TAconsts[i][0],
                                            m_TAconsts[i][1], m_TAconsts[i][2] );
    }

    // Both v^2==v'^2==1 by definition
    FF.m_F[0] = fV[0] + fV[1] * ( 1 + m21 );
    FF.m_F[1] = fV[2] - fV[1];
    FF.m_F[2] = -m21 * ( fV[1] + fV[2] );

    FF.m_G[0] = fA[0] - fA[1] * ( 1 - m21 );
    FF.m_G[1] = fA[2] - fA[1];
    FF.m_G[2] = -m21 * ( +fA[1] + fA[2] );

    FF.m_FT[0] = fTV[1] * ( m1 + m2 ) + fTV[0] * ( q2 / m1 );
    FF.m_FT[1] = +fTV[0] * ( m2 - m1 ) - fTV[1] * m1;
    FF.m_FT[2] = m2 * ( fTV[0] - fTV[1] ) - fTV[0] * m21 * m2;

    FF.m_GT[0] = -fTA[1] * ( m1 - m2 ) + fTA[0] * ( q2 / m1 );
    FF.m_GT[1] = -fTA[1] * m1 + fTA[0] * ( m1 + m2 );
    FF.m_GT[2] = -fTA[0] * m2 * m21 - m2 * ( fTA[0] + fTA[1] );

    return;
}

double EvtRareLbToLllFFGutsche::formFactorParametrization( const double s,
                                                           const double f0,
                                                           const double a,
                                                           const double b ) const
{
    return f0 / ( 1 - a * s + b * s * s );
}
