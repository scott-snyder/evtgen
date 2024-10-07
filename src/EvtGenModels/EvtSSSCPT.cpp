
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

#include "EvtGenModels/EvtSSSCPT.hh"

#include "EvtGenBase/EvtCPUtil.hh"
#include "EvtGenBase/EvtConst.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtReport.hh"

#include <stdlib.h>
#include <string>

std::string EvtSSSCPT::getName() const
{
    return "SSS_CPT";
}

EvtDecayBase* EvtSSSCPT::clone() const
{
    return new EvtSSSCPT;
}

void EvtSSSCPT::init()
{
    // check that there are 8 arguments
    checkNArg( 8 );
    checkNDaug( 2 );

    // Set amplitude coeffs
    setAmpCoeffs();
}

void EvtSSSCPT::setAmpCoeffs()
{
    m_P = EvtComplex( cos( -getArg( 0 ) ), sin( -getArg( 0 ) ) );
    m_Q = EvtComplex( cos( getArg( 0 ) ), sin( getArg( 0 ) ) );
    m_D = EvtComplex( getArg( 6 ) * cos( getArg( 7 ) ),
                      getArg( 6 ) * sin( getArg( 7 ) ) );
    m_Im = EvtComplex( 0.0, 1.0 );

    m_A = EvtComplex( getArg( 2 ) * cos( getArg( 3 ) ),
                      getArg( 2 ) * sin( getArg( 3 ) ) );
    m_Abar = EvtComplex( getArg( 4 ) * cos( getArg( 5 ) ),
                         getArg( 4 ) * sin( getArg( 5 ) ) );
}

void EvtSSSCPT::initProbMax()
{
    const double maxProb = 2.0 * abs2( m_A ) + 4.0 * abs2( m_Abar ) * abs2( m_D );
    setProbMax( maxProb );
}

void EvtSSSCPT::decay( EvtParticle* p )
{
    //added by Lange Jan4,2000
    static const EvtId B0 = EvtPDL::getId( "B0" );
    static const EvtId B0B = EvtPDL::getId( "anti-B0" );

    double t;
    EvtId other_b;

    EvtCPUtil::getInstance()->OtherB( p, t, other_b, 0.5 );

    p->initializePhaseSpace( getNDaug(), getDaugs() );

    EvtComplex amp;

    if ( other_b == B0B ) {
        amp = m_A * cos( getArg( 1 ) * t / ( 2 * EvtConst::c ) ) +
              m_Im * sin( getArg( 1 ) * t / ( 2 * EvtConst::c ) ) *
                  ( m_Q / m_P * m_A + 2.0 * m_D * m_Abar );
    }
    if ( other_b == B0 ) {
        amp = m_Abar * cos( getArg( 1 ) * t / ( 2 * EvtConst::c ) ) +
              m_Im * sin( getArg( 1 ) * t / ( 2 * EvtConst::c ) ) *
                  ( m_P / m_Q * m_A - 2.0 * m_D * m_Abar );
    }

    vertex( amp );

    return;
}
