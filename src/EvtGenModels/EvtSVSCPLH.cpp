
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

#include "EvtGenModels/EvtSVSCPLH.hh"

#include "EvtGenBase/EvtCPUtil.hh"
#include "EvtGenBase/EvtConst.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtVector4C.hh"

#include <stdlib.h>
#include <string>
using std::endl;

std::string EvtSVSCPLH::getName()
{
    return "SVS_CPLH";
}

EvtDecayBase* EvtSVSCPLH::clone()
{
    return new EvtSVSCPLH;
}

void EvtSVSCPLH::init()
{
    // check that there are 8 arguments
    checkNArg( 8 );
    checkNDaug( 2 );

    checkSpinParent( EvtSpinType::SCALAR );

    checkSpinDaughter( 0, EvtSpinType::VECTOR );
    checkSpinDaughter( 1, EvtSpinType::SCALAR );

    static double ctau = EvtPDL::getctau( EvtPDL::getId( "B0" ) );

    // hbar/s
    m_dm = getArg( 0 );
    m_dgamma = EvtConst::c * getArg( 1 ) / ctau;

    m_qop = getArg( 2 ) * exp( EvtComplex( 0.0, getArg( 3 ) ) );

    m_poq = 1.0 / m_qop;

    m_Af = getArg( 4 ) * exp( EvtComplex( 0.0, getArg( 5 ) ) );
    m_Abarf = getArg( 6 ) * exp( EvtComplex( 0.0, getArg( 7 ) ) );

    if ( verbose() ) {
        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << ":EvtSVSCPLH:dm=" << m_dm << endl;
        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << ":EvtSVSCPLH:dGamma=" << m_dgamma << endl;
        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << ":EvtSVSCPLH:q/p=" << m_qop << endl;
        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << ":EvtSVSCPLH:Af=" << m_Af << endl;
        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << ":EvtSVSCPLH:Abarf=" << m_Abarf << endl;
    }
}

void EvtSVSCPLH::initProbMax()
{
    //This is probably not quite right, but it should do as a start...
    //Anders

    setProbMax( 4.0 * ( getArg( 4 ) * getArg( 4 ) + getArg( 6 ) * getArg( 6 ) ) );
}

void EvtSVSCPLH::decay( EvtParticle* p )
{
    p->initializePhaseSpace( getNDaug(), getDaugs() );

    static EvtId B0 = EvtPDL::getId( "B0" );
    static EvtId B0B = EvtPDL::getId( "anti-B0" );

    double t;
    EvtId other_b;

    EvtCPUtil::getInstance()->OtherB( p, t, other_b, 0.5 );

    //convert time from mm to seconds
    t /= EvtConst::c;

    //sign convention is dm=Mheavy-Mlight
    //                   dGamma=Gammalight-Gammaheavy
    //such that in the standard model both of these are positive.
    EvtComplex gp =
        0.5 * ( exp( EvtComplex( 0.25 * t * m_dgamma, -0.5 * t * m_dm ) ) +
                exp( EvtComplex( -0.25 * t * m_dgamma, 0.5 * t * m_dm ) ) );
    EvtComplex gm =
        0.5 * ( exp( EvtComplex( 0.25 * t * m_dgamma, -0.5 * t * m_dm ) ) -
                exp( EvtComplex( -0.25 * t * m_dgamma, 0.5 * t * m_dm ) ) );

    EvtComplex amp;

    if ( other_b == B0B ) {
        amp = gp * m_Af + m_qop * gm * m_Abarf;
    } else if ( other_b == B0 ) {
        amp = gp * m_Abarf + m_poq * gm * m_Af;
    } else {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "other_b was not B0 or B0B!" << endl;
        ::abort();
    }

    EvtVector4R p4_parent = p->getP4Restframe();
    ;

    double norm = p->getDaug( 0 )->mass() /
                  ( p->getDaug( 0 )->getP4().d3mag() * p4_parent.mass() );

    EvtParticle* v = p->getDaug( 0 );

    vertex( 0, amp * norm * ( p4_parent * ( v->epsParent( 0 ) ) ) );
    vertex( 1, amp * norm * ( p4_parent * ( v->epsParent( 1 ) ) ) );
    vertex( 2, amp * norm * ( p4_parent * ( v->epsParent( 2 ) ) ) );

    return;
}
