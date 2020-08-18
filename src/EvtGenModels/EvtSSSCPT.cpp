
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
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"

#include <stdlib.h>
#include <string>

std::string EvtSSSCPT::getName()
{
    return "SSS_CPT";
}

EvtDecayBase* EvtSSSCPT::clone()
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
    P = EvtComplex( cos( -getArg( 0 ) ), sin( -getArg( 0 ) ) );
    Q = EvtComplex( cos( getArg( 0 ) ), sin( getArg( 0 ) ) );
    D = EvtComplex( getArg( 6 ) * cos( getArg( 7 ) ),
                    getArg( 6 ) * sin( getArg( 7 ) ) );
    Im = EvtComplex( 0.0, 1.0 );

    A = EvtComplex( getArg( 2 ) * cos( getArg( 3 ) ),
                    getArg( 2 ) * sin( getArg( 3 ) ) );
    Abar = EvtComplex( getArg( 4 ) * cos( getArg( 5 ) ),
                       getArg( 4 ) * sin( getArg( 5 ) ) );
}

void EvtSSSCPT::initProbMax()
{
    const double maxProb = 2.0 * abs2( A ) + 4.0 * abs2( Abar ) * abs2( D );
    setProbMax( maxProb );
}

void EvtSSSCPT::decay( EvtParticle* p )
{
    //added by Lange Jan4,2000
    static EvtId B0 = EvtPDL::getId( "B0" );
    static EvtId B0B = EvtPDL::getId( "anti-B0" );

    double t;
    EvtId other_b;

    EvtCPUtil::getInstance()->OtherB( p, t, other_b, 0.5 );

    p->initializePhaseSpace( getNDaug(), getDaugs() );

    EvtComplex amp;

    if ( other_b == B0B ) {
        amp = A * cos( getArg( 1 ) * t / ( 2 * EvtConst::c ) ) +
              Im * sin( getArg( 1 ) * t / ( 2 * EvtConst::c ) ) *
                  ( Q / P * A + 2.0 * D * Abar );
    }
    if ( other_b == B0 ) {
        amp = Abar * cos( getArg( 1 ) * t / ( 2 * EvtConst::c ) ) +
              Im * sin( getArg( 1 ) * t / ( 2 * EvtConst::c ) ) *
                  ( P / Q * A - 2.0 * D * Abar );
    }

    vertex( amp );

    return;
}
