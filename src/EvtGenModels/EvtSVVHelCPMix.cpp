
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

#include "EvtGenModels/EvtSVVHelCPMix.hh"

#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtTensor3C.hh"
#include "EvtGenBase/EvtTensor4C.hh"
#include "EvtGenBase/EvtVector3C.hh"
#include "EvtGenBase/EvtVector3R.hh"
#include "EvtGenBase/EvtVector4C.hh"

#include <ctype.h>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <string>

std::string EvtSVVHelCPMix::getName()
{
    return "SVVHELCPMIX";
}

EvtDecayBase* EvtSVVHelCPMix::clone()
{
    return new EvtSVVHelCPMix;
}

void EvtSVVHelCPMix::init()
{
    // check that there are 12 arguments
    checkNArg( 12 );
    checkNDaug( 2 );

    checkSpinParent( EvtSpinType::SCALAR );

    checkSpinDaughter( 0, EvtSpinType::VECTOR );
    checkSpinDaughter( 1, EvtSpinType::VECTOR );

    m_hp = EvtComplex( getArg( 0 ) * cos( getArg( 1 ) ),
                       getArg( 0 ) * sin( getArg( 1 ) ) );
    m_h0 = EvtComplex( getArg( 2 ) * cos( getArg( 3 ) ),
                       getArg( 2 ) * sin( getArg( 3 ) ) );
    m_hm = EvtComplex( getArg( 4 ) * cos( getArg( 5 ) ),
                       getArg( 4 ) * sin( getArg( 5 ) ) );
    m_averageM = getArg( 6 );
    m_deltaM = getArg( 7 );
    m_gamma = getArg( 8 );
    m_deltagamma = getArg( 9 );
    m_weakmixingphase = EvtComplex( cos( getArg( 10 ) ), sin( getArg( 10 ) ) );
    m_weakdirectphase = EvtComplex( cos( getArg( 11 ) ), sin( getArg( 11 ) ) );
}

void EvtSVVHelCPMix::initProbMax()
{
    setProbMax( getArg( 0 ) * getArg( 0 ) + getArg( 2 ) * getArg( 2 ) +
                getArg( 4 ) * getArg( 4 ) );
}

void EvtSVVHelCPMix::decay( EvtParticle* p )
{
    EvtParticle* parent = p;
    EvtAmp& amp = m_amp2;
    EvtId n_v1 = getDaug( 0 );
    EvtId n_v2 = getDaug( 1 );

    //  Routine to decay a vector into a vector and scalar.  Started
    //  by ryd on Oct 17, 1996.
    // Modified by J.Catmore to take account of CP-violation and mixing

    int tndaug = 2;
    EvtId tdaug[2];
    EvtId Bs = EvtPDL::getId( "B_s0" );
    EvtId antiBs = EvtPDL::getId( "anti-B_s0" );
    tdaug[0] = n_v1;
    tdaug[1] = n_v2;

    // Phase space and kinematics

    parent->initializePhaseSpace( tndaug, tdaug );

    EvtParticle *v1, *v2;
    v1 = parent->getDaug( 0 );
    v2 = parent->getDaug( 1 );

    EvtVector4R momv1 = v1->getP4();

    EvtVector3R v1dir( momv1.get( 1 ), momv1.get( 2 ), momv1.get( 3 ) );
    v1dir = v1dir / v1dir.d3mag();

    // Definition of quantities used in construction of complex amplitudes:

    EvtTensor3C M;    // Tensor as defined in EvtGen manual, equ 117
    EvtComplex a, b,
        c;    // Helicity amplitudes; EvtGen manual eqns 126-128, also see Phys Lett B 369 p144-150 eqn 15
    //EvtComplex deltamu = EvtComplex(m_deltaM, -0.5*m_deltagamma); // See Phys Rev D 34 p1404

    // conversion from times in mm/c to natural units [GeV]^-1
    double t = ( ( parent->getLifetime() ) / 2.998e11 ) * 6.58e-25;

    // The following two quantities defined in Phys Rev D 34 p1404
    EvtComplex fplus =
        EvtComplex( cos( m_averageM * t ), -1. * sin( m_averageM * t ) ) *
        exp( -( m_gamma / 2.0 ) * t ) *
        ( cos( 0.5 * m_deltaM * t ) * cosh( 0.25 * m_deltagamma * t ) +
          EvtComplex( 0.0, sin( 0.5 * m_deltaM * t ) *
                               sinh( 0.25 * m_deltagamma * t ) ) );
    EvtComplex fminus =
        EvtComplex( cos( m_averageM * t ), -1. * sin( m_averageM * t ) ) *
        exp( -( m_gamma / 2.0 ) * t ) * EvtComplex( 0.0, 1.0 ) *
        ( sin( 0.5 * m_deltaM * t ) * cosh( 0.25 * m_deltagamma * t ) -
          EvtComplex( 0.0, 1.0 ) * sinh( 0.25 * m_deltagamma * t ) *
              cos( 0.5 * m_deltaM * t ) );

    // See EvtGen manual pp 106-107

    a = -0.5 * ( m_hp + m_hm );
    b = EvtComplex( 0.0, 0.5 ) * ( m_hp - m_hm );
    c = ( m_h0 + 0.5 * ( m_hp + m_hm ) );

    M = a * EvtTensor3C::id() + b * EvtGenFunctions::eps( v1dir ) +
        c * EvtGenFunctions::directProd( v1dir, v1dir );

    EvtVector3C t0 = M.cont1( v1->eps( 0 ).vec().conj() );
    EvtVector3C t1 = M.cont1( v1->eps( 1 ).vec().conj() );
    EvtVector3C t2 = M.cont1( v1->eps( 2 ).vec().conj() );

    EvtVector3C eps0 = v2->eps( 0 ).vec().conj();
    EvtVector3C eps1 = v2->eps( 1 ).vec().conj();
    EvtVector3C eps2 = v2->eps( 2 ).vec().conj();

    // We need two sets of equations, one for mesons which were in the Bs state at t=0, and another
    // for those which were in the antiBs state. Each equation consists of a sum of amplitudes - mod-squaring gives the interference terms.

    EvtComplex amplSum00, amplSum01, amplSum02;
    EvtComplex amplSum10, amplSum11, amplSum12;
    EvtComplex amplSum20, amplSum21, amplSum22;

    // First the Bs state:

    if ( parent->getId() == Bs ) {
        amplSum00 = ( fplus * m_weakdirectphase * t0 * eps0 ) +
                    ( fminus * ( 1.0 / m_weakdirectphase ) * m_weakmixingphase *
                      t0 * eps0 );
        amplSum01 = ( fplus * m_weakdirectphase * t0 * eps1 ) +
                    ( fminus * ( 1.0 / m_weakdirectphase ) * m_weakmixingphase *
                      t0 * eps1 );
        amplSum02 = ( fplus * m_weakdirectphase * t0 * eps2 ) +
                    ( fminus * ( 1.0 / m_weakdirectphase ) * m_weakmixingphase *
                      t0 * eps2 );

        amplSum10 = ( fplus * m_weakdirectphase * t1 * eps0 ) +
                    ( fminus * ( 1.0 / m_weakdirectphase ) * m_weakmixingphase *
                      t1 * eps0 );
        amplSum11 = ( fplus * m_weakdirectphase * t1 * eps1 ) +
                    ( fminus * ( 1.0 / m_weakdirectphase ) * m_weakmixingphase *
                      t1 * eps1 );
        amplSum12 = ( fplus * m_weakdirectphase * t1 * eps2 ) +
                    ( fminus * ( 1.0 / m_weakdirectphase ) * m_weakmixingphase *
                      t1 * eps2 );

        amplSum20 = ( fplus * m_weakdirectphase * t2 * eps0 ) +
                    ( fminus * ( 1.0 / m_weakdirectphase ) * m_weakmixingphase *
                      t2 * eps0 );
        amplSum21 = ( fplus * m_weakdirectphase * t2 * eps1 ) +
                    ( fminus * ( 1.0 / m_weakdirectphase ) * m_weakmixingphase *
                      t2 * eps1 );
        amplSum22 = ( fplus * m_weakdirectphase * t2 * eps2 ) +
                    ( fminus * ( 1.0 / m_weakdirectphase ) * m_weakmixingphase *
                      t2 * eps2 );
    }

    // Now the anti-Bs state:

    if ( parent->getId() == antiBs ) {
        amplSum00 = ( fminus * m_weakdirectphase * ( 1.0 / m_weakmixingphase ) *
                      t0 * eps0 ) +
                    ( fplus * ( 1.0 / m_weakdirectphase ) * t0 * eps0 );
        amplSum01 = ( fminus * m_weakdirectphase * ( 1.0 / m_weakmixingphase ) *
                      t0 * eps1 ) +
                    ( fplus * ( 1.0 / m_weakdirectphase ) * t0 * eps1 );
        amplSum02 = ( fminus * m_weakdirectphase * ( 1.0 / m_weakmixingphase ) *
                      t0 * eps2 ) +
                    ( fplus * ( 1.0 / m_weakdirectphase ) * t0 * eps2 );

        amplSum10 = ( fminus * m_weakdirectphase * ( 1.0 / m_weakmixingphase ) *
                      t1 * eps0 ) +
                    ( fplus * ( 1.0 / m_weakdirectphase ) * t1 * eps0 );
        amplSum11 = ( fminus * m_weakdirectphase * ( 1.0 / m_weakmixingphase ) *
                      t1 * eps1 ) +
                    ( fplus * ( 1.0 / m_weakdirectphase ) * t1 * eps1 );
        amplSum12 = ( fminus * m_weakdirectphase * ( 1.0 / m_weakmixingphase ) *
                      t1 * eps2 ) +
                    ( fplus * ( 1.0 / m_weakdirectphase ) * t1 * eps2 );

        amplSum20 = ( fminus * m_weakdirectphase * ( 1.0 / m_weakmixingphase ) *
                      t2 * eps0 ) +
                    ( fplus * ( 1.0 / m_weakdirectphase ) * t2 * eps0 );
        amplSum21 = ( fminus * m_weakdirectphase * ( 1.0 / m_weakmixingphase ) *
                      t2 * eps1 ) +
                    ( fplus * ( 1.0 / m_weakdirectphase ) * t2 * eps1 );
        amplSum22 = ( fminus * m_weakdirectphase * ( 1.0 / m_weakmixingphase ) *
                      t2 * eps2 ) +
                    ( fplus * ( 1.0 / m_weakdirectphase ) * t2 * eps2 );
    }

    // Now set the amplitudes
    amp.vertex( 0, 0, amplSum00 );
    amp.vertex( 0, 1, amplSum01 );
    amp.vertex( 0, 2, amplSum02 );

    amp.vertex( 1, 0, amplSum10 );
    amp.vertex( 1, 1, amplSum11 );
    amp.vertex( 1, 2, amplSum12 );

    amp.vertex( 2, 0, amplSum20 );
    amp.vertex( 2, 1, amplSum21 );
    amp.vertex( 2, 2, amplSum22 );

    return;
}

std::string EvtSVVHelCPMix::getParamName( int i )
{
    switch ( i ) {
        case 0:
            return "plusHelAmp";
        case 1:
            return "plusHelAmpPhase";
        case 2:
            return "zeroHelAmp";
        case 3:
            return "zeroHelAmpPhase";
        case 4:
            return "minusHelAmp";
        case 5:
            return "minusHelAmpPhase";
        case 6:
            return "averageM";
        case 7:
            return "deltaM";
        case 8:
            return "gamma";
        case 9:
            return "deltaGamma";
        case 10:
            return "weakMixPhase";
        case 11:
            return "weakDirectPhase";
        default:
            return "";
    }
}

std::string EvtSVVHelCPMix::getParamDefault( int i )
{
    switch ( i ) {
        case 0:
            return "1.0";
        case 1:
            return "0.0";
        case 2:
            return "1.0";
        case 3:
            return "0.0";
        case 4:
            return "1.0";
        case 5:
            return "0.0";
        default:
            return "";
    }
}
