
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

#include "EvtGenModels/EvtSVSCPiso.hh"

#include "EvtGenBase/EvtCPUtil.hh"
#include "EvtGenBase/EvtConst.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtVector4C.hh"

#include <stdlib.h>
#include <string>

std::string EvtSVSCPiso::getName()
{
    return "SVS_CP_ISO";
}

EvtDecayBase* EvtSVSCPiso::clone()
{
    return new EvtSVSCPiso;
}

void EvtSVSCPiso::init()
{
    // check that there are 27 arguments
    checkNArg( 27 );
    checkNDaug( 2 );

    checkSpinParent( EvtSpinType::SCALAR );

    checkSpinDaughter( 0, EvtSpinType::VECTOR );
    checkSpinDaughter( 1, EvtSpinType::SCALAR );

    // Set amplitude coefficients
    setAmpCoeffs();
    // Calculate amplitude terms
    calcAmpTerms();
}

void EvtSVSCPiso::setAmpCoeffs()
{
    m_Tp0 = EvtComplex( getArg( 3 ) * cos( getArg( 4 ) ),
                        getArg( 3 ) * sin( getArg( 4 ) ) );
    m_Tp0_bar = EvtComplex( getArg( 5 ) * cos( getArg( 6 ) ),
                            getArg( 5 ) * sin( getArg( 6 ) ) );
    m_T0p = EvtComplex( getArg( 7 ) * cos( getArg( 8 ) ),
                        getArg( 7 ) * sin( getArg( 8 ) ) );
    m_T0p_bar = EvtComplex( getArg( 9 ) * cos( getArg( 10 ) ),
                            getArg( 9 ) * sin( getArg( 10 ) ) );
    m_Tpm = EvtComplex( getArg( 11 ) * cos( getArg( 12 ) ),
                        getArg( 11 ) * sin( getArg( 12 ) ) );
    m_Tpm_bar = EvtComplex( getArg( 13 ) * cos( getArg( 14 ) ),
                            getArg( 13 ) * sin( getArg( 14 ) ) );
    m_Tmp = EvtComplex( getArg( 15 ) * cos( getArg( 16 ) ),
                        getArg( 15 ) * sin( getArg( 16 ) ) );
    m_Tmp_bar = EvtComplex( getArg( 17 ) * cos( getArg( 18 ) ),
                            getArg( 17 ) * sin( getArg( 18 ) ) );
    m_P0 = EvtComplex( getArg( 19 ) * cos( getArg( 20 ) ),
                       getArg( 19 ) * sin( getArg( 20 ) ) );
    m_P0_bar = EvtComplex( getArg( 21 ) * cos( getArg( 22 ) ),
                           getArg( 21 ) * sin( getArg( 22 ) ) );
    m_P1 = EvtComplex( getArg( 23 ) * cos( getArg( 24 ) ),
                       getArg( 23 ) * sin( getArg( 24 ) ) );
    m_P1_bar = EvtComplex( getArg( 25 ) * cos( getArg( 26 ) ),
                           getArg( 25 ) * sin( getArg( 26 ) ) );
}

void EvtSVSCPiso::initProbMax()
{
    const double max1 = abs2( m_A_f ) + abs2( m_Abar_f );
    const double max2 = abs2( m_A_fbar ) + abs2( m_Abar_fbar );
    // Amplitude has momentum normalisation that roughly scales with (parent mass)/2
    // so probability will scale with 0.25 * parenMassSq. Use 0.3 * parMassSq
    // in case we get larger normalisation values
    const double parMass = EvtPDL::getMeanMass( getParentId() );
    const double max = 0.3 * parMass * parMass * ( max1 + max2 );
    setProbMax( max );
}

void EvtSVSCPiso::decay( EvtParticle* p )
{
    //added by Lange Jan4,2000
    static EvtId B0 = EvtPDL::getId( "B0" );
    static EvtId B0B = EvtPDL::getId( "anti-B0" );

    double t;
    EvtId other_b;

    int first_time = 0;
    int flip = 0;
    EvtId ds[2];

    // Randomly generate the tag (B0 or B0B)
    const double tag = EvtRandom::Flat( 0.0, 1.0 );
    if ( tag < 0.5 ) {
        EvtCPUtil::getInstance()->OtherB( p, t, other_b, 1.0 );
        other_b = B0;
    } else {
        EvtCPUtil::getInstance()->OtherB( p, t, other_b, 0.0 );
        other_b = B0B;
    }

    if ( p->getNDaug() == 0 )
        first_time = 1;

    if ( first_time ) {
        if ( EvtRandom::Flat( 0.0, 1.0 ) < getArg( 2 ) )
            flip = 1;
    } else {
        if ( getDaug( 0 ) != p->getDaug( 0 )->getId() )
            flip = 1;
    }

    if ( !flip ) {
        ds[0] = getDaug( 0 );
        ds[1] = getDaug( 1 );
    } else {
        ds[0] = EvtPDL::chargeConj( getDaug( 0 ) );
        ds[1] = EvtPDL::chargeConj( getDaug( 1 ) );
    }

    p->initializePhaseSpace( getNDaug(), ds );

    EvtParticle *v, *s;
    v = p->getDaug( 0 );
    s = p->getDaug( 1 );

    EvtComplex amp;

    if ( m_charged == 0 ) {
        if ( !flip ) {
            if ( other_b == B0B ) {
                amp = m_A_f * cos( getArg( 1 ) * t / ( 2 * EvtConst::c ) ) +
                      EvtComplex( cos( -2.0 * getArg( 0 ) ),
                                  sin( -2.0 * getArg( 0 ) ) ) *
                          EvtComplex( 0.0, 1.0 ) * m_Abar_f *
                          sin( getArg( 1 ) * t / ( 2 * EvtConst::c ) );
            }
            if ( other_b == B0 ) {
                amp = m_A_f *
                          EvtComplex( cos( 2.0 * getArg( 0 ) ),
                                      sin( 2.0 * getArg( 0 ) ) ) *
                          EvtComplex( 0.0, 1.0 ) *
                          sin( getArg( 1 ) * t / ( 2 * EvtConst::c ) ) +
                      m_Abar_f * cos( getArg( 1 ) * t / ( 2 * EvtConst::c ) );
            }
        } else {
            if ( other_b == B0B ) {
                amp = m_A_fbar * cos( getArg( 1 ) * t / ( 2 * EvtConst::c ) ) +
                      EvtComplex( cos( -2.0 * getArg( 0 ) ),
                                  sin( -2.0 * getArg( 0 ) ) ) *
                          EvtComplex( 0.0, 1.0 ) * m_Abar_fbar *
                          sin( getArg( 1 ) * t / ( 2 * EvtConst::c ) );
            }
            if ( other_b == B0 ) {
                amp = m_A_fbar *
                          EvtComplex( cos( 2.0 * getArg( 0 ) ),
                                      sin( 2.0 * getArg( 0 ) ) ) *
                          EvtComplex( 0.0, 1.0 ) *
                          sin( getArg( 1 ) * t / ( 2 * EvtConst::c ) ) +
                      m_Abar_fbar * cos( getArg( 1 ) * t / ( 2 * EvtConst::c ) );
            }
        }

    } else {
        amp = m_A_f;
    }

    const EvtVector4R p4_parent = v->getP4() + s->getP4();
    const double norm = 1.0 / v->getP4().d3mag();

    vertex( 0, amp * norm * p4_parent * ( v->epsParent( 0 ) ) );
    vertex( 1, amp * norm * p4_parent * ( v->epsParent( 1 ) ) );
    vertex( 2, amp * norm * p4_parent * ( v->epsParent( 2 ) ) );

    return;
}

void EvtSVSCPiso::calcAmpTerms()
{
    const int Q1 = EvtPDL::chg3( getDaug( 0 ) );
    const int Q2 = EvtPDL::chg3( getDaug( 1 ) );

    //***********************charged modes****************************

    if ( Q1 > 0 && Q2 == 0 ) {
        //V+ S0, so T+0 + 2 P1

        m_charged = 1;
        m_A_f = m_Tp0 + 2.0 * m_P1;
    }

    if ( Q1 < 0 && Q2 == 0 ) {
        //V- S0, so T+0_bar + 2P1_bar

        m_charged = 1;
        m_A_f = m_Tp0_bar + 2.0 * m_P1_bar;
    }

    if ( Q1 == 0 && Q2 > 0 ) {
        //V0 S+, so T0+ - 2 P1

        m_charged = 1;
        m_A_f = m_T0p - 2.0 * m_P1;
    }

    if ( Q1 == 0 && Q2 < 0 ) {
        //V0 S-, so T0+_bar - 2 P1_bar

        m_charged = 1;
        m_A_f = m_T0p_bar - 2.0 * m_P1_bar;
    }

    //***********************neutral modes***************************

    //V+ S-, so Af = T+- + P1 + P0
    m_Apm = m_Tpm + m_P1 + m_P0;
    m_Apm_bar = m_Tpm_bar + m_P1_bar + m_P0_bar;

    //V- S+, so Af = T-+ - P1 + P0
    m_Amp = m_Tmp - m_P1 + m_P0;
    m_Amp_bar = m_Tmp_bar - m_P1_bar + m_P0;

    if ( Q1 > 0 && Q2 < 0 ) {
        //V+ S-
        m_charged = 0;
        m_A_f = m_Apm;
        m_Abar_f = m_Apm_bar;
        m_A_fbar = m_Amp;
        m_Abar_fbar = m_Amp_bar;
    }

    if ( Q1 < 0 && Q2 > 0 ) {
        //V- S+
        m_charged = 0;
        m_A_f = m_Amp;
        m_Abar_f = m_Amp_bar;
        m_A_fbar = m_Apm;
        m_Abar_fbar = m_Apm_bar;
    }

    if ( Q1 == 0 && Q2 == 0 ) {
        //V0 S0
        m_charged = 0;
        m_A_f = m_T0p + m_Tp0 - m_Tpm - m_Tmp - 2.0 * m_P0;
        m_Abar_f = m_T0p_bar + m_Tp0_bar - m_Tpm_bar - m_Tmp_bar - 2.0 * m_P0_bar;
        m_A_fbar = m_A_f;
        m_Abar_fbar = m_Abar_f;
    }
}
