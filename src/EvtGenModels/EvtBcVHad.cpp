
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

#include "EvtGenModels/EvtBcVHad.hh"

#include "EvtGenBase/EvtIdSet.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtSpinType.hh"
#include "EvtGenBase/EvtTensor4C.hh"
#include "EvtGenBase/EvtVector4R.hh"

#include "EvtGenModels/EvtBCVFF2.hh"
#include "EvtGenModels/EvtWHad.hh"

#include <iostream>

std::string EvtBcVHad::getName()
{
    return "BC_VHAD";
}

EvtDecayBase* EvtBcVHad::clone()
{
    return new EvtBcVHad;
}

//======================================================

void EvtBcVHad::init()
{
    // The following decay m_outCodes are supported:
    // 1:  B_c+ -> V pi+
    // 2:  B_c+ -> V pi+ pi0
    // 3:  B_c+ -> V 2pi+ pi-
    // 4:  B_c+ -> V 2pi+ pi- pi0 (not implemented)
    // 5:  B_c+ -> V 3pi+ 2pi-
    // 6:  B_c+ -> V K+ K- pi+
    // 7:  B_c+ -> V K+ pi+ pi-
    // 8:  B_c+ -> V K_S0 K+
    // 9:  B_c+ -> V K+ K- 2pi+ pi-
    // 10: B_c+ -> V 4pi+ 3pi-
    // 11: B_c+ -> V K+ 2pi+ 2pi-

    checkNArg( 1 );
    checkSpinParent( EvtSpinType::SCALAR );
    checkSpinDaughter( 0, EvtSpinType::VECTOR );
    for ( int i = 1; i <= ( getNDaug() - 1 ); i++ ) {
        checkSpinDaughter( i, EvtSpinType::SCALAR );
    }

    m_idVector = getDaug( 0 ).getId();
    m_whichFit = int( getArg( 0 ) + 0.1 );
    m_FFModel = std::make_unique<EvtBCVFF2>( m_idVector, m_whichFit );

    m_WCurr = std::make_unique<EvtWHad>();

    // Determine the code of the final hadronic state
    parseDecay();
}

//======================================================

void EvtBcVHad::parseDecay()
{
    const EvtIdSet BcPlusID{ "B_c+" }, BcMinusID{ "B_c-" };
    const EvtIdSet theK{ "K+", "K-", "K_S0" };
    const int cMode = BcMinusID.contains( getParentId() );

    const EvtIdSet PiPlusID{ cMode == 0 ? "pi+" : "pi-" };
    const EvtIdSet PiMinusID{ cMode == 0 ? "pi-" : "pi+" };
    const EvtIdSet PiZeroID{ "pi0" };
    const EvtIdSet KPlusID{ cMode == 0 ? "K+" : "K-" };
    const EvtIdSet KMinusID{ cMode == 0 ? "K-" : "K+" };
    EvtGenReport( EVTGEN_INFO, "EvtBcVHad" )
        << "parentId = " << getParentId() << std::endl;

    int PiPlusFound = 0, PiMinusFound = 0, PiZeroFound = 0, KPlusFound = 0,
        KMinusFound = 0;
    for ( int iDaughter = 0; iDaughter < getNDaug(); ++iDaughter ) {
        const EvtId daugId = getDaug( iDaughter );
        EvtGenReport( EVTGEN_INFO, "EvtBcVHad" )
            << "iDaughter = " << iDaughter
            << " id = " << getDaug( iDaughter ).getName() << std::endl;
        if ( PiPlusID.contains( daugId ) && PiPlusFound < 4 ) {
            m_iPiPlus[PiPlusFound] = iDaughter;
            PiPlusFound++;
        } else if ( PiMinusID.contains( daugId ) && PiMinusFound < 4 ) {
            m_iPiMinus[PiMinusFound] = iDaughter;
            PiMinusFound++;
        } else if ( PiZeroID.contains( daugId ) && PiZeroFound < 4 ) {
            m_iPiZero[PiZeroFound] = iDaughter;
            PiZeroFound++;
        } else if ( KPlusID.contains( daugId ) && KPlusFound < 4 ) {
            m_iKPlus[KPlusFound] = iDaughter;
            KPlusFound++;
        } else if ( KMinusID.contains( daugId ) && KMinusFound < 4 ) {
            m_iKMinus[KMinusFound] = iDaughter;
            KMinusFound++;
        }
    }

    if ( getNDaug() == 2 && PiPlusFound == 1 ) {
        m_outCode = 1;    // pi+
    } else if ( getNDaug() == 3 && PiPlusFound == 1 && PiZeroFound == 1 ) {
        m_outCode = 2;    // pi+ pi0
    } else if ( getNDaug() == 4 && PiPlusFound == 2 && PiMinusFound == 1 ) {
        m_outCode = 3;    // pi+ pi+ pi-
    } else if ( getNDaug() == 5 && PiPlusFound == 2 && PiMinusFound == 1 &&
                PiZeroFound == 1 ) {
        m_outCode = 4;    // pi+ pi+ pi- pi0
    } else if ( getNDaug() == 6 && PiPlusFound == 3 && PiMinusFound == 2 ) {
        m_outCode = 5;    // 5pi
    } else if ( getNDaug() == 4 && KPlusFound == 1 && KMinusFound == 1 &&
                PiPlusFound == 1 ) {
        m_outCode = 6;    // KKpi
    } else if ( getNDaug() == 4 && KPlusFound == 1 && PiPlusFound == 1 &&
                PiMinusFound == 1 ) {
        m_outCode = 7;    // K+ pi+ pi-
    } else if ( getNDaug() == 3 && theK.contains( getDaug( 1 ) ) &&
                theK.contains( getDaug( 2 ) ) ) {
        m_outCode = 8;    // KK
    } else if ( getNDaug() == 6 && KPlusFound == 1 && KMinusFound == 1 &&
                PiPlusFound == 2 && PiMinusFound == 1 ) {
        m_outCode = 9;    // K+ K- pi+ pi+ pi-
    } else if ( getNDaug() == 8 && PiPlusFound == 4 && PiMinusFound == 3 ) {
        m_outCode = 10;    // 7pi
    } else if ( getNDaug() == 6 && KPlusFound == 1 && PiPlusFound == 2 &&
                PiMinusFound == 2 ) {
        m_outCode = 11;    // K+ pi+ pi+ pi- pi-
    } else {
        EvtGenReport( EVTGEN_ERROR, "EvtBcVHad" )
            << "Init: unknown decay" << std::endl;
    }

    EvtGenReport( EVTGEN_INFO, "EvtBcVHad" )
        << "m_outCode = " << m_outCode << ", m_whichFit = " << m_whichFit
        << std::endl;
    for ( int i = 0; i < 4; ++i ) {
        EvtGenReport( EVTGEN_INFO, "EvtBcVHad" )
            << " i = " << i << ", m_iPiPlus = " << m_iPiPlus[i]
            << ", m_iPiMinus = " << m_iPiMinus[i]
            << ", m_iPiZero = " << m_iPiZero[i] << ", m_iKPlus = " << m_iKPlus[i]
            << ", m_iKMinus = " << m_iKMinus[i] << std::endl;
    }
    EvtGenReport( EVTGEN_INFO, "EvtBcVHad" )
        << "PiPlusFound = " << PiPlusFound << ", PiMinusFound = " << PiMinusFound
        << ", PiZeroFound = " << PiZeroFound << ", KPlusFound = " << KPlusFound
        << ", KMinusFound = " << KMinusFound << std::endl;
}

//======================================================

void EvtBcVHad::initProbMax()
{
    if ( m_outCode == 1 ) {
        if ( m_idVector == EvtPDL::getId( "J/psi" ).getId() &&
             m_whichFit == 1 && getNDaug() == 2 )
            setProbMax( 500. );
        else if ( m_idVector == EvtPDL::getId( "J/psi" ).getId() &&
                  m_whichFit == 2 && getNDaug() == 2 )
            setProbMax( 300. );
        else if ( m_idVector == EvtPDL::getId( "psi(2S)" ).getId() &&
                  m_whichFit == 1 && getNDaug() == 2 )
            setProbMax( 17. );
        else if ( m_idVector == EvtPDL::getId( "psi(2S)" ).getId() &&
                  m_whichFit == 2 && getNDaug() == 2 )
            setProbMax( 40. );
    } else if ( m_outCode == 2 ) {
        if ( m_idVector == EvtPDL::getId( "J/psi" ).getId() &&
             m_whichFit == 1 && getNDaug() == 3 )
            setProbMax( 10950. );
        else if ( m_idVector == EvtPDL::getId( "J/psi" ).getId() &&
                  m_whichFit == 2 && getNDaug() == 3 )
            setProbMax( 4200. );
        else if ( m_idVector == EvtPDL::getId( "psi(2S)" ).getId() &&
                  m_whichFit == 1 && getNDaug() == 3 )
            setProbMax( 500. );
        else if ( m_idVector == EvtPDL::getId( "psi(2S)" ).getId() &&
                  m_whichFit == 2 && getNDaug() == 3 )
            setProbMax( 700. );
    } else if ( m_outCode == 3 ) {
        if ( m_idVector == EvtPDL::getId( "J/psi" ).getId() &&
             m_whichFit == 1 && getNDaug() == 4 )
            setProbMax( 42000. );
        else if ( m_idVector == EvtPDL::getId( "J/psi" ).getId() &&
                  m_whichFit == 2 && getNDaug() == 4 )
            setProbMax( 90000. );
        else if ( m_idVector == EvtPDL::getId( "psi(2S)" ).getId() &&
                  m_whichFit == 1 && getNDaug() == 4 )
            setProbMax( 1660. );
        else if ( m_idVector == EvtPDL::getId( "psi(2S)" ).getId() &&
                  m_whichFit == 2 && getNDaug() == 4 )
            setProbMax( 2600. );
    } else if ( m_outCode == 5 ) {
        if ( m_idVector == EvtPDL::getId( "J/psi" ).getId() &&
             m_whichFit == 1 && getNDaug() == 6 )
            setProbMax( 720000. );
        else if ( m_idVector == EvtPDL::getId( "J/psi" ).getId() &&
                  m_whichFit == 2 && getNDaug() == 6 )
            setProbMax( 519753. );
        else if ( m_idVector == EvtPDL::getId( "psi(2S)" ).getId() &&
                  m_whichFit == 1 && getNDaug() == 6 )
            setProbMax( 40000. );
        else if ( m_idVector == EvtPDL::getId( "psi(2S)" ).getId() &&
                  m_whichFit == 2 && getNDaug() == 6 )
            setProbMax( 30000. );
    } else if ( m_outCode == 6 ) {
        if ( m_idVector == EvtPDL::getId( "J/psi" ).getId() && m_whichFit == 1 )
            setProbMax( 50000. );
        else if ( m_idVector == EvtPDL::getId( "J/psi" ).getId() &&
                  m_whichFit == 2 )
            setProbMax( 22000.0 );
        else if ( m_idVector == EvtPDL::getId( "psi(2S)" ).getId() &&
                  m_whichFit == 1 )
            setProbMax( 2300.0 );
        else if ( m_idVector == EvtPDL::getId( "psi(2S)" ).getId() &&
                  m_whichFit == 2 )
            setProbMax( 1700.00 );
    } else if ( m_outCode == 7 ) {
        if ( m_idVector == EvtPDL::getId( "J/psi" ).getId() && m_whichFit == 1 )
            setProbMax( 2.2e+06 );
        else if ( m_idVector == EvtPDL::getId( "J/psi" ).getId() &&
                  m_whichFit == 2 )
            setProbMax( 930000 );
        else if ( m_idVector == EvtPDL::getId( "psi(2S)" ).getId() &&
                  m_whichFit == 1 )
            setProbMax( 92000.0 );
        else if ( m_idVector == EvtPDL::getId( "psi(2S)" ).getId() &&
                  m_whichFit == 2 )
            setProbMax( 93000.0 );
    } else if ( m_outCode == 8 ) {
        if ( m_idVector == EvtPDL::getId( "J/psi" ).getId() && m_whichFit == 1 )
            setProbMax( 2e2 );
        else if ( m_idVector == EvtPDL::getId( "J/psi" ).getId() &&
                  m_whichFit == 2 )
            setProbMax( 80 );
        else if ( m_idVector == EvtPDL::getId( "psi(2S)" ).getId() &&
                  m_whichFit == 1 )
            setProbMax( 10 );
        else if ( m_idVector == EvtPDL::getId( "psi(2S)" ).getId() &&
                  m_whichFit == 2 )
            setProbMax( 10 );
    } else if ( m_outCode == 9 ) {
        if ( m_idVector == EvtPDL::getId( "J/psi" ).getId() && m_whichFit == 1 )
            setProbMax( 3e4 );
        else if ( m_idVector == EvtPDL::getId( "J/psi" ).getId() &&
                  m_whichFit == 2 )
            setProbMax( 18540 );
        else if ( m_idVector == EvtPDL::getId( "psi(2S)" ).getId() &&
                  m_whichFit == 1 )
            setProbMax( 0.15 * 1e4 );
        else if ( m_idVector == EvtPDL::getId( "psi(2S)" ).getId() &&
                  m_whichFit == 2 )
            setProbMax( 2 * 500 );
    } else if ( m_outCode == 10 ) {
        if ( m_idVector == EvtPDL::getId( "J/psi" ).getId() && m_whichFit == 1 )
            setProbMax( 2e6 );
        else if ( m_idVector == EvtPDL::getId( "J/psi" ).getId() &&
                  m_whichFit == 2 )
            setProbMax( 5e6 );
        else if ( m_idVector == EvtPDL::getId( "psi(2S)" ).getId() &&
                  m_whichFit == 1 )
            setProbMax( 1.5e5 );
        else if ( m_idVector == EvtPDL::getId( "psi(2S)" ).getId() &&
                  m_whichFit == 2 )
            setProbMax( 1e5 );
    } else if ( m_outCode == 11 ) {
        if ( m_idVector == EvtPDL::getId( "J/psi" ).getId() && m_whichFit == 1 )
            setProbMax( 2.5e8 );
        else if ( m_idVector == EvtPDL::getId( "J/psi" ).getId() &&
                  m_whichFit == 2 )
            setProbMax( 1.4e7 );
        else if ( m_idVector == EvtPDL::getId( "psi(2S)" ).getId() &&
                  m_whichFit == 1 )
            setProbMax( 2e6 );
        else if ( m_idVector == EvtPDL::getId( "psi(2S)" ).getId() &&
                  m_whichFit == 2 )
            setProbMax( 8e4 );
    } else {
        EvtGenReport( EVTGEN_ERROR, "EvtBcHad" )
            << "probmax: Have not yet implemented this final state in BC_VHAD model, m_outCode = "
            << m_outCode << std::endl;
        for ( int id = 0; id < ( getNDaug() - 1 ); id++ ) {
            EvtGenReport( EVTGEN_ERROR, "EvtBcVHad" )
                << "Daug " << id << " " << EvtPDL::name( getDaug( id ) ).c_str()
                << std::endl;
        }
    }
}

//======================================================

EvtVector4C EvtBcVHad::hardCurr( EvtParticle* parent ) const
{
    EvtVector4C hardCur;

    if ( m_outCode == 1 ) {
        // pi+
        hardCur = m_WCurr->WCurrent( parent->getDaug( m_iPiPlus[0] )->getP4() );

    } else if ( m_outCode == 2 ) {
        // pi+ pi0
        hardCur = m_WCurr->WCurrent( parent->getDaug( m_iPiPlus[0] )->getP4(),
                                     parent->getDaug( m_iPiZero[0] )->getP4() );

    } else if ( m_outCode == 3 ) {
        // pi+ pi+ pi-
        hardCur = m_WCurr->WCurrent( parent->getDaug( m_iPiPlus[0] )->getP4(),
                                     parent->getDaug( m_iPiPlus[1] )->getP4(),
                                     parent->getDaug( m_iPiMinus[0] )->getP4() );
    } else if ( m_outCode == 5 ) {
        // Bc -> psi pi+ pi+ pi- pi- pi+ from Kuhn & Was, hep-ph/0602162
        hardCur =
            m_WCurr->WCurrent_5pi( parent->getDaug( m_iPiPlus[0] )->getP4(),
                                   parent->getDaug( m_iPiPlus[1] )->getP4(),
                                   parent->getDaug( m_iPiPlus[2] )->getP4(),
                                   parent->getDaug( m_iPiMinus[0] )->getP4(),
                                   parent->getDaug( m_iPiMinus[1] )->getP4() );

    } else if ( m_outCode == 6 ) {
        // K+ K- pi+
        hardCur =
            m_WCurr->WCurrent_KKP( parent->getDaug( m_iKPlus[0] )->getP4(),
                                   parent->getDaug( m_iKMinus[0] )->getP4(),
                                   parent->getDaug( m_iPiPlus[0] )->getP4() );

    } else if ( m_outCode == 7 ) {
        // K+ pi+ pi-
        hardCur =
            m_WCurr->WCurrent_KPP( parent->getDaug( m_iKPlus[0] )->getP4(),
                                   parent->getDaug( m_iPiPlus[0] )->getP4(),
                                   parent->getDaug( m_iPiMinus[0] )->getP4() );

    } else if ( m_outCode == 8 ) {
        // K_S0 K+
        hardCur = m_WCurr->WCurrent_KSK( parent->getDaug( 1 )->getP4(),
                                         parent->getDaug( 2 )->getP4() );

    } else if ( m_outCode == 9 ) {
        // K+ K- pi+ pi+ pi-
        hardCur = m_WCurr->WCurrent_KKPPP(
            parent->getDaug( m_iKPlus[0] )->getP4(),     // K+
            parent->getDaug( m_iKMinus[0] )->getP4(),    // K-
            parent->getDaug( m_iPiPlus[0] )->getP4(),    // pi+
            parent->getDaug( m_iPiPlus[1] )->getP4(),    // pi+
            parent->getDaug( m_iPiMinus[0] )->getP4()    // pi-
        );
    } else if ( m_outCode == 10 ) {
        // 1=pi+ 2=pi+ 3=pi+ 4=pi+ 5=pi- 6=pi- 7=pi- with symmetrization of identical particles
        hardCur = m_WCurr->WCurrent_7pi(
            parent->getDaug( m_iPiPlus[0] )->getP4(),     // pi+
            parent->getDaug( m_iPiPlus[1] )->getP4(),     // pi+
            parent->getDaug( m_iPiPlus[2] )->getP4(),     // pi+
            parent->getDaug( m_iPiPlus[3] )->getP4(),     // pi+
            parent->getDaug( m_iPiMinus[0] )->getP4(),    // pi-
            parent->getDaug( m_iPiMinus[1] )->getP4(),    // pi-
            parent->getDaug( m_iPiMinus[2] )->getP4()     // pi-
        );
    } else if ( m_outCode == 11 ) {
        // 1=K+ 2 = pi+ 3 = pi+ 4 = pi- 5 = pi- with symmetrization
        hardCur = m_WCurr->WCurrent_K4pi(
            parent->getDaug( m_iKPlus[0] )->getP4(),      // K+
            parent->getDaug( m_iPiPlus[0] )->getP4(),     // pi+
            parent->getDaug( m_iPiPlus[1] )->getP4(),     // pi+
            parent->getDaug( m_iPiMinus[0] )->getP4(),    // pi-
            parent->getDaug( m_iPiMinus[1] )->getP4()     // pi-
        );
    }

    return hardCur;
}

//======================================================

void EvtBcVHad::decay( EvtParticle* parent )
{
    parent->initializePhaseSpace( getNDaug(), getDaugs() );

    // Calculate hadronic current
    const EvtVector4C hardCur = hardCurr( parent );

    EvtParticle* Jpsi = parent->getDaug( 0 );
    //std::cout<<"4th comp: "<<Jpsi->eps(3)<<std::endl;

    const EvtVector4R p4b( parent->mass(), 0., 0., 0. ),    // Bc momentum
        p4meson = Jpsi->getP4(),                            // J/psi momenta
        Q = p4b - p4meson, p4Sum = p4meson + p4b;
    const double Q2 = Q.mass2();

    // Calculate Bc -> V W form-factors
    double a1f( 0.0 ), a2f( 0.0 ), vf( 0.0 ), a0f( 0.0 );

    const double mMeson = Jpsi->mass();
    const double mB = parent->mass();
    const double mVar = mB + mMeson;

    m_FFModel->getvectorff( parent->getId(), Jpsi->getId(), Q2, mMeson, &a1f,
                            &a2f, &vf, &a0f );

    const double a3f = ( mVar / ( 2.0 * mMeson ) ) * a1f -
                       ( ( mB - mMeson ) / ( 2.0 * mMeson ) ) * a2f;

    // Calculate Bc -> V W current
    EvtTensor4C H = a1f * mVar * EvtTensor4C::g();
    H.addDirProd( ( -a2f / mVar ) * p4b, p4Sum );
    H += EvtComplex( 0.0, vf / mVar ) *
         dual( EvtGenFunctions::directProd( p4Sum, Q ) );
    H.addDirProd( ( a0f - a3f ) * 2.0 * ( mMeson / Q2 ) * p4b, Q );
    const EvtVector4C Heps = H.cont2( hardCur );

    for ( int i = 0; i < 3; i++ ) {
        const EvtVector4C eps =
            Jpsi->epsParent( i ).conj();    // psi-meson polarization vector
        const EvtComplex amp = eps * Heps;
        vertex( i, amp );
    }
}
