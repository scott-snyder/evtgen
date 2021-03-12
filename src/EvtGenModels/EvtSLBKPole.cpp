
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

#include "EvtGenModels/EvtSLBKPole.hh"    //modified

#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtSemiLeptonicScalarAmp.hh"
#include "EvtGenBase/EvtSemiLeptonicTensorAmp.hh"
#include "EvtGenBase/EvtSemiLeptonicVectorAmp.hh"

#include "EvtGenModels/EvtSLBKPoleFF.hh"    //modified

#include <string>

std::string EvtSLBKPole::getName() const
{
    return "SLBKPOLE";    //modified
}

EvtDecayBase* EvtSLBKPole::clone() const
{    //modified

    return new EvtSLBKPole;
}

void EvtSLBKPole::decay( EvtParticle* p )
{    //modified

    p->initializePhaseSpace( getNDaug(), getDaugs() );

    m_calcamp->CalcAmp( p, m_amp2, m_SLBKPoleffmodel.get() );    //modified
    return;
}

void EvtSLBKPole::initProbMax()
{
    EvtId parnum, mesnum, lnum, nunum;

    parnum = getParentId();
    mesnum = getDaug( 0 );
    lnum = getDaug( 1 );
    nunum = getDaug( 2 );

    double mymaxprob = m_calcamp->CalcMaxProb( parnum, mesnum, lnum, nunum,
                                               m_SLBKPoleffmodel.get() );    //modified

    setProbMax( mymaxprob );
}

void EvtSLBKPole::init()
{    //modified

    checkNDaug( 3 );

    //We expect the parent to be a scalar
    //and the daughters to be X lepton neutrino

    checkSpinParent( EvtSpinType::SCALAR );
    checkSpinDaughter( 1, EvtSpinType::DIRAC );
    checkSpinDaughter( 2, EvtSpinType::NEUTRINO );

    EvtSpinType::spintype mesontype = EvtPDL::getSpinType( getDaug( 0 ) );

    m_SLBKPoleffmodel = std::make_unique<EvtSLBKPoleFF>( getNArg(),
                                                         getArgs() );    //modified

    switch ( mesontype ) {
        case EvtSpinType::SCALAR:
            m_calcamp = std::make_unique<EvtSemiLeptonicScalarAmp>();
            break;
        case EvtSpinType::VECTOR:
            m_calcamp = std::make_unique<EvtSemiLeptonicVectorAmp>();
            break;
        case EvtSpinType::TENSOR:
            m_calcamp = std::make_unique<EvtSemiLeptonicTensorAmp>();
            break;
        default:;
    }
}
