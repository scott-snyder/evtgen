
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

#include "EvtGenModels/EvtBcSMuNu.hh"

#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtSemiLeptonicScalarAmp.hh"

#include "EvtGenModels/EvtBCSFF.hh"

#include <iostream>
#include <stdlib.h>
#include <string>

using namespace std;

std::string EvtBcSMuNu::getName() const
{
    return "BC_SMN";
}

EvtDecayBase* EvtBcSMuNu::clone() const
{
    return new EvtBcSMuNu;
}

void EvtBcSMuNu::decay( EvtParticle* p )
{
    p->initializePhaseSpace( getNDaug(), getDaugs() );
    m_calcamp->CalcAmp( p, m_amp2, m_ffmodel.get() );
}

void EvtBcSMuNu::init()
{
    checkNArg( 1 );
    checkNDaug( 3 );

    //We expect the parent to be a scalar
    //and the daughters to be X lepton neutrino

    checkSpinParent( EvtSpinType::SCALAR );

    checkSpinDaughter( 0, EvtSpinType::SCALAR );
    checkSpinDaughter( 1, EvtSpinType::DIRAC );
    checkSpinDaughter( 2, EvtSpinType::NEUTRINO );

    m_idScalar = getDaug( 0 ).getId();
    m_whichfit = int( getArg( 0 ) + 0.1 );

    m_ffmodel = std::make_unique<EvtBCSFF>( m_idScalar, m_whichfit );

    m_calcamp = std::make_unique<EvtSemiLeptonicScalarAmp>();
}

void EvtBcSMuNu::initProbMax()
{
    EvtId parId = getParentId();
    EvtId mesonId = getDaug( 0 );
    EvtId lepId = getDaug( 1 );
    EvtId nuId = getDaug( 2 );

    int nQ2Bins = 200;
    double maxProb = m_calcamp->CalcMaxProb( parId, mesonId, lepId, nuId,
                                             m_ffmodel.get(), nQ2Bins );

    if ( verbose() ) {
        EvtGenReport( EVTGEN_INFO, "EvtBcSMuNu" )
            << "Max prob = " << maxProb << endl;
    }

    setProbMax( maxProb );
}
