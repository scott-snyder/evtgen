
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

#ifndef EVTPARTICLEDECAY_HH
#define EVTPARTICLEDECAY_HH

#include "EvtGenBase/EvtDecayBase.hh"

class EvtParticleDecay {
  public:
    EvtParticleDecay()
    {
        m_decay = nullptr;
        m_brfrsum = 0.0;
        m_massmin = 0.0;
    }

    ~EvtParticleDecay()
    {
        if ( m_decay != nullptr )
            delete m_decay;
    }

    void chargeConj( EvtParticleDecay* decay );

    void setDecayModel( EvtDecayBase* decay ) { m_decay = decay; }
    EvtDecayBase* getDecayModel() { return m_decay; }
    double getBrfrSum() { return m_brfrsum; }
    void setBrfrSum( double brfrsum ) { m_brfrsum = brfrsum; }
    double getMassMin() { return m_massmin; }
    void setMassMin( double massmin ) { m_massmin = massmin; }

    void printSummary();

  private:
    EvtDecayBase* m_decay;

    double m_brfrsum;
    double m_massmin;
};

#endif
