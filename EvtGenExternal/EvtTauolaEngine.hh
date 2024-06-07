
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

#ifdef EVTGEN_TAUOLA

#ifndef EVTTAUOLAENGINE_HH
#define EVTTAUOLAENGINE_HH

#include "EvtGenBase/EvtDecayBase.hh"
#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtVector4R.hh"

#include "EvtGenModels/EvtAbsExternalGen.hh"

#ifdef EVTGEN_HEPMC3
#include "HepMC3/Relatives.h"
#include "HepMC3/Units.h"

#include "Tauola/TauolaHepMC3Event.h"
#include "Tauola/TauolaHepMC3Particle.h"
#else
#include "Tauola/TauolaHepMCEvent.h"
#include "Tauola/TauolaHepMCParticle.h"
#include "Tauola/TauolaParticle.h"
#endif
#include "EvtGenBase/EvtHepMCEvent.hh"

#include <map>
#include <vector>

// Description: Interface to the TAUOLA external generator

class EvtTauolaEngine : public EvtAbsExternalGen {
  public:
    EvtTauolaEngine( bool useEvtGenRandom = true );

    bool doDecay( EvtParticle* theMother ) override;

    void initialise() override;

  protected:
  private:
    GenParticlePtr createGenParticle( EvtParticle* theParticle );

    void setUpPossibleTauModes();
    void setOtherParameters();

    int getModeInt( EvtDecayBase* decayModel );

    void decayTauEvent( EvtParticle* tauParticle );

    bool m_initialised{ false };
    // PDG standard code integer ID for tau particle
    int m_tauPDG{ 15 };
    // Number of possible decay modes in Tauola
    int m_nTauolaModes{ 22 };
    // Neutral and charged spin propagator choices
    int m_neutPropType{ 0 };
    int m_posPropType{ 0 };
    int m_negPropType{ 0 };
};

#endif

#endif
