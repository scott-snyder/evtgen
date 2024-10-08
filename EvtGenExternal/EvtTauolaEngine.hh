
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

#include <mutex>

// Description: Interface to the TAUOLA external generator

class EvtTauolaEngine : public EvtAbsExternalGen {
  public:
    EvtTauolaEngine( bool useEvtGenRandom = true, bool seedTauolaFortran = true );

    bool doDecay( EvtParticle* theMother ) override;

    void initialise() override;

  private:
    GenParticlePtr createGenParticle( const EvtParticle* theParticle ) const;

    void setUpPossibleTauModes();
    void setOtherParameters();

    int getModeInt( EvtDecayBase* decayModel ) const;

    void decayTauEvent( EvtParticle* tauParticle );

    bool m_useEvtGenRandom{ true };
    bool m_seedTauolaFortran{ true };

    // PDG standard code integer ID for tau particle
    static constexpr int m_tauPDG{ 15 };
    // Number of possible decay modes in Tauola
    static constexpr int m_nTauolaModes{ 22 };
    // Neutral and charged spin propagator choices
    static int m_neutPropType;
    static int m_posPropType;
    static int m_negPropType;

    static bool m_initialised;
    static std::mutex m_tauola_mutex;
};

#endif

#endif
