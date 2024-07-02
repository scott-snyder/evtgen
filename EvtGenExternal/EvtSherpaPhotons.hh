
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

#ifdef EVTGEN_SHERPA
#ifndef EVTSHERPAPHOTONS_HH
#define EVTSHERPAPHOTONS_HH

#include "EvtGenBase/EvtAbsRadCorr.hh"
#include "EvtGenBase/EvtHepMCEvent.hh"
#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"

#include "SHERPA/Main/Sherpa.H"

#include <string>

class EvtParticle;
class EvtAbsExternalGen;

/* Description: EvtGen's interface to PHOTOS for generation of
 *              QED final state radiation.
 */
class EvtSherpaPhotons : public EvtAbsRadCorr {
  public:
    EvtSherpaPhotons( const bool useEvtGenRandom = true,
                      const double infraredCutOff = 1.0e-7, const int mode = 2,
                      const int useME = 0 );

    void initialise() override;

    void doRadCorr( EvtParticle* p ) override;

  private:
    // Provides vector of pointers to the configuration strings for Sherpa
    std::vector<char*> addParameters();

    // Updates the particle properties table of Sherpa
    void updateParticleLists();

    // Vector containing the configuration strings for Sherpa
    // INIT_ONLY=6 intialises the Sherpa objects without launching simulation.
    std::vector<std::string> m_configs{ "Sherpa", "INIT_ONLY=6" };

    // Use EvtGen's random number generator
    bool m_useEvtGenRandom = true;

    // Default settings for PHOTONS++
    // Minimum photon Energy (infrared cut-off).
    double m_infraredCutOff = 1.0e-7;
    // Switches on (default) the hard emission corrections (mode 2),
    // otherwise follow the soft-photon only approach (mode 1).
    int m_mode = 2;
    // Switches off (default) or on the exact matrix-element corrections (works only for mode 2)
    int m_useME = 0;

    // Default photon properties
    const std::string m_photonType = "gamma";
    EvtId m_gammaId = EvtId( -1, -1 );
    long int m_gammaPDG = 22;

    // The Sherpa instance.
    static std::unique_ptr<SHERPA::Sherpa> m_sherpaGen;
    static bool m_initialised;
    static std::mutex m_sherpa_mutex;
};

#endif

#endif
