
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

#ifdef EVTGEN_PHOTOS
#ifndef EVTPHOTOS_HH
#define EVTPHOTOS_HH

#include "EvtGenBase/EvtAbsRadCorr.hh"
#include "EvtGenBase/EvtHepMCEvent.hh"
#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtVector4R.hh"

#ifdef EVTGEN_HEPMC3
#include "HepMC3/Units.h"

#include "Photos/PhotosHepMC3Event.h"
#include "Photos/PhotosHepMC3Particle.h"
#else
#include "Photos/PhotosHepMCEvent.h"
#include "Photos/PhotosHepMCParticle.h"
#include "Photos/PhotosParticle.h"
#endif

#include <mutex>
#include <string>

class EvtParticle;
class EvtAbsExternalGen;

/*
 * Description: EvtGen's interface to PHOTOS for generation of
 *              QED final-state radiation.
 */

class EvtPHOTOS : public EvtAbsRadCorr {
  public:
    EvtPHOTOS( const std::string& photonType = "gamma",
               const bool useEvtGenRandom = true,
               const double infraredCutOff = 1.0e-7,
               const double maxWtInterference = 64.0 );

    void initialise() override;

    void doRadCorr( EvtParticle* theParticle ) override;

  private:
    GenParticlePtr createGenParticle( const EvtParticle& theParticle,
                                      bool incoming ) const;

    int getNumberOfPhotons( const GenVertexPtr theVertex ) const;

    // Use EvtGen's random number generator
    bool m_useEvtGenRandom = true;

    // Default settings for PHOTOS
    // Minimum photon Energy (infrared cut-off).
    double m_infraredCutOff = 1.0e-7;
    // Maximum value of the interference weight
    double m_maxWtInterference = 64.0;

    // Default photon type, id, pdg number and mass
    std::string m_photonType = "gamma";
    EvtId m_gammaId = EvtId( -1, -1 );
    long int m_gammaPDG = 22;

    static bool m_initialised;
    static std::mutex m_photos_mutex;
};

#endif

#endif
