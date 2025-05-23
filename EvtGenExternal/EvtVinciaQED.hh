
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

#ifdef EVTGEN_VINCIA
#ifndef EVTVINCIAQED_HH
#define EVTVINCIAQED_HH

#include "EvtGenBase/EvtAbsRadCorr.hh"
#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtVector4R.hh"

#include "EvtGenExternal/EvtPythiaRandom.hh"

#include "Pythia8/ParticleData.h"
#include "Pythia8/Pythia.h"
#include "Pythia8/ShowerModel.h"

#include <string>

class EvtParticle;

/*
 * Description: EvtGen's module for generation of
 *              QED final state radiation using VinciaQED.
 */

class EvtVinciaQED : public EvtAbsRadCorr {
  public:
    EvtVinciaQED( std::string xmlDir = "./xmldoc", bool useEvtGenRandom = true,
                  double infraredCutOff = 1.0e-7 );

    void initialise() override;

    void doRadCorr( EvtParticle* theParticle ) override;

  private:
    void storeDaughterInfo( Pythia8::Event& theEvent, EvtParticle* theParticle,
                            int startInt, int ancestorId );

    // Specify if we are going to use the random number generator (engine) from EvtGen for Pythia 8.
    bool m_useEvtGenRandom{ true };

    // Minimum photon Energy (infrared cut-off).
    double m_infraredCutOff{ 1.0e-7 };

    EvtId m_gammaId = EvtId( -1, -1 );
    long int m_gammaPDG{ 22 };

    std::unique_ptr<Pythia8::Pythia> m_vinciaPythiaGen;
    std::shared_ptr<EvtPythiaRandom> m_evtgenRandom;

    bool m_initialised{ false };
};

#endif

#endif
