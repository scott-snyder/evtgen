
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

#include "EvtGenExternal/EvtExternalGenList.hh"

#include "EvtGenModels/EvtNoRadCorr.hh"

#include "EvtGenExternal/EvtExternalGenFactory.hh"
#include "EvtGenExternal/EvtPHOTOS.hh"
#include "EvtGenExternal/EvtPythia.hh"
#include "EvtGenExternal/EvtSherpaPhotons.hh"
#include "EvtGenExternal/EvtTauola.hh"

EvtExternalGenList::EvtExternalGenList( bool convertPythiaCodes,
                                        std::string pythiaXmlDir,
                                        std::string photonType,
                                        bool useEvtGenRandom ) :
    m_photonType{ photonType }, m_useEvtGenRandom{ useEvtGenRandom }
{
    // Instantiate the external generator factory
    EvtExternalGenFactory& extFactory = EvtExternalGenFactory::getInstance();

    if ( pythiaXmlDir.size() < 1 ) {
        // If we have no string defined, check the value of the
        // PYTHIA8DATA environment variable which should be set to the
        // xmldoc Pythia directory
        char* pythiaDataDir = getenv( "PYTHIA8DATA" );
        if ( pythiaDataDir != nullptr ) {
            pythiaXmlDir = pythiaDataDir;
        }
    }

    extFactory.definePythiaGenerator( pythiaXmlDir, convertPythiaCodes,
                                      useEvtGenRandom );

    extFactory.defineTauolaGenerator( useEvtGenRandom );
}

EvtExternalGenList::~EvtExternalGenList()
{
}

#ifdef EVTGEN_PHOTOS
EvtAbsRadCorr* EvtExternalGenList::getPhotosModel( const double infraredCutOff,
                                                   const double maxWtInterference )
{
    // Define the Photos model, which uses the EvtPHOTOS class.
    EvtPHOTOS* photosModel = new EvtPHOTOS( m_photonType, m_useEvtGenRandom,
                                            infraredCutOff, maxWtInterference );
    return photosModel;
}
#else
EvtAbsRadCorr* EvtExternalGenList::getPhotosModel(
    const double /*infraredCutOff*/, const double /*maxWtInterference*/ )
{
    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << " PHOTOS generator has been called for FSR simulation, but it was not switched on during compilation."
        << std::endl;

    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << " The simulation will be generated without FSR." << std::endl;

    return new EvtNoRadCorr{};
}
#endif

#ifdef EVTGEN_SHERPA
EvtAbsRadCorr* EvtExternalGenList::getSherpaPhotonsModel(
    const double infraredCutOff, const int mode, const int useME )
{
    // Define the Photos model, which uses the EvtSherpaPhotonsEngine class.
    EvtSherpaPhotons* sherpaPhotonsModel =
        new EvtSherpaPhotons( m_useEvtGenRandom, infraredCutOff, mode, useME );
    return sherpaPhotonsModel;
}
#else
EvtAbsRadCorr* EvtExternalGenList::getSherpaPhotonsModel(
    const double /*infraredCutOff*/, const int /*mode*/, const int /*useME*/ )
{
    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << " Sherpa's PHOTONS++ generator has been called for FSR simulation, but Sherpa was not switched on during compilation."
        << std::endl;

    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << " The simulation will be generated without FSR." << std::endl;

    return new EvtNoRadCorr{};
}
#endif

std::list<EvtDecayBase*> EvtExternalGenList::getListOfModels()
{
    // Create the Pythia and Tauola models, which use their own engine classes.
    EvtPythia* pythiaModel = new EvtPythia();
    EvtTauola* tauolaModel = new EvtTauola();

    std::list<EvtDecayBase*> extraModels;
    extraModels.push_back( pythiaModel );
    extraModels.push_back( tauolaModel );

    // Return the list of models
    return extraModels;
}
