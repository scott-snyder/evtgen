
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

#include "EvtGenExternal/EvtExternalGenFactory.hh"

#include "EvtGenBase/EvtReport.hh"

#ifdef EVTGEN_PYTHIA
#include "EvtGenExternal/EvtPythiaEngine.hh"
#endif

#ifdef EVTGEN_TAUOLA
#include "EvtGenExternal/EvtTauolaEngine.hh"
#endif

#include <iostream>
using std::endl;

EvtExternalGenFactory& EvtExternalGenFactory::getInstance()
{
    static thread_local EvtExternalGenFactory theFactory;

    return theFactory;
}

// Only define the generator if we have the external ifdef variable set
#ifdef EVTGEN_PYTHIA
void EvtExternalGenFactory::definePythiaGenerator( std::string xmlDir,
                                                   bool convertPhysCodes,
                                                   bool useEvtGenRandom )
{
    EvtGenReport( EVTGEN_INFO, "EvtGen" )
        << "Defining EvtPythiaEngine: data tables defined in " << xmlDir << endl;
    if ( convertPhysCodes == true ) {
        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << "Pythia 6 codes in decay files will be converted to Pythia 8 codes"
            << endl;
    } else {
        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << "Pythia 8 codes need to be used in decay files" << endl;
    }

    if ( useEvtGenRandom == true ) {
        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << "Using EvtGen random engine for Pythia 8 as well" << endl;
    }

    m_extGenMap[GenId::PythiaGenId] = std::make_unique<EvtPythiaEngine>(
        xmlDir, convertPhysCodes, useEvtGenRandom );
}
#else
void EvtExternalGenFactory::definePythiaGenerator( std::string, bool, bool )
{
}
#endif

#ifdef EVTGEN_TAUOLA
void EvtExternalGenFactory::defineTauolaGenerator( bool useEvtGenRandom )
{
    EvtGenReport( EVTGEN_INFO, "EvtGen" ) << "Defining EvtTauolaEngine." << endl;

    m_extGenMap[GenId::TauolaGenId] = std::make_unique<EvtTauolaEngine>(
        useEvtGenRandom );
}
#else
void EvtExternalGenFactory::defineTauolaGenerator( bool )
{
}
#endif

EvtAbsExternalGen* EvtExternalGenFactory::getGenerator( const GenId genId )
{
    ExtGenMap::iterator iter = m_extGenMap.find( genId );

    if ( iter == m_extGenMap.end() ) {
        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << "EvtAbsExternalGen::getGenerator: could not find generator for genId = "
            //FIXME C++23 use std::to_underlying
            << static_cast<std::underlying_type_t<GenId>>( genId ) << endl;
        return nullptr;
    }

    // Retrieve the external generator engine
    auto& theGenerator = iter->second;
    return theGenerator.get();
}

void EvtExternalGenFactory::initialiseAllGenerators()
{
    for ( auto& [id, theGenerator] : m_extGenMap ) {
        if ( theGenerator ) {
            theGenerator->initialise();
        }
    }
}
