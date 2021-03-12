
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

#include "EvtGenBase/EvtRadCorr.hh"

#include "EvtGenBase/EvtAbsRadCorr.hh"
#include "EvtGenBase/EvtReport.hh"

#include <cstdlib>
#include <iostream>

using std::endl;

thread_local EvtAbsRadCorr* EvtRadCorr::m_fsrEngine = nullptr;
bool EvtRadCorr::m_alwaysRadCorr = false;
bool EvtRadCorr::m_neverRadCorr = false;

void EvtRadCorr::setRadCorrEngine( EvtAbsRadCorr* fsrEngine )
{
    m_fsrEngine = fsrEngine;

    if ( m_fsrEngine ) {
        m_fsrEngine->initialise();
    }
}

void EvtRadCorr::doRadCorr( EvtParticle* p )
{
    if ( m_fsrEngine == nullptr ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "No RadCorr model available in "
            << "EvtRadCorr::doRadCorr()." << endl;
        ::abort();
    }

    if ( !m_neverRadCorr )
        m_fsrEngine->doRadCorr( p );
    return;
}

bool EvtRadCorr::alwaysRadCorr()
{
    return m_alwaysRadCorr;
}
bool EvtRadCorr::neverRadCorr()
{
    return m_neverRadCorr;
}

void EvtRadCorr::setAlwaysRadCorr()
{
    m_alwaysRadCorr = true;
    m_neverRadCorr = false;
}
void EvtRadCorr::setNeverRadCorr()
{
    m_alwaysRadCorr = false;
    m_neverRadCorr = true;
}
void EvtRadCorr::setNormalRadCorr()
{
    m_alwaysRadCorr = false;
    m_neverRadCorr = false;
}
