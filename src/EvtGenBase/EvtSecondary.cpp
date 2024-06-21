
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

#include "EvtGenBase/EvtSecondary.hh"

#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtReport.hh"

#include <iostream>
using std::endl;
using std::ostream;

void EvtSecondary::init()
{
    m_npart = 0;
}

int EvtSecondary::getNPart()
{
    return m_npart;
}

void EvtSecondary::createSecondary( int stdhepindex, EvtParticle* prnt )
{
    m_stdhepindex[m_npart] = stdhepindex;
    if ( prnt->getNDaug() == 0 ) {
        m_id1[m_npart] = 0;
        m_id2[m_npart] = 0;
        m_id3[m_npart] = 0;
        m_npart++;
        return;
    }
    if ( prnt->getNDaug() == 1 ) {
        m_id1[m_npart] = EvtPDL::getStdHep( prnt->getDaug( 0 )->getId() );
        m_id2[m_npart] = 0;
        m_id3[m_npart] = 0;
        m_npart++;
        return;
    }
    if ( prnt->getNDaug() == 2 ) {
        m_id1[m_npart] = EvtPDL::getStdHep( prnt->getDaug( 0 )->getId() );
        m_id2[m_npart] = EvtPDL::getStdHep( prnt->getDaug( 1 )->getId() );
        m_id3[m_npart] = 0;
        m_npart++;
        return;
    }
    if ( prnt->getNDaug() == 3 ) {
        m_id1[m_npart] = EvtPDL::getStdHep( prnt->getDaug( 0 )->getId() );
        m_id2[m_npart] = EvtPDL::getStdHep( prnt->getDaug( 1 )->getId() );
        m_id3[m_npart] = EvtPDL::getStdHep( prnt->getDaug( 2 )->getId() );
        m_npart++;
        return;
    }

    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "More than 3 decay products in a secondary particle!" << endl;
}

ostream& operator<<( ostream& s, const EvtSecondary& secondary )
{
    s << endl;
    s << "Secondary decays:" << endl;

    int i;
    for ( i = 0; i < secondary.m_npart; i++ ) {
        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << i << " " << secondary.m_stdhepindex[i] << " "
            << secondary.m_id1[i] << " " << secondary.m_id2[i] << " "
            << secondary.m_id3[i] << endl;
    }

    s << endl;

    return s;
}
