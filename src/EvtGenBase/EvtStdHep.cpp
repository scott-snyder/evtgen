
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

#include "EvtGenBase/EvtStdHep.hh"

#include "EvtGenBase/EvtVector4R.hh"

#include <iomanip>
#include <iostream>
using namespace std;

void EvtStdHep::init()
{
    m_npart = 0;
}

int EvtStdHep::getNPart()
{
    return m_npart;
}

void EvtStdHep::createParticle( EvtVector4R p4, EvtVector4R x, int prntfirst,
                                int prntlast, int id )
{
    m_p4[m_npart] = p4;
    m_x[m_npart] = x;
    m_prntfirst[m_npart] = prntfirst;
    m_prntlast[m_npart] = prntlast;
    m_daugfirst[m_npart] = -1;
    m_dauglast[m_npart] = -1;
    m_id[m_npart] = id;
    m_istat[m_npart] = 1;

    //we also need to fix up the parents pointer to the daughter!

    if ( prntfirst >= 0 ) {
        int i;
        for ( i = prntfirst; i <= prntlast; i++ ) {
            m_istat[i] = 2;
            if ( m_daugfirst[i] == -1 )
                m_daugfirst[i] = m_npart;
            if ( m_dauglast[i] < m_npart )
                m_dauglast[i] = m_npart;
        }
    }

    m_npart++;
}

void EvtStdHep::translate( EvtVector4R d )
{
    int i;
    for ( i = 0; i < m_npart; i++ ) {
        m_x[i] += d;
    }
}

ostream& operator<<( ostream& s, const EvtStdHep& stdhep )
{
    int w = s.width();
    int p = s.precision();
    std::ios::fmtflags f = s.flags();

    s << endl;
    s << "  N      Id Ist   M1   M2   DF   DL      px      py      pz       E       t       x       y       z"
      << endl;
    int i;
    for ( i = 0; i < stdhep.m_npart; i++ ) {
        s.width( 3 );
        s << i << " ";
        s.width( 7 );
        s << stdhep.m_id[i] << " ";
        s.width( 3 );
        s << stdhep.m_istat[i] << " ";
        s.width( 4 );
        s << stdhep.m_prntfirst[i] << " ";
        s.width( 4 );
        s << stdhep.m_prntlast[i] << " ";
        s.width( 4 );
        s << stdhep.m_daugfirst[i] << " ";
        s.width( 4 );
        s << stdhep.m_dauglast[i] << " ";
        s.width( 7 );
        s.precision( 4 );
        s << setiosflags( ios::right | ios::fixed );
        s << stdhep.m_p4[i].get( 1 ) << " ";
        s.width( 7 );
        s.precision( 4 );
        s << setiosflags( ios::right | ios::fixed );
        s << stdhep.m_p4[i].get( 2 ) << " ";
        s.width( 7 );
        s.precision( 4 );
        s << setiosflags( ios::right | ios::fixed );
        s << stdhep.m_p4[i].get( 3 ) << " ";
        s.width( 7 );
        s.precision( 4 );
        s << setiosflags( ios::right | ios::fixed );
        s << stdhep.m_p4[i].get( 0 ) << " ";
        s.width( 7 );
        s.precision( 4 );
        s << setiosflags( ios::right | ios::fixed );
        s << stdhep.m_x[i].get( 0 ) << " ";
        s.width( 7 );
        s.precision( 4 );
        s << setiosflags( ios::right | ios::fixed );
        s << stdhep.m_x[i].get( 1 ) << " ";
        s.width( 7 );
        s.precision( 4 );
        s << setiosflags( ios::right | ios::fixed );
        s << stdhep.m_x[i].get( 2 ) << " ";
        s.width( 7 );
        s.precision( 4 );
        s << setiosflags( ios::right | ios::fixed );
        s << stdhep.m_x[i].get( 3 ) << endl;
        s.width( 0 );
    }

    s << endl;

    s.width( w );
    s.precision( p );
    s.flags( (std::ios::fmtflags)f );

    return s;
}
