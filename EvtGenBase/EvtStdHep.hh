
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

#ifndef EVTSTDHEP_HH
#define EVTSTDHEP_HH

#include "EvtGenBase/EvtVector4R.hh"

#include <iosfwd>

const int EVTSTDHEPLENGTH = 1000;

class EvtStdHep {
  public:
    EvtStdHep() {}
    ~EvtStdHep() {}

    void init();

    int getFirstMother( int i ) { return m_prntfirst[i]; }
    int getLastMother( int i ) { return m_prntlast[i]; }
    int getFirstDaughter( int i ) { return m_daugfirst[i]; }
    int getLastDaughter( int i ) { return m_dauglast[i]; }

    int getStdHepID( int i ) { return m_id[i]; }
    int getIStat( int i ) { return m_istat[i]; }

    EvtVector4R getP4( int i ) { return m_p4[i]; }
    EvtVector4R getX4( int i ) { return m_x[i]; }

    void translate( EvtVector4R d );

    int getNPart();
    void createParticle( EvtVector4R p4, EvtVector4R x, int prntfirst,
                         int prntlast, int id );

    friend std::ostream& operator<<( std::ostream& s, const EvtStdHep& stdhep );

  private:
    int m_npart;
    EvtVector4R m_p4[EVTSTDHEPLENGTH];
    EvtVector4R m_x[EVTSTDHEPLENGTH];
    int m_prntfirst[EVTSTDHEPLENGTH];
    int m_prntlast[EVTSTDHEPLENGTH];
    int m_daugfirst[EVTSTDHEPLENGTH];
    int m_dauglast[EVTSTDHEPLENGTH];
    int m_id[EVTSTDHEPLENGTH];
    int m_istat[EVTSTDHEPLENGTH];
};

#endif
