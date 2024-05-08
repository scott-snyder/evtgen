
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

#ifndef EVTSECONDARY_HH
#define EVTSECONDARY_HH

const int EVTSECONDARYLENGTH = 100;

class EvtParticle;
#include <iosfwd>

class EvtSecondary {
  public:
    EvtSecondary() {}
    ~EvtSecondary() {}

    void init();

    int getStdHepIndex( int i ) { return m_stdhepindex[i]; }
    int getD1( int i ) { return m_id1[i]; }
    int getD2( int i ) { return m_id2[i]; }
    int getD3( int i ) { return m_id3[i]; }

    int getNPart();
    void createSecondary( int stdhepindex, EvtParticle* prnt );

    friend std::ostream& operator<<( std::ostream& s,
                                     const EvtSecondary& secondary );

  private:
    int m_npart;
    int m_stdhepindex[EVTSECONDARYLENGTH];
    int m_id1[EVTSECONDARYLENGTH];
    int m_id2[EVTSECONDARYLENGTH];
    int m_id3[EVTSECONDARYLENGTH];
};

#endif
