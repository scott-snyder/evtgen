
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

#ifndef EVTDALITZDECAYINFO_HH
#define EVTDALITZDECAYINFO_HH

#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtDalitzReso.hh"
#include "EvtGenBase/EvtId.hh"

#include <vector>

// Description: Model to describe a generic dalitz decay

class EvtDalitzDecayInfo final {
  public:
    EvtDalitzDecayInfo( EvtId d1, EvtId d2, EvtId d3 ) :
        m_d1( d1 ), m_d2( d2 ), m_d3( d3 ), m_probMax( 0. )
    {
    }

    void addResonance( EvtComplex amp, EvtDalitzReso res )
    {
        m_resonances.push_back( std::pair<EvtComplex, EvtDalitzReso>( amp, res ) );
    }
    void addResonance( std::pair<EvtComplex, EvtDalitzReso> res )
    {
        m_resonances.push_back( res );
    }
    void setProbMax( double probMax ) { m_probMax = probMax; }

    const std::vector<std::pair<EvtComplex, EvtDalitzReso>>& getResonances() const
    {
        return m_resonances;
    }
    double getProbMax() const { return m_probMax; }

    inline const EvtId& daughter1() const { return m_d1; }
    inline const EvtId& daughter2() const { return m_d2; }
    inline const EvtId& daughter3() const { return m_d3; }

  private:
    EvtId m_d1, m_d2, m_d3;
    std::vector<std::pair<EvtComplex, EvtDalitzReso>> m_resonances;
    double m_probMax;
};

#endif
