
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

#ifndef EVT_TWO_BODY_VERTEX_HH
#define EVT_TWO_BODY_VERTEX_HH

#include "EvtGenBase/EvtBlattWeisskopf.hh"
#include "EvtGenBase/EvtTwoBodyKine.hh"

#include <iostream>
#include <memory>

// Two-body propagator vertex AB->A,B with an attached Blatt-Weisskopf form factor.

class EvtTwoBodyVertex {
  public:
    EvtTwoBodyVertex();
    EvtTwoBodyVertex( double mA, double mB, double mAB, int L );
    EvtTwoBodyVertex( const EvtTwoBodyVertex& other );
    EvtTwoBodyVertex& operator=( const EvtTwoBodyVertex& other );

    double widthFactor( EvtTwoBodyKine x ) const;
    double formFactor( EvtTwoBodyKine x ) const;
    double phaseSpaceFactor( EvtTwoBodyKine x, EvtTwoBodyKine::Index ) const;

    inline int L() const { return m_LL; }
    inline double mA() const { return m_kine.mA(); }
    inline double mB() const { return m_kine.mB(); }
    inline double mAB() const { return m_kine.mAB(); }
    inline double pD() const { return m_p0; }
    void print( std::ostream& os ) const;

    void set_f( double R );

  private:
    EvtTwoBodyKine m_kine;
    int m_LL;
    double m_p0;
    std::unique_ptr<EvtBlattWeisskopf> m_f;    // optional Blatt-Weisskopf form factor
};

std::ostream& operator<<( std::ostream& os, const EvtTwoBodyVertex& v );

#endif
