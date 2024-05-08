
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

#ifndef EVTFLATTE_HH
#define EVTFLATTE_HH

#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtVector4R.hh"

#include <vector>

using std::vector;

// Helper class

class EvtFlatteParam {
  public:
    EvtFlatteParam( double m1, double m2, double g ) :
        m_m1( m1 ), m_m2( m2 ), m_g( g )
    {
    }

    inline double m1() const { return m_m1; }
    inline double m2() const { return m_m2; }
    inline double g() const { return m_g; }

  private:
    double m_m1, m_m2, m_g;
};

//class declaration

class EvtFlatte final {
  public:
    //operator
    EvtFlatte& operator=( const EvtFlatte& );

    //constructor with all information about the resonance
    EvtFlatte( const EvtVector4R& p4_p, const EvtVector4R& p4_d1,
               const EvtVector4R& p4_d2, double ampl, double theta, double mass,
               vector<EvtFlatteParam>& params );

    //accessors
    //return 4-momenta of the particles involved
    inline const EvtVector4R& p4_p() { return m_p4_p; }
    inline const EvtVector4R& p4_d1() { return m_p4_d1; }
    inline const EvtVector4R& p4_d2() { return m_p4_d2; }

    //return amplitude
    inline double amplitude() { return m_ampl; }

    //return theta
    inline double theta() { return m_theta; }

    //return bwm
    inline double mass() { return m_mass; }

    //functions

    //calculate amplitude for this resonance
    EvtComplex resAmpl();

  private:
    inline EvtComplex sqrtCplx( double in )
    {
        return ( in > 0 ) ? EvtComplex( sqrt( in ), 0 )
                          : EvtComplex( 0, sqrt( -in ) );
    }

    EvtVector4R m_p4_p, m_p4_d1, m_p4_d2;
    double m_ampl, m_theta, m_mass;
    vector<EvtFlatteParam> m_params;
};

#endif
