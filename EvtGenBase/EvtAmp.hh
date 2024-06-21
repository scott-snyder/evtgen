
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

#ifndef EVTAMP_HH
#define EVTAMP_HH

#include "EvtGenBase/EvtComplex.hh"

class EvtSpinDensity;
class EvtId;

class EvtAmp {
  public:
    EvtAmp();
    EvtAmp( const EvtAmp& amp );

    void init( EvtId p, int ndaug, EvtId* daug );

    void setAmp( int* ind, const EvtComplex& amp );

    const EvtComplex& getAmp( int* ind ) const;

    EvtSpinDensity getSpinDensity();

    EvtSpinDensity contract( int i, const EvtAmp& a );
    EvtAmp contract( int i, const EvtSpinDensity& rho );

    //sum over the i:th daugther of a1 and contract with parent of a2
    EvtAmp contract( int i, const EvtAmp& a1, const EvtAmp& a2 );

    EvtSpinDensity getForwardSpinDensity( EvtSpinDensity* rho_list, int k );
    EvtSpinDensity getBackwardSpinDensity( EvtSpinDensity* rho_list );

    EvtAmp& operator=( const EvtAmp& amp );

    /**
  * sets the amplitudes calculated in the decay objects
  */
    void vertex( const EvtComplex& amp );

    /**
  * sets the amplitudes calculated in the decay objects
  */
    void vertex( int i1, const EvtComplex& amp );

    /**
  * sets the amplitudes calculated in the decay objects
  */
    void vertex( int i1, int i2, const EvtComplex& amp );

    /**
  * sets the amplitudes calculated in the decay objects
  */
    void vertex( int i1, int i2, int i3, const EvtComplex& amp );

    /**
  * sets the amplitudes calculated in the decay objects
  */
    void vertex( int* i1, const EvtComplex& amp );

    void dump();

  private:
    friend class EvtDecayAmp;

    void setNDaug( int n );
    void setNState( int parent_states, int* daug_states );

    // the amplitudes
    EvtComplex m_amp[125];

    // the number of daughters
    int m_ndaug;

    // the number of states of the parent
    int m_pstates;

    // number of states of the daughter
    int m_dstates[10];

    // the nontrivial index of the daughter
    int m_dnontrivial[10];

    // number of nontrivial daugts+parent
    int m_nontrivial;

    // compact nstates
    int m_nstate[5];
};

#endif
