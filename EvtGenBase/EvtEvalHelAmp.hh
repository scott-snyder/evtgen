
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

#ifndef EVTEVALHELAMP_HH
#define EVTEVALHELAMP_HH

#include "EvtGenBase/EvtDecayAmp.hh"
#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtSpinType.hh"

class EvtParticle;
class EvtAmp;

class EvtEvalHelAmp {
  public:
    EvtEvalHelAmp( EvtId idA, EvtId idB, EvtId idC, EvtComplexPtrPtr HBC );

    virtual ~EvtEvalHelAmp();

    double probMax();

    void evalAmp( EvtParticle* p, EvtAmp& amp );

  private:
    void fillHelicity( int* lambda2, int n, int J2, EvtId id );
    void setUpRotationMatrices( EvtParticle* p, double theta, double phi );
    void applyRotationMatrices();

    //spins states available for particle A, B, and C.
    int m_nA, m_nB, m_nC;

    //helicity amplitudes
    EvtComplexPtrPtr m_HBC;

    //2 times spin for each of the particles
    int m_JA2, m_JB2, m_JC2;

    //2 times the helicity for the states
    int *m_lambdaA2, *m_lambdaB2, *m_lambdaC2;

    //Rotation matrices
    EvtComplexPtrPtr m_RA, m_RB, m_RC;

    //temporary array for amplitudes
    EvtComplexPtrPtrPtr m_amp, m_amp1, m_amp3;
};

#endif
