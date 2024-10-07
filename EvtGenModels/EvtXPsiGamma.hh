
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

#ifndef EVTXPSIGAMMA_HH
#define EVTXPSIGAMMA_HH

#include "EvtGenBase/EvtDecayAmp.hh"
#include "EvtGenBase/EvtId.hh"

#include <fstream>
#include <stdio.h>

class EvtComplex;
class EvtParticle;
class EvtTensor4C;
class EvtVector4C;
class EvtVector4R;

// Description:Implementation of the X3872(2-+) -> J/psi gamma decay
// Description: Routine to implement radiative decay X3872(2-+) -> J/psi gamma
//      according to F. Brazzi et al, arXiv:1103.3155

class EvtXPsiGamma : public EvtDecayAmp {
  public:
    std::string getName() const override;
    EvtDecayBase* clone() const override;

    void init() override;
    void initProbMax() override;
    void decay( EvtParticle* p ) override;

  protected:
    // This function is not declared const because epsParentPhoton is not const
    void calcAmp( EvtParticle& parent, EvtAmp& amp );

  private:
    double calcPstar( double m_parent, double m_1, double m_2 ) const;

    EvtComplex fT2( EvtVector4R p, EvtVector4R q, EvtTensor4C epsPI,
                    EvtVector4C epsEps, EvtVector4C epsEta ) const;
    EvtComplex fT3( EvtVector4R p, EvtVector4R q, EvtTensor4C epsPI,
                    EvtVector4C epsEps, EvtVector4C epsEta ) const;

    EvtId m_ID0;

    double m_gOmega;
    double m_gPOmega;
    double m_gRho;
    double m_gPRho;
    double m_fOmega;
    double m_fRho;
};

#endif
