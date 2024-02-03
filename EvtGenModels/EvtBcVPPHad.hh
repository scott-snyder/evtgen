/***********************************************************************
* Copyright 1998-2023 CERN for the benefit of the EvtGen authors       *
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

#ifndef EvtBcVPPHad_HH
#define EvtBcVPPHad_HH

#include "EvtGenBase/EvtDecayAmp.hh"
#include "EvtGenBase/EvtVector4C.hh"

#include "EvtGenModels/EvtBCVFF2.hh"
#include "EvtGenModels/EvtWHad.hh"

#include <string>

class EvtParticle;

// Description: Module to implement Bc -> psi + p + pbar + pi decays

class EvtBcVPPHad : public EvtDecayAmp {
  public:
    std::string getName() override;
    EvtDecayBase* clone() override;
    void initProbMax() override;
    void init() override;
    void decay( EvtParticle* parent ) override;

  protected:
    // Hadronic current function
    EvtVector4C hardCurrPP( EvtParticle* parent, int i1, int i2 ) const;

  private:
    // Code of the Bc -> VW formfactor set:
    // 1 - SR
    // 2 - PM
    int m_whichFit;

    // Final vector particle code
    int m_idVector;

    // Code of the hadronic final state
    // 1: p+ p- pi+
    int m_outCode;

    std::unique_ptr<EvtBCVFF2> m_FFModel;
    std::unique_ptr<EvtWHad> m_WCurr;
};

#endif
