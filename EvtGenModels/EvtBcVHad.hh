
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

#ifndef EvtBcVHad_HH
#define EvtBcVHad_HH

#include "EvtGenBase/EvtDecayAmp.hh"
#include "EvtGenBase/EvtVector4C.hh"

#include "EvtGenModels/EvtBCVFF2.hh"
#include "EvtGenModels/EvtWHad.hh"

#include <array>
#include <memory>
#include <string>

class EvtParticle;

// Description: Module to implement Bc -> psi + (n pi) + (m K) decays

class EvtBcVHad : public EvtDecayAmp {
  public:
    std::string getName() const override;
    EvtDecayBase* clone() const override;
    void initProbMax() override;
    void init() override;
    void decay( EvtParticle* parent ) override;

  protected:
    // Hadronic current function
    EvtVector4C hardCurr( EvtParticle* parent ) const;
    void parseDecay();

  private:
    // Code of the Bc -> VW formfactor set:
    // 1 - SR
    // 2 - PM
    int m_whichFit;

    // Final vector particle code
    int m_idVector;

    // Code of the hadronic final state
    // 1:  B_c+ -> V pi+
    // 2:  B_c+ -> V pi+ pi0
    // 3:  B_c+ -> V 2pi+ pi-
    // 4:  B_c+ -> V 2pi+ pi- pi0 (not implemented)
    // 5:  B_c+ -> V 3pi+ 2pi-
    // 6:  B_c+ -> V K+ K- pi+
    // 7:  B_c+ -> V K+ pi+ pi-
    // 8:  B_c+ -> V K_S0 K+
    // 9:  B_c+ -> V K+ K- 2pi+ pi-
    // 10: B_c+ -> V 4pi+ 3pi-
    // 11: B_c+ -> V K+ 2pi+ 2pi-
    int m_outCode;

    std::unique_ptr<EvtBCVFF2> m_FFModel;
    std::unique_ptr<EvtWHad> m_WCurr;

    std::array<int, 4> m_iPiPlus = {
        {-1, -1, -1, -1}
    };
    std::array<int, 4> m_iPiMinus = {
        {-1, -1, -1, -1}
    };
    std::array<int, 4> m_iPiZero = {
        {-1, -1, -1, -1}
    };
    std::array<int, 4> m_iKPlus = {
        {-1, -1, -1, -1}
    };
    std::array<int, 4> m_iKMinus = {
        {-1, -1, -1, -1}
    };
};

#endif
