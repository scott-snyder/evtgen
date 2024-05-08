
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

#ifndef __EVTD0GAMMADALITZ_HH__
#define __EVTD0GAMMADALITZ_HH__

#include "EvtGenBase/EvtCyclic3.hh"
#include "EvtGenBase/EvtDalitzReso.hh"
#include "EvtGenBase/EvtDecayAmp.hh"
#include "EvtGenBase/EvtFlatte.hh"
#include "EvtGenBase/EvtSpinType.hh"

#include <vector>

class EvtParticle;

class EvtD0gammaDalitz : public EvtDecayAmp {
  private:
    int m_d1;
    int m_d2;
    int m_d3;

    bool m_isKsPiPi;

    // Useful constants.
    static const EvtSpinType::spintype& m_SCALAR;
    static const EvtSpinType::spintype& m_VECTOR;
    static const EvtSpinType::spintype& m_TENSOR;

    static const EvtDalitzReso::CouplingType& m_EtaPic;
    static const EvtDalitzReso::CouplingType& m_PicPicKK;

    static const EvtDalitzReso::NumType& m_RBW;
    static const EvtDalitzReso::NumType& m_GS;
    static const EvtDalitzReso::NumType& m_KMAT;

    static const EvtCyclic3::Pair& m_AB;
    static const EvtCyclic3::Pair& m_AC;
    static const EvtCyclic3::Pair& m_BC;

    // Values to be read or computed based on values in the evt.pdl file.
    // IDs of the relevant particles.
    EvtId m_BP;
    EvtId m_BM;
    EvtId m_B0;
    EvtId m_B0B;
    EvtId m_D0;
    EvtId m_D0B;
    EvtId m_KM;
    EvtId m_KP;
    EvtId m_K0;
    EvtId m_K0B;
    EvtId m_KL;
    EvtId m_KS;
    EvtId m_PIM;
    EvtId m_PIP;

    // Flavor of the B mother.
    EvtId m_bFlavor;

    // Masses of the relevant particles.
    double m_mD0;
    double m_mKs;
    double m_mPi;
    double m_mK;

    void readPDGValues();
    void reportInvalidAndExit() const;

    EvtComplex dalitzKsPiPi( const EvtDalitzPoint& point ) const;
    EvtComplex dalitzKsKK( const EvtDalitzPoint& point ) const;

  public:
    std::string getName() override;
    EvtDecayBase* clone() override;

    void init() override;
    void initProbMax() override;

    void decay( EvtParticle* p ) override;
};

#endif
