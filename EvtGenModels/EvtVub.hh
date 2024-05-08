
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

#ifndef EVTVUB_HH
#define EVTVUB_HH

#include "EvtGenBase/EvtDecayIncoherent.hh"

#include "EvtGenModels/EvtVubdGamma.hh"

#include <memory>
#include <vector>

class EvtParticle;

// Description:
// Class to generate inclusive B to X_u l nu decays according to various
// decay models. Implemtented are ACCM, parton-model and a QCD model.

class EvtVub : public EvtDecayIncoherent {
  public:
    std::string getName() override;

    EvtDecayBase* clone() override;

    void initProbMax() override;

    void init() override;

    void decay( EvtParticle* p ) override;

  private:
    double m_mb;        // the b-quark pole mass in GeV (try 4.65 to 4.9)
    double m_a;         // Parameter for the Fermi Motion (1.29 is good)
    double m_alphas;    // Strong Coupling at m_b (around 0.24)
    double m_dGMax;     // max dGamma*p2 value;
    int m_nbins;
    int m_storeQplus;
    std::vector<double> m_masses;
    std::vector<double> m_weights;

    std::unique_ptr<EvtVubdGamma> m_dGamma;    // calculates the decay rate
    double findPFermi();
    std::vector<double> m_pf;
};

#endif
