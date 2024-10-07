
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

#ifndef EVTSVSCPISO_HH
#define EVTSVSCPISO_HH

#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtDecayAmp.hh"

class EvtParticle;

// Description: Routine to decay scalar -> vectors scalar
//              with CP violation and isospin amplitudes.
//              More specifically, it is indended to handle
//              decays like B->rho pi and B->a1 pi.

class EvtSVSCPiso : public EvtDecayAmp {
  public:
    std::string getName() const override;
    EvtDecayBase* clone() const override;

    void init() override;
    void initProbMax() override;

    void decay( EvtParticle* p ) override;

  private:
    // Amplitude coefficients
    EvtComplex m_Tp0, m_Tp0_bar, m_T0p, m_T0p_bar;
    EvtComplex m_Tpm, m_Tpm_bar, m_Tmp, m_Tmp_bar;
    EvtComplex m_P1, m_P1_bar, m_P0, m_P0_bar;

    // Amplitudes
    EvtComplex m_A_f, m_Abar_f;
    EvtComplex m_A_fbar, m_Abar_fbar;
    EvtComplex m_Apm, m_Apm_bar, m_Amp, m_Amp_bar;
    // Charged mode flag
    int m_charged{ 0 };

    // Set amplitude coeffs from decay model pars
    void setAmpCoeffs();
    // Calculate amplitude terms
    void calcAmpTerms();
};

#endif
