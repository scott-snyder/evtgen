
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

#ifndef EVTLAMBDACPHH_HH
#define EVTLAMBDACPHH_HH

#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtDecayAmp.hh"
#include "EvtGenBase/EvtResonance2.hh"
#include "EvtGenBase/EvtVector4R.hh"

#include <string>
#include <vector>

class EvtParticle;

// Description: Decay model for Lambda_c -> K- pi+ p using amplitudes
//              from the Fermilab E791 analysis: arXiv:hep-ex/9912003v1

class EvtLambdacPHH : public EvtDecayAmp {
  public:
    EvtLambdacPHH();

    std::string getName() const override;
    EvtDecayBase* clone() const override;

    void init() override;
    void initProbMax() override;
    void decay( EvtParticle* p ) override;

  protected:
    // Resonance enumeration
    enum class LcResLabel
    {
        NonReson = 0,
        Kstar,
        Delta,
        Lambda
    };

    // Amplitude functions
    std::vector<EvtComplex> calcResAmpTerms( EvtLambdacPHH::LcResLabel resIndex,
                                             const EvtResonance2& res,
                                             double norm ) const;

    EvtComplex DecayAmp3( EvtLambdacPHH::LcResLabel resonance, int m,
                          int mprime, double theta_res, double phi_res,
                          double theta_prime_daughter_res,
                          double phi_prime_daughter_res ) const;

    EvtComplex fampl3( double amplitude_res, double phi_res, int spinMother,
                       int m_spinMother, int m_prime_spinMother,
                       double theta_res, float spin_res, float m_spin_res,
                       float m_prime_spin_res, double theta_daughter_res,
                       double phi_prime_daughter_res ) const;

    // Find resonance normalisation factors
    void calcNormalisations();

    void getFitFractions();

    // Inverse cos/sin functions that checks for valid arguments
    double getACos( double num, double denom ) const;
    double getASin( double num, double denom ) const;

  private:
    // Daughter ordering for K-, pi+, p
    int m_d1, m_d2, m_d3;

    // Resonance parameters
    double m_Nplusplus, m_Nplusminus, m_Nminusplus, m_Nminusminus;
    double m_phiNplusplus, m_phiNplusminus, m_phiNminusplus, m_phiNminusminus;
    double m_E1, m_phiE1, m_E2, m_phiE2, m_E3, m_phiE3, m_E4, m_phiE4;
    double m_F1, m_phiF1, m_F2, m_phiF2, m_H1, m_phiH1, m_H2, m_phiH2;

    double m_NRNorm, m_KstarNorm, m_DeltaNorm, m_LambdaNorm;
    double m_KstarM, m_KstarW, m_KstarR;
    double m_DeltaM, m_DeltaW, m_DeltaR;
    double m_LambdaM, m_LambdaW, m_LambdaR;
    double m_Lambda_cR;

    EvtVector4R m_zprime, m_p4_Lambda_c;
    double m_zpMag, m_p4_Lambdac_Mag;
};

#endif
