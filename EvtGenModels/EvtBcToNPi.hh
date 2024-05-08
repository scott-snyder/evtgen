
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

#ifndef EvtBcToNPi_HH
#define EvtBcToNPi_HH

#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtDecayAmp.hh"
#include "EvtGenBase/EvtDecayBase.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtVector4R.hh"

#include <string>

// Description: General decay model for Bc -> V + npi and Bc -> P + npi

class EvtBcToNPi : public EvtDecayAmp {
  public:
    EvtBcToNPi( bool printAuthorInfo = false );

    std::string getName() override;

    EvtDecayBase* clone() override;

    void initProbMax() override;

    void init() override;

    void decay( EvtParticle* p ) override;

  protected:
    int m_nCall;
    double m_maxAmp2;

    // Bc form factors
    double m_maxProb;
    double m_FA0_N, m_FA0_c1, m_FA0_c2;
    double m_FAm_N, m_FAm_c1, m_FAm_c2;
    double m_FAp_N, m_FAp_c1, m_FAp_c2;
    double m_FV_N, m_FV_c1, m_FV_c2;

    double m_Fp_N, m_Fp_c1, m_Fp_c2;
    double m_Fm_N, m_Fm_c1, m_Fm_c2;

    // W -> pi... form factors
    double m_beta;
    double m_mRho;
    double m_gammaRho;
    double m_mRhopr;
    double m_gammaRhopr;
    double m_mA1;
    double m_gammaA1;

    double energy1( double M, double m1, double m2 );
    double mom1( double M, double m1, double m2 );
    EvtComplex Fpi( EvtVector4R q1, EvtVector4R q2 );
    double pi3G( double m2, int dupD );

  private:
    void printAuthorInfo();
};

#endif
