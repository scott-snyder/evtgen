
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

#ifndef EVTBLLNUL_AMP_HH
#define EVTBLLNUL_AMP_HH

#include "EvtGenBase/EvtAmp.hh"
#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtTensor4C.hh"
#include "EvtGenBase/EvtVector4R.hh"

#include <vector>

class EvtParticle;

// Description: Header file for the amplitude calculation for the "BLLNUL"
//              model which generates rare four-leptonic B-decays
//              B^-(p) -> ell^+(k_1) ell^-(k_2) neu (k_3) ell^-(k_4)

class EvtBLLNuLAmp {
  public:
    EvtBLLNuLAmp( double Vub = 4.09e-3 );
    EvtBLLNuLAmp( double qSqMin, double kSqMin, bool symmetry,
                  double Vub = 4.09e-3 );

    void CalcAmp( EvtParticle* parent, EvtAmp& amp ) const;
    void setParameters( double qSqMin, double kSqMin, bool symmetry );

    // Resonance poles
    class ResPole final {
      public:
        ResPole( double mass, double width, double coupling );

        EvtComplex propagator( double qSq, int numForm = 0 ) const;

        double getMass() const { return m_m0; }
        double getMassSq() const { return m_m0Sq; }
        double getWidth() const { return m_w0; }
        double getCoupling() const { return m_c; }

      private:
        double m_m0;    // pole mass
        double m_m0Sq;
        double m_w0;    // width
        double m_c;     // coupling constant
        EvtComplex m_I;
        EvtComplex m_Imw;
    };

  protected:
    EvtTensor4C getHadronTensor( const EvtVector4R& q, const EvtVector4R& k,
                                 const double qSq, const double kSq,
                                 const double MB, const int sign ) const;

    std::vector<EvtComplex> getVMDTerms( double qSq, double kSq, double MB ) const;

    EvtComplex getBStarTerm( double qSq, double kSq, double MB ) const;

    double FF_B2Bstar( double qSq ) const;

    double FF_V( double kSq ) const;

    double FF_A1( double kSq ) const;

    double FF_A2( double kSq ) const;

  private:
    // Kinematic cut-offs
    double m_qSqMin;
    double m_kSqMin;

    // If we have identical charged lepton flavours
    bool m_symmetry;

    // B+, B- Ids
    EvtId m_BpId, m_BnId;

    // Form factor constants
    double m_coupling, m_sqrt2;
    double m_fBu;

    // Resonance poles
    EvtBLLNuLAmp::ResPole m_Bstar, m_Upsilon;

    std::vector<EvtBLLNuLAmp::ResPole> m_resPoles;
    int m_nPoles;

    // Complex number constants
    EvtComplex m_zero, m_unitI;
};

inline void EvtBLLNuLAmp::setParameters( double qSqMin, double kSqMin,
                                         bool symmetry )
{
    m_qSqMin = qSqMin;
    m_kSqMin = kSqMin;
    m_symmetry = symmetry;
}

#endif
