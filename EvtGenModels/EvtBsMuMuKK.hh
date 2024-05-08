
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

#ifndef EvtBsMuMuKK_HH
#define EvtBsMuMuKK_HH

#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtDecayAmp.hh"

#include <string>

class EvtParticle;

// Description: Routine to implement Bs -> J/psi KK

class EvtBsMuMuKK : public EvtDecayAmp {
  public:
    std::string getName() override;
    EvtDecayBase* clone() override;

    void init() override;
    void initProbMax() override;

    void decay( EvtParticle* p ) override;

  protected:
    EvtComplex Flatte( const double m0, const double m ) const;

    EvtComplex GetRho( const double m0, const double m ) const;

    EvtComplex Breit_Wigner( const double Gamma0, const double m0,
                             const double m, const int J, const double q0,
                             const double q ) const;

    double Integral( const double Gamma0, const double m0, const int JR,
                     const int JB, const double q0, const double M_KK_ll,
                     const double M_KK_ul, const int fcntype ) const;

    double X_J( const int J, const double q, const int isB ) const;

    double Wignerd( int J, int l, int alpha, double theta ) const;

    EvtComplex AngularDist( int J, int l, int alpha, double cK, double cL,
                            double chi ) const;

    EvtComplex AmpTime( const int q, const EvtComplex& gplus,
                        const EvtComplex& gminus, const double delta,
                        const double lambda_abs, const double Amp,
                        const double phis, const int eta ) const;

  private:
    double m_MBs, m_MJpsi, m_Mf0, m_Mphi, m_Mf2p, m_MKp, m_MKm, m_MK0, m_Mpip,
        m_Mpi0, m_Mmu;
    double m_Gamma0phi, m_Gamma0f2p;
    double m_kin_lower_limit, m_kin_upper_limit, m_kin_middle;
    double m_p30Kp_mid_CMS, m_p30Kp_ll_CMS, m_p30Kp_phi_CMS, m_p30Kp_f2p_CMS;
    double m_p30Jpsi_mid_CMS, m_p30Jpsi_ll_CMS, m_p30Jpsi_phi_CMS,
        m_p30Jpsi_f2p_CMS;
    double m_int_const_NR, m_int_Flatte_f0, m_int_BW_phi, m_int_BW_f2p;
    double m_f_S_NR, m_f_f0, m_f_phi, m_f_f2p, m_f_phi_0, m_f_phi_perp,
        m_f_f2p_0, m_f_f2p_perp;
    double m_A_S_NR, m_A_f0, m_A_phi_0, m_A_phi_perp, m_A_phi_par, m_A_f2p_0,
        m_A_f2p_perp;
    double m_A_f2p_par;
    double m_delta_S_NR, m_delta_f0, m_delta_phi_0, m_delta_phi_perp,
        m_delta_phi_par;
    double m_delta_f2p_0, m_delta_f2p_perp, m_delta_f2p_par;
    double m_phis_S_NR, m_phis_f0, m_phis_phi_0, m_phis_phi_perp, m_phis_phi_par;
    double m_phis_f2p_0, m_phis_f2p_perp, m_phis_f2p_par;
    double m_lambda_S_NR_abs, m_lambda_f0_abs, m_lambda_phi_0_abs,
        m_lambda_phi_perp_abs;
    double m_lambda_phi_par_abs, m_lambda_f2p_0_abs, m_lambda_f2p_perp_abs;
    double m_lambda_f2p_par_abs;
    double m_Gamma, m_deltaGamma, m_ctau, m_deltaMs;
};

#endif
