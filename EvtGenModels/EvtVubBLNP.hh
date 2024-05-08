
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

#ifndef EVTVUBBLNP_HH
#define EVTVUBBLNP_HH

#include "EvtGenBase/EvtDecayIncoherent.hh"

#include <vector>

class EvtParticle;

// Description: Modeled on Riccardo Faccini's EvtVubNLO module
// tripleDiff from BLNP's notebook (based on BLNP4, hep-ph/0504071)

class EvtVubBLNP : public EvtDecayIncoherent {
  public:
    std::string getName() override;

    EvtDecayBase* clone() override;

    void initProbMax() override;

    void init() override;

    void decay( EvtParticle* Bmeson ) override;

  private:
    // Input parameters
    double m_mBB;
    double m_lambda2;

    // Shape function parameters
    double m_b;
    double m_Lambda;
    double m_Ecut;
    double m_wzero;

    // SF and SSF modes
    int m_itype;
    double m_dtype;
    int m_isubl;

    // flags
    int m_flag1;
    int m_flag2;
    int m_flag3;

    // Quark mass
    double m_mb;

    // Matching scales
    double m_muh;
    double m_mui;
    double m_mubar;

    // Perturbative quantities
    double m_CF;
    double m_CA;

    double m_beta0;
    double m_beta1;
    double m_beta2;

    double m_zeta3;

    double m_Gamma0;
    double m_Gamma1;
    double m_Gamma2;

    double m_gp0;
    double m_gp1;

    double m_Lbar;
    double m_mupisq;
    double m_moment2;

    int m_flagpower;
    int m_flag2loop;

    int m_maxLoop;
    double m_precision;

    std::vector<double> m_gvars;

    double rate3( double Pp, double Pl, double Pm );
    double F1( double Pp, double Pm, double muh, double mui, double mubar,
               double doneJS, double done1 );
    double F2( double Pp, double Pm, double muh, double mui, double mubar,
               double done3 );
    double F3( double Pp, double Pm, double muh, double mui, double mubar,
               double done2 );
    double DoneJS( double Pp, double Pm, double mui );
    double Done1( double Pp, double Pm, double mui );
    double Done2( double Pp, double Pm, double mui );
    double Done3( double Pp, double Pm, double mui );
    static double IntJS( double what, const std::vector<double>& vars );
    static double Int1( double what, const std::vector<double>& vars );
    static double Int2( double what, const std::vector<double>& vars );
    static double Int3( double what, const std::vector<double>& vars );
    static double g1( double w, const std::vector<double>& vars );
    static double g2( double w, const std::vector<double>& vars );
    static double g3( double w, const std::vector<double>& vars );
    static double Shat( double w, const std::vector<double>& vars );
    static double Mzero( double muf, double mu, double mupisq,
                         const std::vector<double>& vars );
    double wS( double w );
    double t( double w );
    double u( double w );
    double v( double w );
    double myfunction( double w, double Lbar, double mom2 );
    double myfunctionBIK( double w, double Lbar, double mom2 );
    double dU1nlo( double muh, double mui );
    double U1lo( double muh, double mui );
    double Sfun( double mu1, double mu2, double epsilon );
    double S0( double a1, double r );
    double S1( double a1, double r );
    double S2( double a1, double r );
    double aGamma( double mu1, double mu2, double epsilon );
    double agp( double mu1, double mu2, double epsilon );
    double alo( double muh, double mui );
    double anlo( double muh, double mui );    // d/depsilon of aGamma
    static double alphas( double mu, const std::vector<double>& vars );
    double PolyLog( double v, double z );
    static double Gamma( double z );
    static double Gamma( double a, double x );
    static double gamser( double a, double x, double LogGamma );
    static double gammcf( double a, double x, double LogGamma );
    double findBLNPWhat();
    std::vector<double> m_pf;
};

#endif
