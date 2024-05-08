
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

#ifndef EVTBBSCALAR_HH
#define EVTBBSCALAR_HH

#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtConst.hh"
#include "EvtGenBase/EvtDecayAmp.hh"
#include "EvtGenBase/EvtDiracParticle.hh"
#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtScalarParticle.hh"
#include "EvtGenBase/EvtVector4C.hh"

#include <bitset>
#include <map>
#include <string>
#include <vector>

// Implementation of the decay B- -> lambda p_bar pi according to
// hep-ph/0204185, hep-ph/0211240
// This model is intended to be applicable to all decays of the type B-> baryon baryon scalar

class EvtBBScalar : public EvtDecayAmp {
  public:
    EvtBBScalar();
    std::string getName() override;
    EvtBBScalar* clone() override;
    void decay( EvtParticle* p ) override;
    void init() override;
    void initProbMax() override;

  private:
    struct FormFactor {
        double m_value;
        double m_sigma1;
        double m_sigma2;
        double m_mV;
    };

    enum Baryons
    {
        Lambda,
        Proton,
        Neutron,
        Sigma0,
        Sigma_minus,
        Xi0,
        Xi_minus,
        nBaryons
    };

    // used values of constants
    static const EvtComplex m_I;
    static const EvtComplex m_V_ub;
    static const EvtComplex m_V_us_star;
    static const EvtComplex m_a1;
    static const EvtComplex m_V_tb;
    static const EvtComplex m_V_ts_star;
    static const EvtComplex m_a4;
    static const EvtComplex m_a6;

    // used parameters in the calculation of the magnetic form factors
    static const double m_x[];
    static const double m_y[];
    // quark masses as used in the model
    static const double m_ms;
    static const double m_mu;
    static const double m_mb;

    // used to choose the right m_value for the form factor depending on the type of scalar
    std::string m_scalarType;
    mutable std::map<std::string, FormFactor> m_f0Map;
    mutable std::map<std::string, FormFactor> m_f1Map;

    // only consider F1+F2 here
    std::bitset<nBaryons> m_baryonCombination;
    void setKnownBaryonTypes( const EvtId& baryon );

    double B_pi_f1( double t ) const;
    double B_pi_f0( double t ) const;
    double baryonF1F2( double t ) const;
    double G_p( double t ) const;
    double G_n( double t ) const;

    double baryon_gA( double t ) const;
    double baryon_hA( double t ) const;
    double baryon_gP( double t ) const;
    double baryon_fS( double t ) const;

    double D_A( double t ) const;
    double F_A( double t ) const;
    double D_P( double t ) const;
    double F_P( double t ) const;
    double D_S( double t ) const;
    double F_S( double t ) const;

    // (mB1 - mB2)/(mq1 - mq1)
    double m_massRatio;
    double m_baryonMassSum;
    double formFactorFit( double t, const std::vector<double>& params ) const;

    static const EvtComplex m_const_B;
    static const EvtComplex m_const_C;
    const EvtVector4C amp_A( const EvtVector4R& p4B, const EvtVector4R& p4Scalar );
    const EvtComplex amp_B( const EvtDiracParticle* baryon1,
                            const EvtDiracSpinor& b1Pol,
                            const EvtDiracParticle* baryon2,
                            const EvtDiracSpinor& b2Pol, int index );
    const EvtComplex amp_B_vectorPart( const EvtDiracParticle* baryon1,
                                       const EvtDiracSpinor& b1Pol,
                                       const EvtDiracParticle* baryon2,
                                       const EvtDiracSpinor& b2Pol, int index );
    const EvtComplex amp_B_axialPart( const EvtDiracParticle* baryon1,
                                      const EvtDiracSpinor& b1Pol,
                                      const EvtDiracParticle* baryon2,
                                      const EvtDiracSpinor& b2Pol, int index );
    const EvtComplex amp_C( const EvtDiracParticle* baryon1,
                            const EvtDiracSpinor& b1Pol,
                            const EvtDiracParticle* baryon2,
                            const EvtDiracSpinor& b2Pol, int index );
    const EvtComplex amp_C_scalarPart( const EvtDiracSpinor& b1Pol,
                                       const EvtDiracSpinor& b2Pol, double t );
    const EvtComplex amp_C_pseudoscalarPart( const EvtDiracSpinor& b1Pol,
                                             const EvtDiracSpinor& b2Pol,
                                             double t );

    // initialize phasespace and calculate the amplitude for one (i=0,1) state of the photon
    EvtComplex calcAmpliude( const EvtParticle* p, const unsigned int polState );
};

#endif
