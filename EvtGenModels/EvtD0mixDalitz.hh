
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

#ifndef __EVTD0MIXDALITZ_HH__
#define __EVTD0MIXDALITZ_HH__

#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtCyclic3.hh"
#include "EvtGenBase/EvtDalitzPoint.hh"
#include "EvtGenBase/EvtDalitzReso.hh"
#include "EvtGenBase/EvtDecayAmp.hh"
#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtSpinType.hh"

// Description:
//   The D0mixDalitz model, with many resonances and mixing implemented.

class EvtD0mixDalitz : public EvtDecayAmp {
  private:
    int m_d1;
    int m_d2;
    int m_d3;

    // Mixing parameters.
    double m_x;
    double m_y;

    // q/p CP violation in the mixing.
    EvtComplex m_qp;

    // Checker of the decay mode.
    bool m_isKsPiPi;
    bool m_isRBWmodel;

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

    // Masses of the relevant particles.
    double m_mD0;
    double m_mKs;
    double m_mPi;
    double m_mK;

    // Life time and decay rate.
    double m_ctau;
    double m_gamma;

    // Some useful integrals over the Dalitz plot.
    EvtComplex m_iChi;
    EvtComplex m_iChi2;

    void readPDGValues();
    EvtComplex dalitzKsPiPi( const EvtDalitzPoint& point );
    EvtComplex dalitzKsKK( const EvtDalitzPoint& point );

    // Time evolution functions for hamiltonian eigenstates.
    //    Negative exponential part removed.
    EvtComplex h1( const double& ct ) const;
    EvtComplex h2( const double& ct ) const;

    void reportInvalidAndExit() const
    {
        EvtGenReport( EVTGEN_ERROR, "EvtD0mixDalitz" )
            << "EvtD0mixDalitz: Invalid mode." << std::endl;
        exit( 1 );
    }

  public:
    EvtD0mixDalitz() :
        m_d1( 0 ),
        m_d2( 0 ),
        m_d3( 0 ),
        m_x( 0. ),
        m_y( 0. ),
        m_qp( 1. ),
        m_isKsPiPi( false ),
        m_isRBWmodel( true )
    {
    }

    // One-line inline functions.
    std::string getName() override { return "D0MIXDALITZ"; }
    EvtDecayBase* clone() override { return new EvtD0mixDalitz; }
    void initProbMax() override { setProbMax( 5200. ); }

    void init() override;
    void decay( EvtParticle* p ) override;
};

#endif
