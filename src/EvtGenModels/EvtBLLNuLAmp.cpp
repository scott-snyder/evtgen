
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

#include "EvtGenModels/EvtBLLNuLAmp.hh"

#include "EvtGenBase/EvtConst.hh"
#include "EvtGenBase/EvtDiracSpinor.hh"
#include "EvtGenBase/EvtIdSet.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtVector4C.hh"

#include <cmath>

EvtBLLNuLAmp::EvtBLLNuLAmp( double Vub ) :
    m_qSqMin( 0.0 ),
    m_kSqMin( 0.0 ),
    m_symmetry( false ),
    m_BpId( EvtPDL::getId( "B+" ) ),
    m_BnId( EvtPDL::getId( "B-" ) ),
    m_coupling( 0.0 ),
    m_sqrt2( sqrt( 2.0 ) ),
    m_fBu( 0.191 ),    // leptonic constant (GeV)
    m_Bstar( EvtBLLNuLAmp::ResPole( 5.32, 0.00658, 0.183 / 3.0 ) ),
    m_Upsilon( EvtBLLNuLAmp::ResPole( 9.64, 0.0, 0.0 ) ),
    m_resPoles(),
    m_nPoles( 0 ),
    m_zero( EvtComplex( 0.0, 0.0 ) ),
    m_unitI( EvtComplex( 0.0, 1.0 ) )
{
    double GF = 1.166371e-5;    // GeV^{-2}
    double alphaEM = 1.0 / 137.0;

    // Normalisation constant, multiplied by 1e4 to increase probability scale
    m_coupling = 400.0 * GF * EvtConst::pi * alphaEM * Vub * 1e4 / m_sqrt2;

    // Define VMD resonance poles using PDG 2016 values with constants from
    // D.Melikhov, N.Nikitin and K.Toms, Phys. Atom. Nucl. 68, 1842 (2005)

    // Rho and omega resonances
    EvtBLLNuLAmp::ResPole rho = EvtBLLNuLAmp::ResPole( 0.77526, 0.1491,
                                                       1.0 / 5.04 );
    m_resPoles.push_back( rho );

    EvtBLLNuLAmp::ResPole omega = EvtBLLNuLAmp::ResPole( 0.78265, 0.00849,
                                                         1.0 / 17.1 );
    m_resPoles.push_back( omega );

    m_nPoles = m_resPoles.size();
}

EvtBLLNuLAmp::EvtBLLNuLAmp( double qSqMin, double kSqMin, bool symmetry,
                            double Vub ) :
    m_qSqMin( qSqMin ),
    m_kSqMin( kSqMin ),
    m_symmetry( symmetry ),
    m_BpId( EvtPDL::getId( "B+" ) ),
    m_BnId( EvtPDL::getId( "B-" ) ),
    m_coupling( 0.0 ),
    m_sqrt2( sqrt( 2.0 ) ),
    m_fBu( 0.191 ),    // leptonic constant (GeV)
    m_Bstar( EvtBLLNuLAmp::ResPole( 5.32, 0.00658, 0.183 / 3.0 ) ),
    m_Upsilon( EvtBLLNuLAmp::ResPole( 9.64, 0.0, 0.0 ) ),
    m_resPoles(),
    m_nPoles( 0 ),
    m_zero( EvtComplex( 0.0, 0.0 ) ),
    m_unitI( EvtComplex( 0.0, 1.0 ) )
{
    double GF = 1.166371e-5;    // GeV^{-2}
    double alphaEM = 1.0 / 137.0;

    // Normalisation constant, multiplied by 1e4 to increase probability scale
    m_coupling = 400.0 * GF * EvtConst::pi * alphaEM * Vub * 1e4 / m_sqrt2;

    // Define VMD resonance poles using PDG 2016 values with constants from
    // D.Melikhov, N.Nikitin and K.Toms, Phys. Atom. Nucl. 68, 1842 (2005)

    // Rho and omega resonances
    EvtBLLNuLAmp::ResPole rho = EvtBLLNuLAmp::ResPole( 0.77526, 0.1491,
                                                       1.0 / 5.04 );
    m_resPoles.push_back( rho );

    EvtBLLNuLAmp::ResPole omega = EvtBLLNuLAmp::ResPole( 0.78265, 0.00849,
                                                         1.0 / 17.1 );
    m_resPoles.push_back( omega );

    m_nPoles = m_resPoles.size();
}

// Storing resonance pole information
EvtBLLNuLAmp::ResPole::ResPole( double mass, double width, double coupling ) :
    m_m0( mass ),
    m_m0Sq( mass * mass ),
    m_w0( width ),
    m_c( coupling ),
    m_I( EvtComplex( 0.0, 1.0 ) ),
    m_Imw( m_I * mass * width )
{
}

EvtComplex EvtBLLNuLAmp::ResPole::propagator( double qSq, int numForm ) const
{
    // Numerator term: mass-squared (default) or mass
    double num( m_m0Sq );
    if ( numForm == 1 ) {
        num = m_m0;
    }

    EvtComplex result = num * m_c / ( ( qSq - m_m0Sq ) + m_Imw );
    return result;
}

// Amplitude calculation
void EvtBLLNuLAmp::CalcAmp( EvtParticle* parent, EvtAmp& amp ) const
{
    // Check for 4 daughters and an existing parent
    if ( !parent || parent->getNDaug() != 4 ) {
        return;
    }

    // The first two charged leptons. The 2nd one will have
    // the same charge as the 3rd charged lepton
    EvtParticle* lepA = parent->getDaug( 0 );
    EvtParticle* lepB = parent->getDaug( 1 );
    // The neutrino
    EvtParticle* neu = parent->getDaug( 2 );
    // The third charged lepton
    EvtParticle* lepC = parent->getDaug( 3 );

    // Kinematics
    double MB = parent->mass();    // B-meson mass, GeV

    // 4-momenta of the leptons in the B rest frame. The daughters will already
    // be in the correct order since this check is done in EvtBLLNuL::init()
    // when initialising the model using the decay file
    EvtVector4R p1 = lepA->getP4();
    EvtVector4R p2 = lepB->getP4();
    EvtVector4R p3 = neu->getP4();
    EvtVector4R p4 = lepC->getP4();

    // 4-momenta sums
    EvtVector4R q12 = p1 + p2;
    EvtVector4R k34 = p3 + p4;

    // Mandelstam variables: q^2 and k^2
    double q12Sq = q12.mass2();
    double k34Sq = k34.mass2();

    // Check if we are above mass thresholds
    bool threshold( true );
    if ( q12Sq < m_qSqMin || k34Sq < m_kSqMin ) {
        threshold = false;
    }

    // For the symmetric terms when we exchange the
    // 2nd and 3rd charged leptons: p2 <-> p4
    EvtVector4R q14, k23;
    double q14Sq( 0.0 ), k23Sq( 0.0 );
    if ( m_symmetry ) {
        q14 = p1 + p4;
        k23 = p2 + p3;
        q14Sq = q14.mass2();
        k23Sq = k23.mass2();

        if ( q14Sq < m_qSqMin || k23Sq < m_kSqMin ) {
            threshold = false;
        }
    }

    // B meson id
    EvtId parId = parent->getId();
    // B+ or B- decays
    int sign( 1 );
    if ( parId == m_BnId ) {
        sign = -1;
    }

    // Hadronic tensors
    EvtTensor4C THadronA = getHadronTensor( q12, k34, q12Sq, k34Sq, MB, sign );

    // When we need to include the symmetric terms
    EvtTensor4C THadronB;
    if ( m_symmetry ) {
        THadronB = getHadronTensor( q14, k23, q14Sq, k23Sq, MB, sign );
    }

    // Leptonic currents: A for normal terms, B for symmetric terms
    EvtVector4C L1A, L2A, L1B, L2B;

    int leptonSpins[4];    // array for saving the leptonic spin configuration

    // Loop over lepton spin states
    for ( int i2 = 0; i2 < 2; i2++ ) {
        leptonSpins[0] = i2;

        for ( int i1 = 0; i1 < 2; i1++ ) {
            leptonSpins[1] = i1;

            if ( sign == -1 ) {
                // B- currents
                // L2^{\nu} = \bar mu(k_2) \gamma^{\nu} mu(- k_1)
                L2A = EvtLeptonVCurrent( lepB->spParent( i2 ),
                                         lepA->spParent( i1 ) );

                if ( m_symmetry ) {
                    // Swapping the 2nd and 3rd charged leptons
                    L1B = EvtLeptonVACurrent( lepB->spParent( i2 ),
                                              neu->spParentNeutrino() );
                }

            } else {
                // B+ currents
                // L2^{\nu} = \bar mu(k_1) \gamma^{\nu} mu(- k_2)
                L2A = EvtLeptonVCurrent( lepA->spParent( i1 ),
                                         lepB->spParent( i2 ) );

                if ( m_symmetry ) {
                    // Swapping the 2nd and 3rd charged leptons
                    L1B = EvtLeptonVACurrent( neu->spParentNeutrino(),
                                              lepB->spParent( i2 ) );
                }
            }

            // Production: Tfi^{\mu} = THadron^{\mu \nu} L_{2 \nu}
            EvtVector4C THL2A = THadronA.cont2( L2A );

            for ( int i4 = 0; i4 < 2; i4++ ) {
                leptonSpins[2] = i4;
                leptonSpins[3] = 0;    // neutrino handedness

                if ( sign == -1 ) {
                    // B- currents
                    // L1^{\mu} = \bar e(k_4) \gamma^{\mu} (1 - \gamma^5) nu_e(- k_3)
                    L1A = EvtLeptonVACurrent( lepC->spParent( i4 ),
                                              neu->spParentNeutrino() );

                    if ( m_symmetry ) {
                        // Swapping the 2nd and 3rd charged leptons
                        L2B = EvtLeptonVCurrent( lepC->spParent( i4 ),
                                                 lepA->spParent( i1 ) );
                    }

                } else {
                    // B+ currents
                    // L1^{\mu} = \bar nu_e(k_3) \gamma^{\mu} (1 - \gamma^5) e(- k_4)
                    L1A = EvtLeptonVACurrent( neu->spParentNeutrino(),
                                              lepC->spParent( i4 ) );

                    if ( m_symmetry ) {
                        // Swapping the 2nd and 3rd charged leptons
                        L2B = EvtLeptonVCurrent( lepA->spParent( i1 ),
                                                 lepC->spParent( i4 ) );
                    }
                }

                if ( threshold == false ) {
                    // Below kinematic thresholds
                    amp.vertex( leptonSpins, m_zero );

                } else {
                    // Decay amplitude calculation: L_1^{\mu} Tfi_{\mu}
                    EvtComplex decAmp = L1A * THL2A;

                    // If we also need to swap the 2nd and 3rd charged leptons
                    if ( m_symmetry ) {
                        // Hadronic current production term. L2B depends on i4 so we need
                        // it here instead of inside the i2 loop as was the case for THL2A
                        EvtVector4C THL2B = THadronB.cont2( L2B );

                        // The symmetric amplitude
                        EvtComplex ampB = L1B * THL2B;

                        // Subtract this from the total amplitude
                        decAmp -= ampB;
                    }

                    amp.vertex( leptonSpins, decAmp );
                }

            }    // i4 loop

        }    // i1 loop

    }    // i2 loop
}

EvtTensor4C EvtBLLNuLAmp::getHadronTensor( const EvtVector4R& q,
                                           const EvtVector4R& k,
                                           const double qSq, const double kSq,
                                           const double MB, const int sign ) const
{
    // Hadronic tensor calculation

    EvtTensor4C epskq = dual( EvtGenFunctions::directProd( k, q ) );
    EvtTensor4C qk = EvtGenFunctions::directProd( q, k );

    EvtComplex BstarAmp = getBStarTerm( qSq, kSq, MB );
    std::vector<EvtComplex> VMDAmps = getVMDTerms( qSq, kSq, MB );

    EvtComplex FF_ekq = BstarAmp + VMDAmps[0];
    EvtComplex FF_g = VMDAmps[1] - m_fBu;
    EvtComplex FF_qk = VMDAmps[2];

    // Full hadronic tensor
    EvtTensor4C THadron = sign * 2.0 * FF_ekq * epskq +
                          m_unitI *
                              ( 2.0 * FF_qk * qk - FF_g * EvtTensor4C::g() );

    // Kinematic cuts
    double coeffcut( 0.0 );
    if ( qSq > m_qSqMin && kSq > m_kSqMin ) {
        coeffcut = 1.0 / qSq;
    }

    // Normalisation constant
    THadron *= coeffcut * m_coupling;

    return THadron;
}

std::vector<EvtComplex> EvtBLLNuLAmp::getVMDTerms( double qSq, double kSq,
                                                   double MB ) const
{
    // Find the 3 VMD form factors: epsilon*k*q, g(uv) and q*k terms
    EvtComplex VMD1( 0.0, 0.0 ), VMD2( 0.0, 0.0 ), VMD3( 0.0, 0.0 );

    // Loop over the VMD poles
    for ( int iPole = 0; iPole < m_nPoles; iPole++ ) {
        auto pole = m_resPoles[iPole];

        // Propagator term, common for all factors
        EvtComplex prop = pole.propagator( qSq );

        double mSum = MB + pole.getMass();

        VMD1 += prop / mSum;
        VMD2 += mSum * prop;
    }

    // Third pole summation term is the same as the first one
    VMD3 = VMD1;

    // Multiply by couplings for the given kSq
    VMD1 *= FF_V( kSq );
    VMD2 *= FF_A1( kSq );
    VMD3 *= FF_A2( kSq );

    // Return the factors as a vector
    std::vector<EvtComplex> factors;
    factors.push_back( VMD1 );
    factors.push_back( VMD2 );
    factors.push_back( VMD3 );

    return factors;
}

EvtComplex EvtBLLNuLAmp::getBStarTerm( double qSq, double kSq, double MB ) const
{
    EvtComplex amplitude = m_Bstar.propagator( kSq, 1 ) * FF_B2Bstar( qSq ) /
                           ( MB + m_Bstar.getMass() );
    return amplitude;
}

double EvtBLLNuLAmp::FF_B2Bstar( double qSq ) const
{
    // Electromagnetic FF for B -> B* transition, when gamma is emitted from the b quark
    // D.Melikhov, private communication
    double y = qSq / m_Upsilon.getMassSq();
    double denom = ( 1.0 - y ) * ( 1.0 - 0.81 * y );

    double V( 0.0 );
    if ( fabs( denom ) > 1e-10 ) {
        V = 1.044 / denom;
    }

    return V;
}

double EvtBLLNuLAmp::FF_V( double kSq ) const
{
    // D. Melikhov and B. Stech, PRD 62, 014006 (2000) Table XV
    double y = kSq / m_Bstar.getMassSq();
    double denom = m_sqrt2 * ( 1.0 - y ) * ( 1.0 - 0.59 * y );

    double V( 0.0 );
    if ( fabs( denom ) > 1e-10 ) {
        V = 0.31 / denom;
    }

    return V;
}

double EvtBLLNuLAmp::FF_A1( double kSq ) const
{
    // D. Melikhov and B. Stech, PRD 62, 014006 (2000) Table XV
    double y = kSq / m_Bstar.getMassSq();
    double denom = ( ( 0.1 * y - 0.73 ) * y + 1.0 ) * m_sqrt2;

    double A1( 0.0 );
    if ( fabs( denom ) > 1e-10 ) {
        A1 = 0.26 / denom;
    }

    return A1;
}

double EvtBLLNuLAmp::FF_A2( double kSq ) const
{
    // D. Melikhov and B. Stech, PRD 62, 014006 (2000) Table XV
    double y = kSq / m_Bstar.getMassSq();
    double denom = ( ( 0.5 * y - 1.4 ) * y + 1.0 ) * m_sqrt2;

    double A2( 0.0 );
    if ( fabs( denom ) > 1e-10 ) {
        A2 = 0.24 / denom;
    }

    return A2;
}
