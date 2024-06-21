
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

#include "EvtGenModels/EvtPartWave.hh"

#include "EvtGenBase/EvtCGCoefSingle.hh"
#include "EvtGenBase/EvtConst.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtTensor4C.hh"
#include "EvtGenBase/EvtVector4C.hh"

#include <algorithm>
#include <stdlib.h>
#include <string>
using std::endl;

std::string EvtPartWave::getName()
{
    return "PARTWAVE";
}

EvtDecayBase* EvtPartWave::clone()
{
    return new EvtPartWave;
}

void EvtPartWave::init()
{
    checkNDaug( 2 );

    //find out how many states each particle have
    int nA = EvtSpinType::getSpinStates( EvtPDL::getSpinType( getParentId() ) );
    int nB = EvtSpinType::getSpinStates( EvtPDL::getSpinType( getDaug( 0 ) ) );
    int nC = EvtSpinType::getSpinStates( EvtPDL::getSpinType( getDaug( 1 ) ) );

    if ( verbose() ) {
        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << "nA,nB,nC:" << nA << "," << nB << "," << nC << endl;
    }

    //find out what 2 times the spin is
    int JA2 = EvtSpinType::getSpin2( EvtPDL::getSpinType( getParentId() ) );
    int JB2 = EvtSpinType::getSpin2( EvtPDL::getSpinType( getDaug( 0 ) ) );
    int JC2 = EvtSpinType::getSpin2( EvtPDL::getSpinType( getDaug( 1 ) ) );

    if ( verbose() ) {
        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << "JA2,JB2,JC2:" << JA2 << "," << JB2 << "," << JC2 << endl;
    }

    //allocate memory
    int* lambdaA2 = new int[nA];
    int* lambdaB2 = new int[nB];
    int* lambdaC2 = new int[nC];

    EvtComplexPtr* HBC = new EvtComplexPtr[nB];
    for ( int ib = 0; ib < nB; ib++ ) {
        HBC[ib] = new EvtComplex[nC];
    }

    //find the allowed helicities (actually 2*times the helicity!)

    fillHelicity( lambdaA2, nA, JA2 );
    fillHelicity( lambdaB2, nB, JB2 );
    fillHelicity( lambdaC2, nC, JC2 );

    if ( verbose() ) {
        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << "Helicity states of particle A:" << endl;
        for ( int i = 0; i < nA; i++ ) {
            EvtGenReport( EVTGEN_INFO, "EvtGen" ) << lambdaA2[i] << endl;
        }

        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << "Helicity states of particle B:" << endl;
        for ( int i = 0; i < nB; i++ ) {
            EvtGenReport( EVTGEN_INFO, "EvtGen" ) << lambdaB2[i] << endl;
        }

        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << "Helicity states of particle C:" << endl;
        for ( int i = 0; i < nC; i++ ) {
            EvtGenReport( EVTGEN_INFO, "EvtGen" ) << lambdaC2[i] << endl;
        }

        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << "Will now figure out the valid (M_LS) states:" << endl;
    }

    int Lmin = std::max( JA2 - JB2 - JC2,
                         std::max( JB2 - JA2 - JC2, JC2 - JA2 - JB2 ) );
    if ( Lmin < 0 )
        Lmin = 0;
    int Lmax = JA2 + JB2 + JC2;

    int nPartialWaveAmp = 0;

    int nL[50];
    int nS[50];

    for ( int L = Lmin; L <= Lmax; L += 2 ) {
        int Smin = abs( L - JA2 );
        if ( Smin < abs( JB2 - JC2 ) )
            Smin = abs( JB2 - JC2 );
        int Smax = L + JA2;
        if ( Smax > abs( JB2 + JC2 ) )
            Smax = abs( JB2 + JC2 );
        for ( int S = Smin; S <= Smax; S += 2 ) {
            nL[nPartialWaveAmp] = L;
            nS[nPartialWaveAmp] = S;

            nPartialWaveAmp++;
            if ( verbose() ) {
                EvtGenReport( EVTGEN_INFO, "EvtGen" )
                    << "M[" << L << "][" << S << "]" << endl;
            }
        }
    }

    checkNArg( nPartialWaveAmp * 2 );

    int argcounter = 0;

    EvtComplex M[50];

    double partampsqtot = 0.0;

    for ( int i = 0; i < nPartialWaveAmp; i++ ) {
        M[i] = getArg( argcounter ) *
               exp( EvtComplex( 0.0, getArg( argcounter + 1 ) ) );
        argcounter += 2;
        partampsqtot += abs2( M[i] );
        if ( verbose() ) {
            EvtGenReport( EVTGEN_INFO, "EvtGen" )
                << "M[" << nL[i] << "][" << nS[i] << "]=" << M[i] << endl;
        }
    }

    //Now calculate the helicity amplitudes

    double helampsqtot = 0.0;

    for ( int ib = 0; ib < nB; ib++ ) {
        for ( int ic = 0; ic < nC; ic++ ) {
            HBC[ib][ic] = 0.0;
            if ( abs( lambdaB2[ib] - lambdaC2[ic] ) <= JA2 ) {
                for ( int i = 0; i < nPartialWaveAmp; i++ ) {
                    int L = nL[i];
                    int S = nS[i];
                    int lambda2 = lambdaB2[ib];
                    int lambda3 = lambdaC2[ic];
                    int s1 = JA2;
                    int s2 = JB2;
                    int s3 = JC2;
                    int m1 = lambda2 - lambda3;
                    EvtCGCoefSingle c1( s2, s3 );
                    EvtCGCoefSingle c2( L, S );

                    if ( verbose() ) {
                        EvtGenReport( EVTGEN_INFO, "EvtGen" )
                            << "s2,lambda2:" << s2 << " " << lambda2 << endl;
                    }
                    //fkw changes to satisfy KCC
                    double fkwTmp = ( L + 1.0 ) / ( s1 + 1.0 );

                    if ( S >= abs( m1 ) ) {
                        EvtComplex tmp = sqrt( fkwTmp ) *
                                         c1.coef( S, m1, s2, s3, lambda2,
                                                  -lambda3 ) *
                                         c2.coef( s1, m1, L, S, 0, m1 ) * M[i];
                        HBC[ib][ic] += tmp;
                    }
                }
                if ( verbose() ) {
                    EvtGenReport( EVTGEN_INFO, "EvtGen" )
                        << "HBC[" << ib << "][" << ic << "]=" << HBC[ib][ic]
                        << endl;
                }
            }
            helampsqtot += abs2( HBC[ib][ic] );
        }
    }

    if ( fabs( helampsqtot - partampsqtot ) / ( helampsqtot + partampsqtot ) >
         1e-6 ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "In EvtPartWave for decay " << EvtPDL::name( getParentId() )
            << " -> " << EvtPDL::name( getDaug( 0 ) ) << " "
            << EvtPDL::name( getDaug( 1 ) ) << std::endl;
        EvtGenReport( EVTGEN_ERROR, "EvtGen" ) << "With arguments: " << std::endl;
        for ( int i = 0; i * 2 < getNArg(); i++ ) {
            EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                << "M(" << nL[i] << "," << nS[i]
                << ")="
                //				 <<getArg(2*i)<<" "<<getArg(2*i+1)<<std::endl;
                << M[i] << std::endl;
        }
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "The total probability in the partwave basis is: " << partampsqtot
            << std::endl;
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "The total probability in the helamp basis is: " << helampsqtot
            << std::endl;
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Most likely this is because the specified partwave amplitudes "
            << std::endl;
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "project onto unphysical helicities of photons or neutrinos. "
            << std::endl;
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Seriously consider if your specified amplitudes are correct. "
            << std::endl;
    }

    m_evalHelAmp = std::make_unique<EvtEvalHelAmp>( getParentId(), getDaug( 0 ),
                                                    getDaug( 1 ), HBC );
}

void EvtPartWave::initProbMax()
{
    double maxprob = m_evalHelAmp->probMax();

    if ( verbose() ) {
        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << "Calculated probmax" << maxprob << endl;
    }

    setProbMax( maxprob );
}

void EvtPartWave::decay( EvtParticle* p )
{
    //first generate simple phase space
    p->initializePhaseSpace( getNDaug(), getDaugs() );

    m_evalHelAmp->evalAmp( p, m_amp2 );

    return;
}

void EvtPartWave::fillHelicity( int* lambda2, int n, int J2 )
{
    //photon is special case!
    if ( n == 2 && J2 == 2 ) {
        lambda2[0] = 2;
        lambda2[1] = -2;
        return;
    }

    assert( n == J2 + 1 );

    for ( int i = 0; i < n; i++ ) {
        lambda2[i] = n - i * 2 - 1;
    }

    return;
}
