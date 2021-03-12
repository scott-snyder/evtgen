
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

#include "EvtGenBase/EvtEvalHelAmp.hh"

#include "EvtGenBase/EvtAmp.hh"
#include "EvtGenBase/EvtConst.hh"
#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtVector4C.hh"
#include "EvtGenBase/EvtdFunction.hh"

#include <cstdlib>
#include <iostream>

using std::endl;

EvtEvalHelAmp::~EvtEvalHelAmp()
{
    //deallocate memory
    delete[] m_lambdaA2;
    delete[] m_lambdaB2;
    delete[] m_lambdaC2;

    int ia, ib, ic;
    for ( ib = 0; ib < m_nB; ib++ ) {
        delete[] m_HBC[ib];
    }

    delete[] m_HBC;

    for ( ia = 0; ia < m_nA; ia++ ) {
        delete[] m_RA[ia];
    }
    delete[] m_RA;

    for ( ib = 0; ib < m_nB; ib++ ) {
        delete[] m_RB[ib];
    }
    delete[] m_RB;

    for ( ic = 0; ic < m_nC; ic++ ) {
        delete[] m_RC[ic];
    }
    delete[] m_RC;

    for ( ia = 0; ia < m_nA; ia++ ) {
        for ( ib = 0; ib < m_nB; ib++ ) {
            delete[] m_amp[ia][ib];
            delete[] m_amp1[ia][ib];
            delete[] m_amp3[ia][ib];
        }
        delete[] m_amp[ia];
        delete[] m_amp1[ia];
        delete[] m_amp3[ia];
    }

    delete[] m_amp;
    delete[] m_amp1;
    delete[] m_amp3;
}

EvtEvalHelAmp::EvtEvalHelAmp( EvtId idA, EvtId idB, EvtId idC,
                              EvtComplexPtrPtr HBC )
{
    EvtSpinType::spintype typeA = EvtPDL::getSpinType( idA );
    EvtSpinType::spintype typeB = EvtPDL::getSpinType( idB );
    EvtSpinType::spintype typeC = EvtPDL::getSpinType( idC );

    //find out how many states each particle have
    m_nA = EvtSpinType::getSpinStates( typeA );
    m_nB = EvtSpinType::getSpinStates( typeB );
    m_nC = EvtSpinType::getSpinStates( typeC );

    //find out what 2 times the spin is
    m_JA2 = EvtSpinType::getSpin2( typeA );
    m_JB2 = EvtSpinType::getSpin2( typeB );
    m_JC2 = EvtSpinType::getSpin2( typeC );

    //allocate memory
    m_lambdaA2 = new int[m_nA];
    m_lambdaB2 = new int[m_nB];
    m_lambdaC2 = new int[m_nC];

    m_HBC = new EvtComplexPtr[m_nB];
    int ia, ib, ic;
    for ( ib = 0; ib < m_nB; ib++ ) {
        m_HBC[ib] = new EvtComplex[m_nC];
    }

    m_RA = new EvtComplexPtr[m_nA];
    for ( ia = 0; ia < m_nA; ia++ ) {
        m_RA[ia] = new EvtComplex[m_nA];
    }
    m_RB = new EvtComplexPtr[m_nB];
    for ( ib = 0; ib < m_nB; ib++ ) {
        m_RB[ib] = new EvtComplex[m_nB];
    }
    m_RC = new EvtComplexPtr[m_nC];
    for ( ic = 0; ic < m_nC; ic++ ) {
        m_RC[ic] = new EvtComplex[m_nC];
    }

    m_amp = new EvtComplexPtrPtr[m_nA];
    m_amp1 = new EvtComplexPtrPtr[m_nA];
    m_amp3 = new EvtComplexPtrPtr[m_nA];
    for ( ia = 0; ia < m_nA; ia++ ) {
        m_amp[ia] = new EvtComplexPtr[m_nB];
        m_amp1[ia] = new EvtComplexPtr[m_nB];
        m_amp3[ia] = new EvtComplexPtr[m_nB];
        for ( ib = 0; ib < m_nB; ib++ ) {
            m_amp[ia][ib] = new EvtComplex[m_nC];
            m_amp1[ia][ib] = new EvtComplex[m_nC];
            m_amp3[ia][ib] = new EvtComplex[m_nC];
        }
    }

    //find the allowed helicities (actually 2*times the helicity!)

    fillHelicity( m_lambdaA2, m_nA, m_JA2, idA );
    fillHelicity( m_lambdaB2, m_nB, m_JB2, idB );
    fillHelicity( m_lambdaC2, m_nC, m_JC2, idC );

    for ( ib = 0; ib < m_nB; ib++ ) {
        for ( ic = 0; ic < m_nC; ic++ ) {
            m_HBC[ib][ic] = HBC[ib][ic];
        }
    }
}

double EvtEvalHelAmp::probMax()
{
    double c = 1.0 / sqrt( 4 * EvtConst::pi / ( m_JA2 + 1 ) );

    int ia, ib, ic;

    double theta;
    int itheta;

    double maxprob = 0.0;

    for ( itheta = -10; itheta <= 10; itheta++ ) {
        theta = acos( 0.099999 * itheta );
        for ( ia = 0; ia < m_nA; ia++ ) {
            double prob = 0.0;
            for ( ib = 0; ib < m_nB; ib++ ) {
                for ( ic = 0; ic < m_nC; ic++ ) {
                    m_amp[ia][ib][ic] = 0.0;
                    if ( abs( m_lambdaB2[ib] - m_lambdaC2[ic] ) <= m_JA2 ) {
                        m_amp[ia][ib][ic] =
                            c * m_HBC[ib][ic] *
                            EvtdFunction::d( m_JA2, m_lambdaA2[ia],
                                             m_lambdaB2[ib] - m_lambdaC2[ic],
                                             theta );
                        prob += real( m_amp[ia][ib][ic] *
                                      conj( m_amp[ia][ib][ic] ) );
                    }
                }
            }

            prob *= sqrt( 1.0 * m_nA );

            if ( prob > maxprob )
                maxprob = prob;
        }
    }

    return maxprob;
}

void EvtEvalHelAmp::evalAmp( EvtParticle* p, EvtAmp& amp )
{
    //find theta and phi of the first daughter

    EvtVector4R pB = p->getDaug( 0 )->getP4();

    double theta = acos( pB.get( 3 ) / pB.d3mag() );
    double phi = atan2( pB.get( 2 ), pB.get( 1 ) );

    double c = sqrt( ( m_JA2 + 1 ) / ( 4 * EvtConst::pi ) );

    int ia, ib, ic;

    double prob1 = 0.0;

    for ( ia = 0; ia < m_nA; ia++ ) {
        for ( ib = 0; ib < m_nB; ib++ ) {
            for ( ic = 0; ic < m_nC; ic++ ) {
                m_amp[ia][ib][ic] = 0.0;
                if ( abs( m_lambdaB2[ib] - m_lambdaC2[ic] ) <= m_JA2 ) {
                    double dfun = EvtdFunction::d( m_JA2, m_lambdaA2[ia],
                                                   m_lambdaB2[ib] - m_lambdaC2[ic],
                                                   theta );

                    m_amp[ia][ib][ic] =
                        c * m_HBC[ib][ic] *
                        exp( EvtComplex( 0.0,
                                         phi * 0.5 *
                                             ( m_lambdaA2[ia] - m_lambdaB2[ib] +
                                               m_lambdaC2[ic] ) ) ) *
                        dfun;
                }
                prob1 += real( m_amp[ia][ib][ic] * conj( m_amp[ia][ib][ic] ) );
            }
        }
    }

    setUpRotationMatrices( p, theta, phi );

    applyRotationMatrices();

    double prob2 = 0.0;

    for ( ia = 0; ia < m_nA; ia++ ) {
        for ( ib = 0; ib < m_nB; ib++ ) {
            for ( ic = 0; ic < m_nC; ic++ ) {
                prob2 += real( m_amp[ia][ib][ic] * conj( m_amp[ia][ib][ic] ) );
                if ( m_nA == 1 ) {
                    if ( m_nB == 1 ) {
                        if ( m_nC == 1 ) {
                            amp.vertex( m_amp[ia][ib][ic] );
                        } else {
                            amp.vertex( ic, m_amp[ia][ib][ic] );
                        }
                    } else {
                        if ( m_nC == 1 ) {
                            amp.vertex( ib, m_amp[ia][ib][ic] );
                        } else {
                            amp.vertex( ib, ic, m_amp[ia][ib][ic] );
                        }
                    }
                } else {
                    if ( m_nB == 1 ) {
                        if ( m_nC == 1 ) {
                            amp.vertex( ia, m_amp[ia][ib][ic] );
                        } else {
                            amp.vertex( ia, ic, m_amp[ia][ib][ic] );
                        }
                    } else {
                        if ( m_nC == 1 ) {
                            amp.vertex( ia, ib, m_amp[ia][ib][ic] );
                        } else {
                            amp.vertex( ia, ib, ic, m_amp[ia][ib][ic] );
                        }
                    }
                }
            }
        }
    }

    if ( fabs( prob1 - prob2 ) > 0.000001 * prob1 ) {
        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << "prob1,prob2:" << prob1 << " " << prob2 << endl;
        ::abort();
    }

    return;
}

void EvtEvalHelAmp::fillHelicity( int* lambda2, int n, int J2, EvtId id )
{
    int i;

    //photon is special case!
    if ( n == 2 && J2 == 2 ) {
        lambda2[0] = 2;
        lambda2[1] = -2;
        return;
    }

    //and so is the neutrino!
    if ( n == 1 && J2 == 1 ) {
        if ( EvtPDL::getStdHep( id ) > 0 ) {
            //particle i.e. lefthanded
            lambda2[0] = -1;
        } else {
            //anti particle i.e. righthanded
            lambda2[0] = 1;
        }
        return;
    }

    assert( n == J2 + 1 );

    for ( i = 0; i < n; i++ ) {
        lambda2[i] = n - i * 2 - 1;
    }

    return;
}

void EvtEvalHelAmp::setUpRotationMatrices( EvtParticle* p, double theta,
                                           double phi )
{
    switch ( m_JA2 ) {
        case 0:
        case 1:
        case 2:
        case 3:
        case 4:
        case 5:
        case 6:
        case 7:
        case 8:

        {
            EvtSpinDensity R = p->rotateToHelicityBasis();

            int i, j, n;

            n = R.getDim();

            assert( n == m_nA );

            for ( i = 0; i < n; i++ ) {
                for ( j = 0; j < n; j++ ) {
                    m_RA[i][j] = R.get( i, j );
                }
            }

        }

        break;

        default:
            EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                << "Spin2(m_JA2)=" << m_JA2 << " not supported!" << endl;
            ::abort();
    }

    switch ( m_JB2 ) {
        case 0:
        case 1:
        case 2:
        case 3:
        case 4:
        case 5:
        case 6:
        case 7:
        case 8:

        {
            int i, j, n;

            EvtSpinDensity R =
                p->getDaug( 0 )->rotateToHelicityBasis( phi, theta, -phi );

            n = R.getDim();

            assert( n == m_nB );

            for ( i = 0; i < n; i++ ) {
                for ( j = 0; j < n; j++ ) {
                    m_RB[i][j] = conj( R.get( i, j ) );
                }
            }

        }

        break;

        default:
            EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                << "Spin2(m_JB2)=" << m_JB2 << " not supported!" << endl;
            ::abort();
    }

    switch ( m_JC2 ) {
        case 0:
        case 1:
        case 2:
        case 3:
        case 4:
        case 5:
        case 6:
        case 7:
        case 8:

        {
            int i, j, n;

            EvtSpinDensity R = p->getDaug( 1 )->rotateToHelicityBasis(
                phi, EvtConst::pi + theta, phi - EvtConst::pi );

            n = R.getDim();

            assert( n == m_nC );

            for ( i = 0; i < n; i++ ) {
                for ( j = 0; j < n; j++ ) {
                    m_RC[i][j] = conj( R.get( i, j ) );
                }
            }

        }

        break;

        default:
            EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                << "Spin2(m_JC2)=" << m_JC2 << " not supported!" << endl;
            ::abort();
    }
}

void EvtEvalHelAmp::applyRotationMatrices()
{
    int ia, ib, ic, i;

    EvtComplex temp;

    for ( ia = 0; ia < m_nA; ia++ ) {
        for ( ib = 0; ib < m_nB; ib++ ) {
            for ( ic = 0; ic < m_nC; ic++ ) {
                temp = 0;
                for ( i = 0; i < m_nC; i++ ) {
                    temp += m_RC[i][ic] * m_amp[ia][ib][i];
                }
                m_amp1[ia][ib][ic] = temp;
            }
        }
    }

    for ( ia = 0; ia < m_nA; ia++ ) {
        for ( ic = 0; ic < m_nC; ic++ ) {
            for ( ib = 0; ib < m_nB; ib++ ) {
                temp = 0;
                for ( i = 0; i < m_nB; i++ ) {
                    temp += m_RB[i][ib] * m_amp1[ia][i][ic];
                }
                m_amp3[ia][ib][ic] = temp;
            }
        }
    }

    for ( ib = 0; ib < m_nB; ib++ ) {
        for ( ic = 0; ic < m_nC; ic++ ) {
            for ( ia = 0; ia < m_nA; ia++ ) {
                temp = 0;
                for ( i = 0; i < m_nA; i++ ) {
                    temp += m_RA[i][ia] * m_amp3[i][ib][ic];
                }
                m_amp[ia][ib][ic] = temp;
            }
        }
    }
}
