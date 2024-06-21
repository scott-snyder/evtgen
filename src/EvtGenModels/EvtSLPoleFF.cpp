
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

#include "EvtGenModels/EvtSLPoleFF.hh"

#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtReport.hh"

#include <math.h>
#include <stdlib.h>
#include <string>

EvtSLPoleFF::EvtSLPoleFF( int numarg, double* arglist )
{
    //arg - maybe ignore the last argument - if odd ... Sigh
    m_numSLPoleargs = numarg - ( numarg % 2 );
    for ( int i = 0; i < m_numSLPoleargs; i++ ) {
        m_SLPoleargs[i] = arglist[i];
    }

    return;
}

void EvtSLPoleFF::getscalarff( EvtId parent, EvtId, double t, double,
                               double* fpf, double* f0f )
{
    // Form factors have a general form, with parameters passed in
    // from the arguements.

    if ( m_numSLPoleargs != 8 ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Problem in EvtSLPoleFF::getscalarff\n";
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "wrong number of arguements!!!\n";
    }

    double mb = EvtPDL::getMeanMass( parent );
    double mb2 = mb * mb;

    double f0, af, bf, powf;

    f0 = m_SLPoleargs[0];
    af = m_SLPoleargs[1];
    bf = m_SLPoleargs[2];
    powf = m_SLPoleargs[3];
    *fpf = f0 /
           ( pow( 1.0 + ( af * t / mb2 ) + ( bf * ( ( t / mb2 ) * ( t / mb2 ) ) ),
                  powf ) );

    f0 = m_SLPoleargs[4];
    af = m_SLPoleargs[5];
    bf = m_SLPoleargs[6];
    powf = m_SLPoleargs[7];

    *f0f = f0 /
           ( pow( 1.0 + ( af * t / mb2 ) + ( bf * ( ( t / mb2 ) * ( t / mb2 ) ) ),
                  powf ) );

    return;
}

void EvtSLPoleFF::getvectorff( EvtId parent, EvtId, double t, double,
                               double* a1f, double* a2f, double* vf, double* a0f )
{
    if ( m_numSLPoleargs != 16 ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Problem in EvtSLPoleFF::getvectorff\n";
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "wrong number of arguements!!!\n";
        EvtGenReport( EVTGEN_ERROR, "EvtGen" ) << m_numSLPoleargs << "\n";
    }

    double mb = EvtPDL::getMeanMass( parent );
    double mb2 = mb * mb;

    double f0, af, bf, powf;

    f0 = m_SLPoleargs[0];
    af = m_SLPoleargs[1];
    bf = m_SLPoleargs[2];
    powf = m_SLPoleargs[3];
    *a1f = f0 /
           ( pow( 1.0 + ( af * t / mb2 ) + ( bf * ( ( t / mb2 ) * ( t / mb2 ) ) ),
                  powf ) );

    f0 = m_SLPoleargs[4];
    af = m_SLPoleargs[5];
    bf = m_SLPoleargs[6];
    powf = m_SLPoleargs[7];

    *a2f = f0 /
           ( pow( 1.0 + ( af * t / mb2 ) + ( bf * ( ( t / mb2 ) * ( t / mb2 ) ) ),
                  powf ) );

    f0 = m_SLPoleargs[8];
    af = m_SLPoleargs[9];
    bf = m_SLPoleargs[10];
    powf = m_SLPoleargs[11];

    *vf = f0 /
          ( pow( 1.0 + ( af * t / mb2 ) + ( bf * ( ( t / mb2 ) * ( t / mb2 ) ) ),
                 powf ) );

    f0 = m_SLPoleargs[12];
    af = m_SLPoleargs[13];
    bf = m_SLPoleargs[14];
    powf = m_SLPoleargs[15];

    *a0f = f0 /
           ( pow( 1.0 + ( af * t / mb2 ) + ( bf * ( ( t / mb2 ) * ( t / mb2 ) ) ),
                  powf ) );
    return;
}

void EvtSLPoleFF::gettensorff( EvtId parent, EvtId, double t, double,
                               double* hf, double* kf, double* bpf, double* bmf )
{
    if ( m_numSLPoleargs != 16 ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Problem in EvtSLPoleFF::gettensorff\n";
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "wrong number of arguements!!!\n";
    }

    double mb = EvtPDL::getMeanMass( parent );
    double mb2 = mb * mb;

    double f0, af, bf, powf;

    f0 = m_SLPoleargs[0];
    af = m_SLPoleargs[1];
    bf = m_SLPoleargs[2];
    powf = m_SLPoleargs[3];
    *hf = f0 /
          ( pow( 1.0 + ( af * t / mb2 ) + ( bf * ( ( t / mb2 ) * ( t / mb2 ) ) ),
                 powf ) );

    f0 = m_SLPoleargs[4];
    af = m_SLPoleargs[5];
    bf = m_SLPoleargs[6];
    powf = m_SLPoleargs[7];

    *kf = f0 /
          ( pow( 1.0 + ( af * t / mb2 ) + ( bf * ( ( t / mb2 ) * ( t / mb2 ) ) ),
                 powf ) );

    f0 = m_SLPoleargs[8];
    af = m_SLPoleargs[9];
    bf = m_SLPoleargs[10];
    powf = m_SLPoleargs[11];

    *bpf = f0 /
           ( pow( 1.0 + ( af * t / mb2 ) + ( bf * ( ( t / mb2 ) * ( t / mb2 ) ) ),
                  powf ) );

    f0 = m_SLPoleargs[12];
    af = m_SLPoleargs[13];
    bf = m_SLPoleargs[14];
    powf = m_SLPoleargs[15];

    *bmf = f0 /
           ( pow( 1.0 + ( af * t / mb2 ) + ( bf * ( ( t / mb2 ) * ( t / mb2 ) ) ),
                  powf ) );
    return;
}

void EvtSLPoleFF::getbaryonff( EvtId, EvtId, double, double, double*, double*,
                               double*, double* )
{
    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "Not implemented :getbaryonff in EvtSLPoleFF.\n";
    ::abort();
}

void EvtSLPoleFF::getdiracff( EvtId, EvtId, double, double, double*, double*,
                              double*, double*, double*, double* )
{
    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "Not implemented :getdiracff in EvtSLPoleFF.\n";
    ::abort();
}

void EvtSLPoleFF::getraritaff( EvtId, EvtId, double, double, double*, double*,
                               double*, double*, double*, double*, double*,
                               double* )
{
    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "Not implemented :getraritaff in EvtSLPoleFF.\n";
    ::abort();
}
