
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

#include "EvtGenModels/EvtBtoXsll.hh"

#include "EvtGenBase/EvtConst.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtReport.hh"

#include "EvtGenModels/EvtBtoXsllUtil.hh"
#include "EvtGenModels/EvtbTosllAmp.hh"

#include <stdlib.h>
using std::endl;

std::string EvtBtoXsll::getName() const
{
    return "BTOXSLL";
}

EvtDecayBase* EvtBtoXsll::clone() const
{
    return new EvtBtoXsll;
}

void EvtBtoXsll::init()
{
    // check that there are no arguments

    checkNArg( 0, 4, 5 );

    checkNDaug( 3 );

    // Check that the two leptons are the same type

    EvtId lepton1type = getDaug( 1 );
    EvtId lepton2type = getDaug( 2 );

    int etyp = 0;
    int mutyp = 0;
    int tautyp = 0;
    if ( lepton1type == EvtPDL::getId( "e+" ) ||
         lepton1type == EvtPDL::getId( "e-" ) ) {
        etyp++;
    }
    if ( lepton2type == EvtPDL::getId( "e+" ) ||
         lepton2type == EvtPDL::getId( "e-" ) ) {
        etyp++;
    }
    if ( lepton1type == EvtPDL::getId( "mu+" ) ||
         lepton1type == EvtPDL::getId( "mu-" ) ) {
        mutyp++;
    }
    if ( lepton2type == EvtPDL::getId( "mu+" ) ||
         lepton2type == EvtPDL::getId( "mu-" ) ) {
        mutyp++;
    }
    if ( lepton1type == EvtPDL::getId( "tau+" ) ||
         lepton1type == EvtPDL::getId( "tau-" ) ) {
        tautyp++;
    }
    if ( lepton2type == EvtPDL::getId( "tau+" ) ||
         lepton2type == EvtPDL::getId( "tau-" ) ) {
        tautyp++;
    }

    if ( etyp != 2 && mutyp != 2 && tautyp != 2 ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Expect two leptons of the same type in EvtBtoXsll.cc\n";
        ::abort();
    }

    // Check that the second and third entries are leptons with positive
    // and negative charge, respectively

    int lpos = 0;
    int lneg = 0;
    if ( lepton1type == EvtPDL::getId( "e+" ) ||
         lepton1type == EvtPDL::getId( "mu+" ) ||
         lepton1type == EvtPDL::getId( "tau+" ) ) {
        lpos++;
    }
    if ( lepton2type == EvtPDL::getId( "e-" ) ||
         lepton2type == EvtPDL::getId( "mu-" ) ||
         lepton2type == EvtPDL::getId( "tau-" ) ) {
        lneg++;
    }

    if ( lpos != 1 || lneg != 1 ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Expect 2nd and 3rd particles to be positive and negative leptons in EvtBtoXsll.cc\n";
        ::abort();
    }

    m_mb = 4.8;
    m_ms = 0.2;
    m_mq = 0.;
    m_pf = 0.41;
    m_mxmin = 1.1;
    if ( getNArg() == 4 ) {
        // b-quark mass
        m_mb = getArg( 0 );
        // s-quark mass
        m_ms = getArg( 1 );
        // spectator quark mass
        m_mq = getArg( 2 );
        // Fermi motion parameter
        m_pf = getArg( 3 );
    }
    if ( getNArg() == 5 ) {
        m_mxmin = getArg( 4 );
    }

    m_calcprob = std::make_unique<EvtBtoXsllUtil>();

    double ml = EvtPDL::getMeanMass( getDaug( 1 ) );

    // determine the maximum probability density from dGdsProb

    int i, j;
    int nsteps = 100;
    double s = 0.0;
    double smin = 4.0 * ml * ml;
    double smax = ( m_mb - m_ms ) * ( m_mb - m_ms );
    double probMax = -10000.0;
    double sProbMax = -10.0;
    double uProbMax = -10.0;

    for ( i = 0; i < nsteps; i++ ) {
        s = smin + ( i + 0.002 ) * ( smax - smin ) / (double)nsteps;
        double prob = m_calcprob->dGdsProb( m_mb, m_ms, ml, s );
        if ( prob > probMax ) {
            sProbMax = s;
            probMax = prob;
        }
    }

    m_dGdsProbMax = probMax;

    if ( verbose() ) {
        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << "dGdsProbMax = " << probMax << " for s = " << sProbMax << endl;
    }

    // determine the maximum probability density from dGdsdupProb

    probMax = -10000.0;
    sProbMax = -10.0;

    for ( i = 0; i < nsteps; i++ ) {
        s = smin + ( i + 0.002 ) * ( smax - smin ) / (double)nsteps;
        double umax = sqrt( ( s - ( m_mb + m_ms ) * ( m_mb + m_ms ) ) *
                            ( s - ( m_mb - m_ms ) * ( m_mb - m_ms ) ) );
        for ( j = 0; j < nsteps; j++ ) {
            double u = -umax + ( j + 0.002 ) * ( 2.0 * umax ) / (double)nsteps;
            double prob = m_calcprob->dGdsdupProb( m_mb, m_ms, ml, s, u );
            if ( prob > probMax ) {
                sProbMax = s;
                uProbMax = u;
                probMax = prob;
            }
        }
    }

    m_dGdsdupProbMax = 2.0 * probMax;

    if ( verbose() ) {
        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << "dGdsdupProbMax = " << probMax << " for s = " << sProbMax
            << " and u = " << uProbMax << endl;
    }
}

void EvtBtoXsll::initProbMax()
{
    noProbMax();
}

void EvtBtoXsll::decay( EvtParticle* p )
{
    p->makeDaughters( getNDaug(), getDaugs() );

    EvtParticle* xhadron = p->getDaug( 0 );
    EvtParticle* leptonp = p->getDaug( 1 );
    EvtParticle* leptonn = p->getDaug( 2 );

    double mass[3];

    findMasses( p, getNDaug(), getDaugs(), mass );

    double mB = p->mass();
    double ml = mass[1];
    double pb( 0. );

    static thread_local int nmsg = 0;

    double xhadronMass = -999.0;

    EvtVector4R p4xhadron;
    EvtVector4R p4leptonp;
    EvtVector4R p4leptonn;

    // require the hadronic system has mass greater than that of a Kaon pion pair

    //  while (xhadronMass < 0.6333)
    // the above minimum value of K+pi mass appears to be too close
    // to threshold as far as JETSET is concerned
    // (JETSET gets caught in an infinite loop)
    // so we choose a lightly larger value for the threshold
    while ( xhadronMass < m_mxmin ) {
        // Apply Fermi motion and determine effective b-quark mass

        // Old BaBar MC parameters
        //    double pf = 0.25;
        //    double ms = 0.2;
        //    double mq = 0.3;

        double mb = 0.0;

        double xbox, ybox;

        while ( mb <= 0.0 ) {
            pb = m_calcprob->FermiMomentum( m_pf );

            // effective b-quark mass
            mb = mB * mB + m_mq * m_mq - 2.0 * mB * sqrt( pb * pb + m_mq * m_mq );
            if ( mb > 0. && sqrt( mb ) - m_ms < 2.0 * ml )
                mb = -10.;
        }
        mb = sqrt( mb );

        //    cout << "b-quark momentum = " << pb << " mass = " <<  mb << endl;

        // generate a dilepton invariant mass

        double s = 0.0;
        double smin = 4.0 * ml * ml;
        double smax = ( mb - m_ms ) * ( mb - m_ms );

        while ( s == 0.0 ) {
            xbox = EvtRandom::Flat( smin, smax );
            ybox = EvtRandom::Flat( m_dGdsProbMax );
            double prob = m_calcprob->dGdsProb( mb, m_ms, ml, xbox );
            if ( !( prob >= 0.0 ) && !( prob <= 0.0 ) ) {
                //	EvtGenReport(EVTGEN_INFO,"EvtGen") << "nan from dGdsProb " << prob << " " << mb << " " << m_ms << " " << ml << " " << xbox << std::endl;
            }
            if ( ybox < prob )
                s = xbox;
        }

        //    cout << "dGdsProb(s) = " << m_calcprob->dGdsProb(mb, m_ms, ml, s)
        //         << " for s = " << s << endl;

        // two-body decay of b quark at rest into s quark and dilepton pair:
        // b -> s (ll)

        EvtVector4R p4sdilep[2];

        double msdilep[2];
        msdilep[0] = m_ms;
        msdilep[1] = sqrt( s );

        EvtGenKine::PhaseSpace( 2, msdilep, p4sdilep, mb );

        // generate dilepton decay with the expected asymmetry: (ll) -> l+ l-

        EvtVector4R p4ll[2];

        double mll[2];
        mll[0] = ml;
        mll[1] = ml;

        double tmp = 0.0;

        while ( tmp == 0.0 ) {
            // (ll) -> l+ l- decay in dilepton rest frame

            EvtGenKine::PhaseSpace( 2, mll, p4ll, msdilep[1] );

            // boost to b-quark rest frame

            p4ll[0] = boostTo( p4ll[0], p4sdilep[1] );
            p4ll[1] = boostTo( p4ll[1], p4sdilep[1] );

            // compute kinematical variable u

            EvtVector4R p4slp = p4sdilep[0] + p4ll[0];
            EvtVector4R p4sln = p4sdilep[0] + p4ll[1];

            double u = p4slp.mass2() - p4sln.mass2();

            ybox = EvtRandom::Flat( m_dGdsdupProbMax );

            double prob = m_calcprob->dGdsdupProb( mb, m_ms, ml, s, u );
            if ( !( prob >= 0.0 ) && !( prob <= 0.0 ) ) {
                EvtGenReport( EVTGEN_INFO, "EvtGen" )
                    << "nan from dGdsProb " << prob << " " << mb << " " << m_ms
                    << " " << ml << " " << s << " " << u << std::endl;
            }
            if ( prob > m_dGdsdupProbMax && nmsg < 20 ) {
                EvtGenReport( EVTGEN_INFO, "EvtGen" )
                    << "d2gdsdup GT d2gdsdup_max:" << prob << " "
                    << m_dGdsdupProbMax << " for s = " << s << " u = " << u
                    << " mb = " << mb << endl;
                nmsg++;
            }
            if ( ybox < prob ) {
                tmp = 1.0;
                //        cout << "dGdsdupProb(s) = " << prob
                //             << " for u = " << u << endl;
            }
        }

        // assign 4-momenta to valence quarks inside B meson in B rest frame

        double phi = EvtRandom::Flat( EvtConst::twoPi );
        double costh = EvtRandom::Flat( -1.0, 1.0 );
        double sinth = sqrt( 1.0 - costh * costh );

        // b-quark four-momentum in B meson rest frame

        EvtVector4R p4b( sqrt( mb * mb + pb * pb ), pb * sinth * sin( phi ),
                         pb * sinth * cos( phi ), pb * costh );

        // B meson in its rest frame
        //
        //    EvtVector4R p4B(mB, 0.0, 0.0, 0.0);
        //
        // boost B meson to b-quark rest frame
        //
        //    p4B = boostTo(p4B, p4b);
        //
        //    cout << " B meson mass in b-quark rest frame = " << p4B.mass() << endl;

        // boost s, l+ and l- to B meson rest frame

        //    EvtVector4R p4s = boostTo(p4sdilep[0], p4B);
        //    p4leptonp       = boostTo(p4ll[0],     p4B);
        //    p4leptonn       = boostTo(p4ll[1],     p4B);

        EvtVector4R p4s = boostTo( p4sdilep[0], p4b );
        p4leptonp = boostTo( p4ll[0], p4b );
        p4leptonn = boostTo( p4ll[1], p4b );

        // spectator quark in B meson rest frame

        EvtVector4R p4q( sqrt( pb * pb + m_mq * m_mq ), -p4b.get( 1 ),
                         -p4b.get( 2 ), -p4b.get( 3 ) );

        // hadron system in B meson rest frame

        p4xhadron = p4s + p4q;
        xhadronMass = p4xhadron.mass();
    }

    // initialize the decay products

    xhadron->init( getDaug( 0 ), p4xhadron );

    // For B-bar mesons (i.e. containing a b quark) we have the normal
    // order of leptons
    if ( p->getId() == EvtPDL::getId( "anti-B0" ) ||
         p->getId() == EvtPDL::getId( "B-" ) ) {
        leptonp->init( getDaug( 1 ), p4leptonp );
        leptonn->init( getDaug( 2 ), p4leptonn );
    }
    // For B mesons (i.e. containing a b-bar quark) we need to flip the
    // role of the positive and negative leptons in order to produce the
    // correct forward-backward asymmetry between the two leptons
    else {
        leptonp->init( getDaug( 1 ), p4leptonn );
        leptonn->init( getDaug( 2 ), p4leptonp );
    }

    return;
}
