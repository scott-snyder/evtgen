
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

#ifndef EVT_INTERVAL_DECAY_AMP
#define EVT_INTERVAL_DECAY_AMP

#define VERBOSE true
#include "EvtGenBase/EvtAmpFactory.hh"
#include "EvtGenBase/EvtAmpPdf.hh"
#include "EvtGenBase/EvtCPUtil.hh"
#include "EvtGenBase/EvtCyclic3.hh"
#include "EvtGenBase/EvtDecayAmp.hh"
#include "EvtGenBase/EvtMacros.hh"
#include "EvtGenBase/EvtMultiChannelParser.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPdf.hh"
#include "EvtGenBase/EvtReport.hh"

#include <iostream>
#include <string>
#include <vector>

// Decay model that uses the "amplitude on an interval"
// templatization

template <class T>
class EvtIntervalDecayAmp : public EvtDecayAmp {
  public:
    EvtIntervalDecayAmp() : m_probMax( 0. ), m_nScan( 0 ), m_fact( nullptr ) {}

    EvtIntervalDecayAmp( const EvtIntervalDecayAmp<T>& other ) :
        m_probMax( other.m_probMax ), m_nScan( other.m_nScan ), COPY_PTR( m_fact )
    {
    }

    virtual ~EvtIntervalDecayAmp() { delete m_fact; }

    // Initialize model

    void init() override
    {
        // Collect model parameters and parse them

        std::vector<std::string> args;
        int i;
        for ( i = 0; i < getNArg(); i++ )
            args.push_back( getArgStr( i ) );
        EvtMultiChannelParser parser;
        parser.parse( args );

        // Create factory and interval

        if ( VERBOSE )
            EvtGenReport( EVTGEN_INFO, "EvtGen" )
                << "Create factory and interval" << std::endl;
        m_fact = createFactory( parser );

        // Maximum PDF value over the Dalitz plot can be specified, or a scan
        // can be performed.

        m_probMax = parser.pdfMax();
        m_nScan = parser.nScan();
        if ( VERBOSE )
            EvtGenReport( EVTGEN_INFO, "EvtGen" )
                << "Pdf maximum " << m_probMax << std::endl;
        if ( VERBOSE )
            EvtGenReport( EVTGEN_INFO, "EvtGen" )
                << "Scan number " << m_nScan << std::endl;
    }

    void initProbMax() override
    {
        if ( 0 == m_nScan ) {
            if ( m_probMax > 0 )
                setProbMax( m_probMax );
            else
                assert( 0 );
        } else {
            double factor = 1.2;    // increase maximum probability by 20%
            EvtAmpPdf<T> pdf( *m_fact->getAmp() );
            EvtPdfSum<T>* pc = m_fact->getPC();
            EvtPdfDiv<T> pdfdiv( pdf, *pc );
            printf( "Sampling %d points to find maximum\n", m_nScan );
            EvtPdfMax<T> x = pdfdiv.findMax( *pc, m_nScan );
            m_probMax = factor * x.value();
            printf( "Found maximum %f\n", x.value() );
            printf( "Increase to   %f\n", m_probMax );
            setProbMax( m_probMax );
        }
    }

    void decay( EvtParticle* p ) override
    {
        // Set things up in most general way

        static EvtId B0 = EvtPDL::getId( "B0" );
        static EvtId B0B = EvtPDL::getId( "anti-B0" );
        double t;
        EvtId other_b;
        EvtComplex ampl( 0., 0. );

        // Sample using pole-compensator pdf

        EvtPdfSum<T>* pc = getPC();
        m_x = pc->randomPoint();

        if ( m_fact->isCPModel() ) {
            // Time-dependent Dalitz plot changes
            // Dec 2005 (ddujmic@slac.stanford.edu)

            EvtComplex A = m_fact->getAmp()->evaluate( m_x );
            EvtComplex Abar = m_fact->getAmpConj()->evaluate( m_x );

            EvtCPUtil::getInstance()->OtherB( p, t, other_b );

            double dm = m_fact->dm();
            double mixAmpli = m_fact->mixAmpli();
            double mixPhase = m_fact->mixPhase();
            EvtComplex qoverp( cos( mixPhase ) * mixAmpli,
                               sin( mixPhase ) * mixAmpli );
            EvtComplex poverq( cos( mixPhase ) / mixAmpli,
                               -sin( mixPhase ) / mixAmpli );

            if ( other_b == B0B )
                ampl = A * cos( dm * t / ( 2 * EvtConst::c ) ) +
                       EvtComplex( 0., 1. ) * Abar *
                           sin( dm * t / ( 2 * EvtConst::c ) ) * qoverp;
            if ( other_b == B0 )
                ampl = Abar * cos( dm * t / ( 2 * EvtConst::c ) ) +
                       EvtComplex( 0., 1. ) * A *
                           sin( dm * t / ( 2 * EvtConst::c ) ) * poverq;

        } else {
            ampl = amplNonCP( m_x );
        }

        // Pole-compensate

        double comp = sqrt( pc->evaluate( m_x ) );
        assert( comp > 0 );
        vertex( ampl / comp );

        // Now generate random angles, rotate and setup
        // the daughters

        std::vector<EvtVector4R> v = initDaughters( m_x );

        size_t N = p->getNDaug();
        if ( v.size() != N ) {
            EvtGenReport( EVTGEN_INFO, "EvtGen" )
                << "Number of daughters " << N << std::endl;
            EvtGenReport( EVTGEN_INFO, "EvtGen" )
                << "Momentum vector size " << v.size() << std::endl;
            assert( 0 );
        }

        for ( size_t i = 0; i < N; i++ ) {
            p->getDaug( i )->init( getDaugs()[i], v[i] );
        }
    }

    virtual EvtAmpFactory<T>* createFactory(
        const EvtMultiChannelParser& parser ) = 0;
    virtual std::vector<EvtVector4R> initDaughters( const T& p ) const = 0;

    // provide access to the decay point and to the amplitude of any decay point.
    // this is used by EvtBtoKD3P:
    const T& x() const { return m_x; }
    EvtComplex amplNonCP( const T& x )
    {
        return m_fact->getAmp()->evaluate( x );
    }
    EvtPdfSum<T>* getPC() { return m_fact->getPC(); }

  protected:
    double m_probMax;    // Maximum probability
    int m_nScan;         // Number of points for max prob DP scan
    T m_x;               // Decay point

    EvtAmpFactory<T>* m_fact;    // factory
};

#endif
