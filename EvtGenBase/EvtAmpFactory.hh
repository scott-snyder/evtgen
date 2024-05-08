
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

#ifndef EVT_AMP_FACTORY_HH
#define EVT_AMP_FACTORY_HH

#include "EvtGenBase/EvtAmpPdf.hh"
#include "EvtGenBase/EvtAmplitudeSum.hh"
#include "EvtGenBase/EvtMacros.hh"
#include "EvtGenBase/EvtMultiChannelParser.hh"
#include "EvtGenBase/EvtPdfMax.hh"
#include "EvtGenBase/EvtPdfSum.hh"

#include <stdio.h>
#include <string>
#include <vector>

// Abstract amplitude factory parameterized by a vector of
// strings. Derived classes construct the amplitude, and PDFs for sampling
// points.

template <class T>
class EvtAmpFactory {
  public:
    EvtAmpFactory() = default;

  protected:
    EvtAmpFactory( EvtAmpFactory<T>&& ) = default;
    EvtAmpFactory( const EvtAmpFactory<T>& other ) :
        m_amp( other.m_amp ? other.m_amp->clone() : nullptr ),
        m_ampConj( other.m_ampConj ? other.m_ampConj->clone() : nullptr ),
        m_pc( other.m_pc ? other.m_pc->clone() : nullptr ),
        m_names( other.m_names ),
        m_dm( other.m_dm ),
        m_mixPhase( other.m_mixPhase ),
        m_verbose( other.m_verbose )
    {
    }

  public:
    virtual ~EvtAmpFactory() = default;

    virtual EvtAmpFactory<T>* clone() const = 0;

    virtual void build( const EvtMultiChannelParser& parser, int nItg )
    {
        m_amp = std::make_unique<EvtAmplitudeSum<T>>();
        m_ampConj = std::make_unique<EvtAmplitudeSum<T>>();
        m_pc = std::make_unique<EvtPdfSum<T>>();
        m_dm = parser.dm();
        m_mixAmpli = parser.mixAmpli();
        m_mixPhase = parser.mixPhase();

        printf( "Amplitude with %d terms\n", parser.getNAmp() );
        int i;
        for ( i = 0; i < parser.getNAmp(); i++ ) {
            std::vector<std::string> v = parser.amp( i );
            EvtComplex c = parser.ampCoef( i );
            processAmp( c, v );
        }

        printf( "Conj. amplitude with %d terms\n", parser.getNAmpConj() );
        for ( i = 0; i < parser.getNAmpConj(); i++ ) {
            std::vector<std::string> v = parser.ampConj( i );
            EvtComplex c = parser.ampConjCoef( i );
            processAmp( c, v, true );
        }

        printf( "Calculating pole compensator integrals %d steps\n", nItg );
        if ( nItg > 0 )
            m_pc->getItg( nItg );

        printf( "End build\n" );
    }

    virtual void processAmp( EvtComplex c, std::vector<std::string> v,
                             bool conj = false ) = 0;

    inline bool isCPModel() const
    {
        return ( m_ampConj->nTerms() > 0 ? true : false );
    }
    inline double dm() const { return m_dm; }
    inline double mixAmpli() const { return m_mixAmpli; }
    inline double mixPhase() const { return m_mixPhase; }

    void setVerbose() { m_verbose = true; }

    EvtAmplitudeSum<T>* getAmp() const { return m_amp.get(); }
    EvtAmplitudeSum<T>* getAmpConj() const { return m_ampConj.get(); }
    EvtPdfSum<T>* getPC() const { return m_pc.get(); }
    EvtAmplitude<T>* getAmp( int i ) const { return m_amp->getTerm( i ); }
    EvtPdf<T>* getPC( int i ) const { return m_pc->getPdf( i ); }
    const char* compName( int i ) const { return m_names[i].c_str(); }

    EvtComplex getCoeff( int i ) const { return m_amp->c( i ); }

    double getTermCoeff( int i ) const { return abs2( m_amp->c( i ) ); }
    double getTermCoeff( int type, int i, int j ) const
    {
        switch ( type ) {
            case 0:
                return 2 * real( m_amp->c( i ) * conj( m_amp->c( j ) ) );    //posre
            case 1:
                return -2 *
                       real( m_amp->c( i ) * conj( m_amp->c( j ) ) );    //negre
            case 2:
                return -2 *
                       imag( m_amp->c( i ) * conj( m_amp->c( j ) ) );    //posim
            case 3:
                return 2 * imag( m_amp->c( i ) * conj( m_amp->c( j ) ) );    //negim
            default:
                assert( 0 );
        }
    }

  protected:
    std::unique_ptr<EvtAmplitudeSum<T>> m_amp;        // _owned_ amplitude
    std::unique_ptr<EvtAmplitudeSum<T>> m_ampConj;    // _owned_ conjugate amplitude
    std::unique_ptr<EvtPdfSum<T>> m_pc;    // _owned_ pole compensator
    std::vector<std::string> m_names;      // names of partial amplitudes

    double m_dm = 0;      // Mass difference for conjugate amplitude
    double m_mixPhase;    // mixing phase
    double m_mixAmpli;    // cpv in mixing
    bool m_verbose = false;
};

#endif
