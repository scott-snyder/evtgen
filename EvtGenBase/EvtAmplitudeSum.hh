
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

#ifndef EVT_AMPLITUDE_SUM_HH
#define EVT_AMPLITUDE_SUM_HH

#include "EvtGenBase/EvtAmplitude.hh"

#include <assert.h>
#include <memory>
#include <vector>

template <class T>
class EvtAmplitudeSum : public EvtAmplitude<T> {
  public:
    EvtAmplitudeSum() {}
    EvtAmplitudeSum( const EvtAmplitudeSum<T>& other ) :
        EvtAmplitude<T>( other )
    {
        int i;
        for ( i = 0; i < other.nTerms(); i++ ) {
            EvtComplex c = other.c( i );
            m_c.push_back( c );
            EvtAmplitude<T>* amp = other.getTerm( i );
            assert( amp );
            EvtAmplitude<T>* amp1 = amp->clone();
            assert( amp1 );
            m_term.push_back( amp1 );
        }
    }

    virtual ~EvtAmplitudeSum()
    {
        for ( size_t i = 0; i < m_term.size(); i++ ) {
            delete m_term[i];
        }
    }

    EvtAmplitudeSum<T>* clone() const override
    {
        return new EvtAmplitudeSum<T>( *this );
    }

    void addTerm( EvtComplex c, const EvtAmplitude<T>& amp )
    {
        m_c.push_back( c );
        m_term.push_back( amp.clone() );
    }

    void addOwnedTerm( EvtComplex c, std::unique_ptr<EvtAmplitude<T>> amp )
    {
        assert( amp );
        m_c.push_back( c );
        m_term.push_back( amp.release() );
    }

    int nTerms() const { return m_term.size(); }    // number of terms

    void print() const
    {
        int N = nTerms();
        printf( "Amplitude has %d terms\n", N );
        int i;
        for ( i = 0; i < N; i++ ) {
            printf( "c%d = (%f,%f)\n", i, real( m_c[i] ), imag( m_c[i] ) );
            assert( m_term[i] );
        }
    }

    inline EvtComplex c( int i ) const { return m_c[i]; }
    inline EvtAmplitude<T>* getTerm( int i ) const { return m_term[i]; }

  protected:
    EvtComplex amplitude( const T& p ) const override
    {
        if ( m_term.size() == 0 )
            printf( "Warning: amplitude sum has zero terms\n" );

        EvtComplex value = 0.;

        for ( size_t i = 0; i < m_term.size(); i++ ) {
            value += m_c[i] * m_term[i]->evaluate( p );
        }
        return value;
    }

  private:
    std::vector<EvtComplex> m_c;             // coefficients
    std::vector<EvtAmplitude<T>*> m_term;    // pointers to amplitudes
};

#endif
