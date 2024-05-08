
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

#ifndef EVT_PDF_HH
#define EVT_PDF_HH

#include "EvtGenBase/EvtMacros.hh"
#include "EvtGenBase/EvtPdfMax.hh"
#include "EvtGenBase/EvtPredGen.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtStreamInputIterator.hh"
#include "EvtGenBase/EvtValError.hh"

#include <assert.h>
#include <stdio.h>

/*
 *  All classes are templated on the point type T
 *
 * EvtPdf:
 *
 * Probability density function defined on an interval of phase-space.
 * Integral over the interval can be calculated by Monte Carlo integration.
 * Some (but not all) PDFs are analytic in the sense that they can be integrated
 * by numeric quadrature and distributions can be generated according to them.
 *
 * EvtPdfGen:
 *
 * Generator adaptor. Can be used to generate random points
 * distributed according to the PDF for analytic PDFs.
 *
 * EvtPdfPred:
 *
 * Predicate adaptor for PDFs. Can be used for generating random points distributed
 * according to the PDF for any PDF using rejection method. (See "Numerical Recipes").
 *
 * EvtPdfUnary:
 *
 * Adapter for generic algorithms. Evaluates the PDF and returns the value
 *
 * EvtPdfDiv:
 *
 * PDF obtained by division of one PDF by another. Because the two PDFs are
 * arbitrary this PDF is not analytic. When importance sampling is used the
 * original PDF is divided by the analytic comparison function. EvtPdfDiv is
 * used to represent the modified PDF.
 */

template <class T>
class EvtPdfPred;
template <class T>
class EvtPdfGen;

template <class T>
class EvtPdf {
  public:
    EvtPdf() {}
    EvtPdf( const EvtPdf& other ) : m_itg( other.m_itg ) {}
    virtual ~EvtPdf() {}
    virtual EvtPdf<T>* clone() const = 0;

    double evaluate( const T& p ) const
    {
        if ( p.isValid() )
            return pdf( p );
        else
            return 0.;
    }

    // Find PDF maximum. Points are sampled according to pc

    EvtPdfMax<T> findMax( const EvtPdf<T>& pc, int N );

    // Find generation efficiency.

    EvtValError findGenEff( const EvtPdf<T>& pc, int N, int nFindMax );

    // Analytic integration. Calls cascade down until an overridden
    // method is called.

    void setItg( EvtValError itg ) { m_itg = itg; }

    EvtValError getItg() const
    {
        if ( !m_itg.valueKnown() )
            m_itg = compute_integral();
        return m_itg;
    }
    EvtValError getItg( int N ) const
    {
        if ( !m_itg.valueKnown() )
            m_itg = compute_integral( N );
        return m_itg;
    }

    virtual EvtValError compute_integral() const
    {
        printf( "Analytic integration of PDF is not defined\n" );
        assert( 0 );
        return EvtValError{};
    }
    virtual EvtValError compute_integral( int ) const
    {
        return compute_integral();
    }

    //  Monte Carlo integration.

    EvtValError compute_mc_integral( const EvtPdf<T>& pc, int N );

    // Generation. Create predicate accept-reject generators.
    // nMax iterations will be used to find the maximum of the accept-reject predicate

    EvtPredGen<EvtPdfGen<T>, EvtPdfPred<T>> accRejGen( const EvtPdf<T>& pc,
                                                       int nMax,
                                                       double factor = 1. );

    virtual T randomPoint();

  protected:
    virtual double pdf( const T& ) const = 0;
    mutable EvtValError m_itg;
};

template <class T>
class EvtPdfGen {
  public:
    typedef T result_type;

    EvtPdfGen() : m_pdf( 0 ) {}
    EvtPdfGen( const EvtPdfGen<T>& other ) :
        m_pdf( other.m_pdf ? other.m_pdf->clone() : nullptr )
    {
    }
    EvtPdfGen( const EvtPdf<T>& pdf ) : m_pdf( pdf.clone() ) {}
    ~EvtPdfGen() { delete m_pdf; }

    result_type operator()() { return m_pdf->randomPoint(); }

  private:
    EvtPdf<T>* m_pdf;
};

template <class T>
class EvtPdfPred {
  public:
    typedef T argument_type;
    typedef bool result_type;

    EvtPdfPred() {}
    EvtPdfPred( const EvtPdf<T>& thePdf ) : m_pdf( thePdf.clone() ) {}
    EvtPdfPred( const EvtPdfPred& other ) :
        COPY_PTR( m_pdf ), COPY_MEM( m_pdfMax )
    {
    }
    ~EvtPdfPred() { delete m_pdf; }

    result_type operator()( argument_type p )
    {
        assert( m_pdf );
        assert( m_pdfMax.valueKnown() );

        double random = EvtRandom::Flat( 0., m_pdfMax.value() );
        return ( random <= m_pdf->evaluate( p ) );
    }

    EvtPdfMax<T> getMax() const { return m_pdfMax; }
    void setMax( const EvtPdfMax<T>& max ) { m_pdfMax = max; }
    template <class InputIterator>
    void compute_max( InputIterator it, InputIterator end, double factor = 1. )
    {
        T p = *it++;
        m_pdfMax = EvtPdfMax<T>( p, m_pdf->evaluate( p ) * factor );

        while ( !( it == end ) ) {
            T pp = *it++;
            double val = m_pdf->evaluate( pp ) * factor;
            if ( val > m_pdfMax.value() )
                m_pdfMax = EvtPdfMax<T>( pp, val );
        }
    }

  private:
    EvtPdf<T>* m_pdf;
    EvtPdfMax<T> m_pdfMax;
};

template <class T>
class EvtPdfUnary {
  public:
    typedef double result_type;
    typedef T argument_type;

    EvtPdfUnary() {}
    EvtPdfUnary( const EvtPdf<T>& thePdf ) : m_pdf( thePdf.clone() ) {}
    EvtPdfUnary( const EvtPdfUnary& other ) : COPY_PTR( m_pdf ) {}
    ~EvtPdfUnary() { delete m_pdf; }

    result_type operator()( argument_type p )
    {
        assert( m_pdf );
        double ret = m_pdf->evaluate( p );
        return ret;
    }

  private:
    EvtPdf<T>* m_pdf;
};

template <class T>
class EvtPdfDiv : public EvtPdf<T> {
  public:
    EvtPdfDiv() : m_num( 0 ), m_den( 0 ) {}
    EvtPdfDiv( const EvtPdf<T>& theNum, const EvtPdf<T>& theDen ) :
        EvtPdf<T>(), m_num( theNum.clone() ), m_den( theDen.clone() )
    {
    }
    EvtPdfDiv( const EvtPdfDiv<T>& other ) :
        EvtPdf<T>( other ), COPY_PTR( m_num ), COPY_PTR( m_den )
    {
    }
    virtual ~EvtPdfDiv()
    {
        delete m_num;
        delete m_den;
    }
    EvtPdf<T>* clone() const override { return new EvtPdfDiv( *this ); }

    double pdf( const T& p ) const override
    {
        double num = m_num->evaluate( p );
        double den = m_den->evaluate( p );
        assert( den != 0 );
        return num / den;
    }

  private:
    EvtPdf<T>* m_num;    // numerator
    EvtPdf<T>* m_den;    // denominator
};

template <class T>
EvtPdfMax<T> EvtPdf<T>::findMax( const EvtPdf<T>& pc, int N )
{
    EvtPdfPred<T> pred( *this );
    EvtPdfGen<T> gen( pc );
    pred.compute_max( iter( gen, N ), iter( gen ) );
    EvtPdfMax<T> p = pred.getMax();
    return p;
}

template <class T>
EvtValError EvtPdf<T>::findGenEff( const EvtPdf<T>& pc, int N, int nFindMax )
{
    assert( N > 0 || nFindMax > 0 );
    EvtPredGen<EvtPdfGen<T>, EvtPdfPred<T>> gen = accRejGen( pc, nFindMax );
    int i;
    for ( i = 0; i < N; i++ )
        gen();
    double eff = double( gen.getPassed() ) / double( gen.getTried() );
    double err = sqrt( double( gen.getPassed() ) ) / double( gen.getTried() );
    return EvtValError( eff, err );
}

template <class T>
EvtValError EvtPdf<T>::compute_mc_integral( const EvtPdf<T>& pc, int N )
{
    assert( N > 0 );

    EvtPdfDiv<T> pdfdiv( *this, pc );
    EvtPdfUnary<T> unary( pdfdiv );

    EvtPdfGen<T> gen( pc );
    EvtStreamInputIterator<T> begin = iter( gen, N );
    EvtStreamInputIterator<T> end;

    double sum = 0.;
    double sum2 = 0.;
    while ( !( begin == end ) ) {
        double value = pdfdiv.evaluate( *begin++ );
        sum += value;
        sum2 += value * value;
    }

    EvtValError x;
    if ( N > 0 ) {
        double av = sum / ( (double)N );
        if ( N > 1 ) {
            double dev2 = ( sum2 - av * av * N ) / ( (double)( N - 1 ) );
            // Due to numerical precision dev2 may sometimes be negative
            if ( dev2 < 0. )
                dev2 = 0.;
            double error = sqrt( dev2 / ( (double)N ) );
            x = EvtValError( av, error );
        } else
            x = EvtValError( av );
    }
    m_itg = x * pc.getItg();
    return m_itg;
}

template <class T>
T EvtPdf<T>::randomPoint()
{
    printf( "Function defined for analytic PDFs only\n" );
    assert( 0 );
    T temp;
    return temp;
}

template <class T>
EvtPredGen<EvtPdfGen<T>, EvtPdfPred<T>> EvtPdf<T>::accRejGen( const EvtPdf<T>& pc,
                                                              int nMax,
                                                              double factor )
{
    EvtPdfGen<T> gen( pc );
    EvtPdfDiv<T> pdfdiv( *this, pc );
    EvtPdfPred<T> pred( pdfdiv );
    pred.compute_max( iter( gen, nMax ), iter( gen ), factor );
    return EvtPredGen<EvtPdfGen<T>, EvtPdfPred<T>>( gen, pred );
}

#endif
