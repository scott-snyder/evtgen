
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

#ifndef EVT_PDF_SUM_HH
#define EVT_PDF_SUM_HH

#include "EvtGenBase/EvtPdf.hh"

#include <cassert>
#include <vector>

// Sum of PDF functions.

template <class T>
class EvtPdfSum : public EvtPdf<T> {
  public:
    EvtPdfSum() {}
    EvtPdfSum( const EvtPdfSum<T>& other );
    virtual ~EvtPdfSum();
    EvtPdfSum* clone() const override { return new EvtPdfSum( *this ); }

    // Manipulate terms and coefficients

    void addTerm( double c, const EvtPdf<T>& pdf )
    {
        assert( c >= 0. );
        m_c.push_back( c );
        m_term.push_back( pdf.clone() );
    }

    void addOwnedTerm( double c, std::unique_ptr<EvtPdf<T>> pdf )
    {
        m_c.push_back( c );
        m_term.push_back( pdf.release() );
    }

    size_t nTerms() const { return m_term.size(); }    // number of terms

    inline double c( int i ) const { return m_c[i]; }
    inline EvtPdf<T>* getPdf( int i ) const { return m_term[i]; }

    // Integrals

    EvtValError compute_integral() const override;
    EvtValError compute_integral( int N ) const override;
    T randomPoint() override;

  protected:
    double pdf( const T& p ) const override;

    std::vector<double> m_c;           // coefficients
    std::vector<EvtPdf<T>*> m_term;    // pointers to pdfs
};

template <class T>
EvtPdfSum<T>::EvtPdfSum( const EvtPdfSum<T>& other ) : EvtPdf<T>( other )
{
    for ( size_t i = 0; i < other.nTerms(); i++ ) {
        m_c.push_back( other.m_c[i] );
        m_term.push_back( other.m_term[i]->clone() );
    }
}

template <class T>
EvtPdfSum<T>::~EvtPdfSum()
{
    for ( size_t i = 0; i < m_c.size(); i++ ) {
        delete m_term[i];
    }
}

template <class T>
double EvtPdfSum<T>::pdf( const T& p ) const
{
    double ret = 0.;
    for ( size_t i = 0; i < m_c.size(); i++ ) {
        ret += m_c[i] * m_term[i]->evaluate( p );
    }
    return ret;
}

/*
 * Compute the sum integral by summing all term integrals.
 */

template <class T>
EvtValError EvtPdfSum<T>::compute_integral() const
{
    EvtValError itg( 0.0, 0.0 );
    for ( size_t i = 0; i < nTerms(); i++ ) {
        itg += m_c[i] * m_term[i]->getItg();
    }
    return itg;
}

template <class T>
EvtValError EvtPdfSum<T>::compute_integral( int N ) const
{
    EvtValError itg( 0.0, 0.0 );
    for ( size_t i = 0; i < nTerms(); i++ )
        itg += m_c[i] * m_term[i]->getItg( N );
    return itg;
}

/*
 * Sample points randomly according to the sum of PDFs. First throw a random number uniformly
 * between zero and the value of the sum integral. Using this random number select one
 * of the PDFs. The generate a random point according to that PDF.
 */

template <class T>
T EvtPdfSum<T>::randomPoint()
{
    if ( !this->m_itg.valueKnown() )
        this->m_itg = compute_integral();

    double max = this->m_itg.value();
    double rnd = EvtRandom::Flat( 0, max );

    double sum = 0.;
    size_t i;
    for ( i = 0; i < nTerms(); i++ ) {
        double itg = m_term[i]->getItg().value();
        sum += m_c[i] * itg;
        if ( sum > rnd )
            break;
    }

    return m_term[i]->randomPoint();
}

#endif
