
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

#ifndef EVT_MULTI_CHANNEL_PARSER
#define EVT_MULTI_CHANNEL_PARSER

#include "EvtGenBase/EvtComplex.hh"

#include <string>
#include <vector>

// Parse decay file tokens from a file or a vector

class EvtDecayMode;

enum
{
    POLAR_RAD,
    POLAR_DEG,
    CARTESIAN
};

class EvtMultiChannelParser {
  public:
    EvtMultiChannelParser() : m_pdfMax( -1. ), m_nScan( 0 ), m_dm( 0. ) {}
    ~EvtMultiChannelParser() {}

    static EvtDecayMode getDecayMode( const char* file );

    void parse( const char* file, const char* model );
    void parse( const std::vector<std::string>& v );

    static void parseComplexCoef( size_t& i, const std::vector<std::string>& v,
                                  EvtComplex& c, int& format );
    static double parseRealCoef( int& i, const std::vector<std::string>& v );
    static bool isKeyword( const std::string& s );

    inline double pdfMax() const { return m_pdfMax; }
    inline int nScan() const { return m_nScan; }
    inline double dm() const { return m_dm; }
    inline double mixPhase() const { return m_mixPhase; }
    inline double mixAmpli() const { return m_mixAmpli; }

    inline std::vector<std::string> amp( int i ) const { return m_amp[i]; }
    inline std::vector<std::string> ampConj( int i ) const
    {
        return m_ampConj[i];
    }
    inline EvtComplex ampCoef( int i ) const { return m_ampCoef[i]; }
    inline EvtComplex ampConjCoef( int i ) const { return m_ampConjCoef[i]; }

    inline int coefFormat( int i ) const { return m_coefFormat[i]; }
    inline int coefConjFormat( int i ) const { return m_coefConjFormat[i]; }

    inline int getNAmp() const { return m_amp.size(); }
    inline int getNAmpConj() const { return m_ampConj.size(); }

  private:
    double m_pdfMax;
    int m_nScan;
    double m_dm;
    double m_mixPhase;
    double m_mixAmpli;

    std::vector<std::vector<std::string>> m_amp;
    std::vector<std::vector<std::string>> m_ampConj;
    std::vector<EvtComplex> m_ampCoef;
    std::vector<int> m_coefFormat;
    std::vector<EvtComplex> m_ampConjCoef;
    std::vector<int> m_coefConjFormat;
};

#endif
