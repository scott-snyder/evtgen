
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

#ifndef EVT_PRED_GEN_HH
#define EVT_PRED_GEN_HH

#include <cstdio>

// A predicate is applied to a generator to get another generator.
// Accept-reject can be implemented in this way.
//
//           Predicate
// Generator    ->     Generator

template <class Generator, class Predicate>
class EvtPredGen {
  public:
    typedef typename Generator::result_type result_type;

    EvtPredGen() : m_tried( 0 ), m_passed( 0 ) {}

    EvtPredGen( Generator gen, Predicate pred ) :
        m_gen( gen ), m_pred( pred ), m_tried( 0 ), m_passed( 0 )
    {
    }

    EvtPredGen( const EvtPredGen& other ) :
        m_gen( other.m_gen ),
        m_pred( other.m_pred ),
        m_tried( other.m_tried ),
        m_passed( other.m_passed )
    {
    }

    ~EvtPredGen() {}

    result_type operator()()
    {
        int i = 0;
        int MAX = 10000;
        while ( i++ < MAX ) {
            m_tried++;
            result_type point = m_gen();
            if ( m_pred( point ) ) {
                m_passed++;
                return point;
            }
        }

        printf( "No random point generated after %d attempts\n", MAX );
        printf( "Sharp peak? Consider using pole compensation.\n" );
        printf( "I will now pick a point at random to return.\n" );
        return m_gen();
    }

    inline int getTried() const { return m_tried; }
    inline int getPassed() const { return m_passed; }

  protected:
    Generator m_gen;
    Predicate m_pred;
    int m_tried;
    int m_passed;
};

#endif
