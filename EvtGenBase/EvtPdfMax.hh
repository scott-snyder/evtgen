
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

#ifndef EVT_PDF_MAX_HH
#define EVT_PDF_MAX_HH

#include "EvtGenBase/EvtMacros.hh"

#include <cassert>

// Pdf maximum and its location

template <class Point>
class EvtPdfMax {
  public:
    EvtPdfMax() : m_value( -1 ), m_valueKnown( false ), m_locKnown( false ) {}
    EvtPdfMax( double value ) :
        m_value( value ), m_valueKnown( true ), m_locKnown( false )
    {
    }
    EvtPdfMax( Point p, double value ) :
        m_value( value ), m_valueKnown( true ), m_locKnown( true ), m_loc( p )
    {
    }

    bool valueKnown() const { return m_valueKnown; }
    double value() const
    {
        assert( m_valueKnown );
        return m_value;
    }
    bool locKnown() const { return m_locKnown; }
    Point loc() const
    {
        assert( m_locKnown );
        return m_loc;
    }

  private:
    double m_value;
    bool m_valueKnown;
    bool m_locKnown;
    Point m_loc;
};

#endif
