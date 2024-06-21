
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

#include "EvtGenBase/EvtComplex.hh"

#include <iostream>
#include <math.h>
using std::ostream;

ostream& operator<<( ostream& s, const EvtComplex& c )
{
    s << "(" << c.m_rpart << "," << c.m_ipart << ")";
    return s;
}

EvtComplex& EvtComplex::operator*=( EvtComplex c )
{
    double r = m_rpart * c.m_rpart - m_ipart * c.m_ipart;
    double i = m_rpart * c.m_ipart + m_ipart * c.m_rpart;

    m_rpart = r;
    m_ipart = i;

    return *this;
}

EvtComplex& EvtComplex::operator/=( EvtComplex c )
{
    double inv = 1.0 / ( c.m_rpart * c.m_rpart + c.m_ipart * c.m_ipart );

    double r = inv * ( m_rpart * c.m_rpart + m_ipart * c.m_ipart );
    double i = inv * ( m_ipart * c.m_rpart - m_rpart * c.m_ipart );

    m_rpart = r;
    m_ipart = i;

    return *this;
}
