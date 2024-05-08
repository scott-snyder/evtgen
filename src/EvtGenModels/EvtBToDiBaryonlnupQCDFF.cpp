
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

#include "EvtGenModels/EvtBToDiBaryonlnupQCDFF.hh"

#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"

EvtBToDiBaryonlnupQCDFF::EvtBToDiBaryonlnupQCDFF() : m_DPars(), m_nDPars( 0 )
{
}

EvtBToDiBaryonlnupQCDFF::EvtBToDiBaryonlnupQCDFF( std::vector<double>& DParameters ) :
    m_DPars( DParameters ), m_nDPars( DParameters.size() )
{
}

void EvtBToDiBaryonlnupQCDFF::getFF( EvtParticle*, double dibaryonMass,
                                     EvtBToDiBaryonlnupQCDFF::FormFactors& FF ) const
{
    if ( m_nDPars == 6 && dibaryonMass > 0.0 ) {
        // 5/3*[1/M^2]^3
        double t = 5.0 / ( 3.0 * pow( dibaryonMass, 6.0 ) );

        double Dp = m_DPars[0];
        double Dpb = m_DPars[1];
        double D2 = m_DPars[2];
        double D3 = m_DPars[3];
        double D4 = m_DPars[4];
        double D5 = m_DPars[5];

        FF.m_F1 = ( Dp + 0.2 * Dpb ) * t;
        FF.m_F2 = -D2 * t;
        FF.m_F3 = -D3 * t;
        FF.m_F4 = -D4 * t;
        FF.m_F5 = -D5 * t;

        FF.m_G1 = ( Dp - 0.2 * Dpb ) * t;
        FF.m_G2 = -FF.m_F2;
        FF.m_G3 = -FF.m_F3;
        FF.m_G4 = -FF.m_F4;
        FF.m_G5 = -FF.m_F5;
    }
}

void EvtBToDiBaryonlnupQCDFF::getDiracFF(
    EvtParticle* parent, double dibaryonMass,
    EvtBToDiBaryonlnupQCDFF::FormFactors& FF ) const
{
    this->getFF( parent, dibaryonMass, FF );
}

void EvtBToDiBaryonlnupQCDFF::getRaritaFF(
    EvtParticle* parent, double dibaryonMass,
    EvtBToDiBaryonlnupQCDFF::FormFactors& FF ) const
{
    this->getFF( parent, dibaryonMass, FF );
}
