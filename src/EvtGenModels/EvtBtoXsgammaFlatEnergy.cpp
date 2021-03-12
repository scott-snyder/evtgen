
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

#include "EvtGenModels/EvtBtoXsgammaFlatEnergy.hh"

#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtReport.hh"

#include "EvtGenModels/EvtBtoXsgamma.hh"

#include <cmath>
#include <cstdlib>

using std::endl;

void EvtBtoXsgammaFlatEnergy::init( int nArg, double* args )
{
    if ( ( nArg ) > 3 || ( nArg > 1 && nArg < 3 ) ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "EvtBtoXsgamma generator model "
            << "EvtBtoXsgammaFlatEnergy expected "
            << "either 1(default config) or two arguments but found: " << nArg
            << endl;
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Will terminate execution!" << endl;
        ::abort();
    }
    m_mB0 = 5.2794;
    double mPi = 0.140;
    double mK = 0.494;
    if ( nArg == 1 ) {
        m_eMin = 1.7;
        //Invariant mass of Xsd must be greater the m_pi+m_K leads to
        //Egamma < (m_B**2-(m_pi+m_k)**2)/(2m_B)
        m_eMax = ( pow( m_mB0, 2 ) - pow( mPi + mK, 2 ) ) / ( 2.0 * m_mB0 );
    } else {
        m_eMin = args[1];
        m_eMax = args[2];
    }
    if ( m_eMax > ( pow( m_mB0, 2 ) - pow( mPi + mK, 2 ) ) / ( 2.0 * m_mB0 ) ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Emax greater than Kinematic limit" << endl;
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Reset to the kinematic limit" << endl;
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "(m_B**2-(m_pi+m_k)**2)/(2m_B)" << endl;
        m_eMax = ( pow( m_mB0, 2 ) - pow( mPi + mK, 2 ) ) / ( 2.0 * m_mB0 );
    }
    m_eRange = m_eMax - m_eMin;
}

double EvtBtoXsgammaFlatEnergy::GetMass( int /*Xscode*/ )
{
    double eGamma = EvtRandom::Flat( m_eRange ) + m_eMin;
    double mH = sqrt( pow( m_mB0, 2 ) - 2.0 * m_mB0 * eGamma );
    return mH;
}
