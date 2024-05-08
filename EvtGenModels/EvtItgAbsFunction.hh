
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

#ifndef EVTITGABSFUNCTION_HH
#define EVTITGABSFUNCTION_HH

//-------------
// C Headers --
//-------------
extern "C" {
}

// Description:
//      Abstraction of a generic function for use in integration methods elsewhere
//      in this package. (Stolen and modified from the BaBar IntegrationUtils package
//      - author: Phil Strother).

class EvtItgAbsFunction {
  public:
    // Constructors

    EvtItgAbsFunction( double lowerRange, double upperRange );

    // Destructor
    virtual ~EvtItgAbsFunction() = default;

    virtual double value( double x ) const;

    virtual double operator()( double x ) const;

    // Selectors (const)

    inline double upperRange() const { return m_upperRange; }
    inline double lowerRange() const { return m_lowerRange; }
    inline void getRange( double& lower, double& upper ) const
    {
        lower = m_lowerRange;
        upper = m_upperRange;
    }
    virtual void setCoeff( int, int, double ) = 0;
    virtual double getCoeff( int, int ) = 0;

  protected:
    virtual double myFunction( double x ) const = 0;
    void setRange( double x1, double x2 )
    {
        m_lowerRange = x1;
        m_upperRange = x2;
    };

  private:
    double m_upperRange;
    double m_lowerRange;
};

#endif    // EVTITGABSFUNCTION_HH
