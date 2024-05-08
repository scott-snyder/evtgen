
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

#ifndef EVTABSLINESHAPE_HH
#define EVTABSLINESHAPE_HH

#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtSpinType.hh"

#include <vector>

class EvtAbsLineShape {
  public:
    EvtAbsLineShape() = default;
    EvtAbsLineShape( double mass, double width, double maxRange,
                     EvtSpinType::spintype sp );
    virtual ~EvtAbsLineShape() = default;
    EvtAbsLineShape& operator=( const EvtAbsLineShape& x );
    EvtAbsLineShape( const EvtAbsLineShape& x );

    double getMass() { return m_mass; }
    double getMassMin() { return m_massMin; }
    double getMassMax() { return m_massMax; }
    double getMaxRange() { return m_maxRange; }
    double getWidth() { return m_width; }
    EvtSpinType::spintype getSpinType() { return m_spin; }
    virtual double rollMass();
    virtual EvtAbsLineShape* clone();

    void reSetMass( double mass ) { m_mass = mass; }
    void reSetWidth( double width ) { m_width = width; }
    void reSetMassMin( double mass ) { m_massMin = mass; }
    void reSetMassMax( double mass ) { m_massMax = mass; }
    virtual void reSetBlatt( double /*blatt*/ ){};
    virtual void reSetBlattBirth( double /*blatt*/ ){};
    void includeBirthFactor( bool yesno ) { m_includeBirthFact = yesno; }
    void includeDecayFactor( bool yesno ) { m_includeDecayFact = yesno; }
    void setPWForDecay( int spin, EvtId d1, EvtId d2 )
    {
        m_userSetPW.push_back( spin );
        m_userSetPWD1.push_back( d1 );
        m_userSetPWD2.push_back( d2 );
    }
    void setPWForBirthL( int spin, EvtId par, EvtId othD )
    {
        m_userSetBirthPW.push_back( spin );
        m_userSetBirthOthD.push_back( othD );
        m_userSetBirthPar.push_back( par );
    }

    virtual double getRandMass( EvtId* parId, int nDaug, EvtId* dauId,
                                EvtId* othDaugId, double maxMass,
                                double* dauMasses );
    virtual double getMassProb( double mass, double massPar, int nDaug,
                                double* massDau );

  protected:
    bool m_includeDecayFact;
    bool m_includeBirthFact;
    double m_mass;
    double m_massMin;
    double m_massMax;
    double m_width;
    double m_maxRange;

    // allow for special cases where the default method of picking the
    //lowest allowed partial wave for a decay is not the right answer.
    // string is "<spin> <daughter1> <daughter2>"
    //new 9/12/2003 Lange
    std::vector<EvtId> m_userSetPWD1, m_userSetPWD2;
    std::vector<int> m_userSetPW;

    // also do it for birth factors
    std::vector<EvtId> m_userSetBirthPar, m_userSetBirthOthD;
    std::vector<int> m_userSetBirthPW;

    EvtSpinType::spintype m_spin;
};

#endif
