
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

#ifndef EVTPARTPROP_HH
#define EVTPARTPROP_HH

#include "EvtGenBase/EvtAbsLineShape.hh"
#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtSpinType.hh"

#include <memory>
#include <string>

class EvtPartProp {
  public:
    EvtPartProp();
    EvtPartProp( const EvtPartProp& x );

    double getMass() { return m_lineShape->getMass(); }
    double getMassMin() { return m_lineShape->getMassMin(); }
    double getMassMax() { return m_lineShape->getMassMax(); }
    double getMaxRange() { return m_lineShape->getMaxRange(); }
    double getWidth() { return m_lineShape->getWidth(); }

    double getRandMass( EvtId* parId, int nDaug, EvtId* dauId, EvtId* othDauId,
                        double maxMass, double* dauMasses )
    {
        return m_lineShape->getRandMass( parId, nDaug, dauId, othDauId, maxMass,
                                         dauMasses );
    }
    double getMassProb( double mass, double massPar, int nDaug, double* massDau )
    {
        return m_lineShape->getMassProb( mass, massPar, nDaug, massDau );
    }

    double getctau() { return m_ctau; }
    void setctau( double tau ) { m_ctau = tau; }

    int getChg3() { return m_chg3; }
    void setChg3( int c3 ) { m_chg3 = c3; }

    EvtSpinType::spintype getSpinType() { return m_spintype; }
    void setSpinType( EvtSpinType::spintype stype ) { m_spintype = stype; }

    const std::string& getName() { return m_name; }
    void setName( std::string pname );

    EvtId getId() const { return m_id; }
    void setId( EvtId id ) { m_id = id; }

    EvtId getIdChgConj() { return m_idchgconj; }
    void setIdChgConj( EvtId idchgconj ) { m_idchgconj = idchgconj; }

    int getStdHep() const { return m_stdhep; }
    void setStdHep( int stdhep ) { m_stdhep = stdhep; }

    int getLundKC() { return m_lundkc; }
    void setLundKC( int lundkc ) { m_lundkc = lundkc; }

    EvtAbsLineShape* getLineShape() { return m_lineShape.get(); }
    void initLineShape( double mass, double width, double maxRange );
    //  void initLineShape(double mass, double width, double maxRange, double mDaug1, double mDaug2, int l);

    // setLineShape takes ownership of l
    void setLineShape( EvtAbsLineShape* l ) { m_lineShape.reset( l ); }
    double rollMass() { return m_lineShape->rollMass(); }

    EvtPartProp& operator=( const EvtPartProp& x );

    void reSetMass( double mass );
    void reSetWidth( double width );

    void reSetMassMin( double mass );
    void reSetMassMax( double mass );
    void reSetBlatt( double blatt );
    void reSetBlattBirth( double blatt );
    void includeBirthFactor( bool yesno );
    void includeDecayFactor( bool yesno );
    void newLineShape( std::string type );
    void setPWForDecay( int spin, EvtId d1, EvtId d2 );
    void setPWForBirthL( int spin, EvtId par, EvtId othD );

  private:
    std::unique_ptr<EvtAbsLineShape> m_lineShape;

    double m_ctau;
    EvtId m_id;
    EvtId m_idchgconj;
    EvtSpinType::spintype m_spintype;
    int m_chg3;
    int m_stdhep;
    int m_lundkc;
    std::string m_name;
};

#endif
