
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

#include "EvtGenBase/EvtPartProp.hh"

#include "EvtGenBase/EvtAbsLineShape.hh"
#include "EvtGenBase/EvtFlatLineShape.hh"
#include "EvtGenBase/EvtManyDeltaFuncLineShape.hh"
#include "EvtGenBase/EvtRelBreitWignerBarrierFact.hh"

#include <cstdlib>
#include <string>

EvtPartProp::EvtPartProp() :
    m_id( -1, -1 ), m_idchgconj( -1, -1 ), m_chg3( 0 ), m_stdhep( 0 ), m_lundkc( 0 )
{
    m_ctau = 0.0;
    m_name = "*******";
    m_spintype = EvtSpinType::SCALAR;
}

EvtPartProp::EvtPartProp( const EvtPartProp& x )
{
    m_lineShape.reset( x.m_lineShape ? x.m_lineShape->clone() : nullptr );
    m_ctau = x.m_ctau;
    m_name = x.m_name;
    m_spintype = x.m_spintype;
    m_id = x.m_id;
    m_idchgconj = x.m_idchgconj;
    m_chg3 = x.m_chg3;
    m_stdhep = x.m_stdhep;
    m_lundkc = x.m_lundkc;
}

void EvtPartProp::setName( std::string pname )
{
    m_name = pname;
}

EvtPartProp& EvtPartProp::operator=( const EvtPartProp& x )
{
    m_lineShape.reset( x.m_lineShape ? x.m_lineShape->clone() : nullptr );

    m_ctau = x.m_ctau;
    m_name = x.m_name;
    m_chg3 = x.m_chg3;
    m_spintype = x.m_spintype;
    return *this;
}

void EvtPartProp::initLineShape( double mass, double width, double maxRange )
{
    m_lineShape = std::make_unique<EvtRelBreitWignerBarrierFact>( mass, width,
                                                                  maxRange,
                                                                  m_spintype );
}

void EvtPartProp::newLineShape( std::string type )
{
    double m = m_lineShape->getMass();
    double w = m_lineShape->getWidth();
    double mR = m_lineShape->getMaxRange();
    EvtSpinType::spintype st = m_lineShape->getSpinType();
    if ( type == "RELBW" ) {
        m_lineShape = std::make_unique<EvtRelBreitWignerBarrierFact>( m, w, mR,
                                                                      st );
    } else if ( type == "NONRELBW" ) {
        m_lineShape = std::make_unique<EvtAbsLineShape>( m, w, mR, st );
    } else if ( type == "FLAT" ) {
        m_lineShape = std::make_unique<EvtFlatLineShape>( m, w, mR, st );
    } else if ( type == "MANYDELTAFUNC" ) {
        m_lineShape = std::make_unique<EvtManyDeltaFuncLineShape>( m, w, mR, st );
    } else {
        m_lineShape.reset();
    }
}

void EvtPartProp::reSetMass( double mass )
{
    if ( !m_lineShape )
        ::abort();
    m_lineShape->reSetMass( mass );
}
void EvtPartProp::reSetWidth( double width )
{
    if ( !m_lineShape )
        ::abort();
    m_lineShape->reSetWidth( width );
}

void EvtPartProp::setPWForDecay( int spin, EvtId d1, EvtId d2 )
{
    if ( !m_lineShape )
        ::abort();
    m_lineShape->setPWForDecay( spin, d1, d2 );
}

void EvtPartProp::setPWForBirthL( int spin, EvtId par, EvtId othD )
{
    if ( !m_lineShape )
        ::abort();
    m_lineShape->setPWForBirthL( spin, par, othD );
}

void EvtPartProp::reSetMassMin( double mass )
{
    if ( !m_lineShape )
        ::abort();
    m_lineShape->reSetMassMin( mass );
}
void EvtPartProp::reSetMassMax( double mass )
{
    if ( !m_lineShape )
        ::abort();
    m_lineShape->reSetMassMax( mass );
}
void EvtPartProp::reSetBlatt( double blatt )
{
    if ( !m_lineShape )
        ::abort();
    m_lineShape->reSetBlatt( blatt );
}
void EvtPartProp::reSetBlattBirth( double blatt )
{
    if ( !m_lineShape )
        ::abort();
    m_lineShape->reSetBlattBirth( blatt );
}
void EvtPartProp::includeBirthFactor( bool yesno )
{
    if ( !m_lineShape )
        ::abort();
    m_lineShape->includeBirthFactor( yesno );
}
void EvtPartProp::includeDecayFactor( bool yesno )
{
    if ( !m_lineShape )
        ::abort();
    m_lineShape->includeDecayFactor( yesno );
}
