
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

#ifndef EVT_MASSAMP_HH
#define EVT_MASSAMP_HH

#include "EvtGenBase/EvtAmplitude.hh"
#include "EvtGenBase/EvtPoint1D.hh"
#include "EvtGenBase/EvtPropBreitWignerRel.hh"
#include "EvtGenBase/EvtTwoBodyVertex.hh"

#include <memory>
// Relativistic lineshape for a two-body decay of a resonance to two
// pseudoscalars. The mass dependence of the width and the vertex factors
// are included in the calculation.

class EvtMassAmp : public EvtAmplitude<EvtPoint1D> {
  public:
    EvtMassAmp( const EvtPropBreitWignerRel& prop, const EvtTwoBodyVertex& vd );
    EvtMassAmp( const EvtMassAmp& other );
    EvtMassAmp& operator=( const EvtMassAmp& other );

    EvtComplex amplitude( const EvtPoint1D& p ) const override;

    EvtAmplitude<EvtPoint1D>* clone() const override
    {
        return new EvtMassAmp( *this );
    }

    void setBirthVtx( const EvtTwoBodyVertex& vb )
    {
        m_vb = std::make_unique<EvtTwoBodyVertex>( vb );
    }

    void addBirthFact() { m_useBirthFact = true; }
    void addDeathFact() { m_useDeathFact = true; }
    void addBirthFactFF() { m_useBirthFactFF = true; }
    void addDeathFactFF() { m_useDeathFactFF = true; }

  private:
    EvtPropBreitWignerRel m_prop;
    EvtTwoBodyVertex m_vd;
    std::unique_ptr<EvtTwoBodyVertex> m_vb;

    bool m_useBirthFact;
    bool m_useDeathFact;
    bool m_useBirthFactFF;
    bool m_useDeathFactFF;
};

#endif
