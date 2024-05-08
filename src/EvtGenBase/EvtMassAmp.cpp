
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

#include "EvtGenBase/EvtMassAmp.hh"

#include "EvtGenBase/EvtPatches.hh"

EvtMassAmp::EvtMassAmp( const EvtPropBreitWignerRel& prop,
                        const EvtTwoBodyVertex& vd ) :
    EvtAmplitude<EvtPoint1D>(),
    m_prop( prop ),
    m_vd( vd ),
    m_useBirthFact( false ),
    m_useDeathFact( false ),
    m_useBirthFactFF( false ),
    m_useDeathFactFF( false )
{
}

EvtMassAmp::EvtMassAmp( const EvtMassAmp& other ) :
    EvtAmplitude<EvtPoint1D>( other ),
    m_prop( other.m_prop ),
    m_vd( other.m_vd ),
    m_vb( other.m_vb ? new EvtTwoBodyVertex( *other.m_vb ) : nullptr ),
    m_useBirthFact( other.m_useBirthFact ),
    m_useDeathFact( other.m_useDeathFact ),
    m_useBirthFactFF( other.m_useBirthFactFF ),
    m_useDeathFactFF( other.m_useDeathFactFF )
{
}

EvtMassAmp& EvtMassAmp::operator=( const EvtMassAmp& other )
{
    EvtAmplitude<EvtPoint1D>::operator=( other );
    m_prop = other.m_prop;
    m_vd = other.m_vd;
    m_vb.reset( other.m_vb ? new EvtTwoBodyVertex( *other.m_vb ) : nullptr );
    m_useBirthFact = other.m_useBirthFact;
    m_useDeathFact = other.m_useDeathFact;
    m_useBirthFactFF = other.m_useBirthFactFF;
    m_useDeathFactFF = other.m_useDeathFactFF;
    return *this;
}

EvtComplex EvtMassAmp::amplitude( const EvtPoint1D& p ) const
{
    // Modified vertex

    double m = p.value();
    // keep things from crashing..

    if ( m < ( m_vd.mA() + m_vd.mB() ) )
        return EvtComplex( 0., 0. );

    EvtTwoBodyKine vd( m_vd.mA(), m_vd.mB(), m );

    // Compute mass-dependent width for relativistic propagator

    EvtPropBreitWignerRel bw( m_prop.m0(), m_prop.g0() * m_vd.widthFactor( vd ) );
    EvtComplex amp = bw.evaluate( m );

    // Birth vertex factors

    if ( m_useBirthFact ) {
        assert( m_vb );
        if ( ( m + m_vb->mB() ) < m_vb->mAB() ) {
            EvtTwoBodyKine vb( m, m_vb->mB(), m_vb->mAB() );
            amp *= m_vb->phaseSpaceFactor( vb, EvtTwoBodyKine::AB );
            amp *= sqrt( ( vb.p() / m_vb->pD() ) );

            if ( m_useBirthFactFF ) {
                assert( m_vb );
                amp *= m_vb->formFactor( vb );
            }
        } else {
            if ( m_vb->L() != 0 )
                amp = 0.;
        }
    }

    // Decay vertex factors

    if ( m_useDeathFact ) {
        amp *= m_vd.phaseSpaceFactor( vd, EvtTwoBodyKine::AB );
        amp *= sqrt( ( vd.p() / m_vd.pD() ) );
    }
    if ( m_useDeathFactFF )
        amp *= m_vd.formFactor( vd );

    return amp;
}
