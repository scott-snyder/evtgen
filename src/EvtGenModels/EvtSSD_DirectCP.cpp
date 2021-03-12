
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

#include "EvtGenModels/EvtSSD_DirectCP.hh"

#include "EvtGenBase/EvtCPUtil.hh"
#include "EvtGenBase/EvtConst.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtScalarParticle.hh"
#include "EvtGenBase/EvtTensor4C.hh"
#include "EvtGenBase/EvtVector4C.hh"

#include <cstdlib>
#include <string>

std::string EvtSSD_DirectCP::getName() const
{
    return "SSD_DirectCP";
}

EvtDecayBase* EvtSSD_DirectCP::clone() const
{
    return new EvtSSD_DirectCP;
}

void EvtSSD_DirectCP::init()
{
    // check that there is 1 argument and 2-body decay

    checkNArg( 1 );
    checkNDaug( 2 );

    const EvtSpinType::spintype d1type = EvtPDL::getSpinType( getDaug( 0 ) );
    const EvtSpinType::spintype d2type = EvtPDL::getSpinType( getDaug( 1 ) );

    if ( ( !( d1type == EvtSpinType::SCALAR || d2type == EvtSpinType::SCALAR ) ) ||
         ( !( ( d2type == EvtSpinType::SCALAR ) ||
              ( d2type == EvtSpinType::VECTOR ) ||
              ( d2type == EvtSpinType::TENSOR ) ) ) ||
         ( !( ( d1type == EvtSpinType::SCALAR ) ||
              ( d1type == EvtSpinType::VECTOR ) ||
              ( d1type == EvtSpinType::TENSOR ) ) ) ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "EvtSSD_DirectCP generator expected "
            << "one of the daugters to be a scalar, "
            << "the other either scalar, vector, or tensor, "
            << "found:" << EvtPDL::name( getDaug( 0 ) ).c_str() << " and "
            << EvtPDL::name( getDaug( 1 ) ).c_str() << std::endl;
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Will terminate execution!" << std::endl;
        ::abort();
    }

    m_acp = getArg( 0 );    // A_CP defined as A_CP = (BR(fbar)-BR(f))/(BR(fbar)+BR(f))
}

void EvtSSD_DirectCP::initProbMax()
{
    double theProbMax = 1.;

    const EvtSpinType::spintype d2type = EvtPDL::getSpinType( getDaug( 1 ) );
    const EvtSpinType::spintype d1type = EvtPDL::getSpinType( getDaug( 0 ) );

    if ( d1type != EvtSpinType::SCALAR || d2type != EvtSpinType::SCALAR ) {
        // Create a scalar parent at rest and initialize it
        // Use noLifeTime() cludge to avoid generating random numbers

        EvtScalarParticle parent{};
        parent.noLifeTime();
        parent.init( getParentId(),
                     EvtVector4R( EvtPDL::getMass( getParentId() ), 0, 0, 0 ) );
        parent.setDiagonalSpinDensity();

        // Create daughters and initialize amplitude
        EvtAmp amp;
        EvtId daughters[2] = { getDaug( 0 ), getDaug( 1 ) };
        amp.init( getParentId(), 2, daughters );
        parent.makeDaughters( 2, daughters );

        const int scalarDaughterIndex = d1type == EvtSpinType::SCALAR ? 0 : 1;
        const int nonScalarDaughterIndex = d1type == EvtSpinType::SCALAR ? 1 : 0;

        EvtParticle* scalarDaughter = parent.getDaug( scalarDaughterIndex );
        EvtParticle* nonScalarDaughter = parent.getDaug( nonScalarDaughterIndex );

        scalarDaughter->noLifeTime();
        nonScalarDaughter->noLifeTime();

        EvtSpinDensity rho;
        rho.setDiag( parent.getSpinStates() );

        // Momentum of daughters in parent's frame
        const double parentMass = EvtPDL::getMass( getParentId() );
        const double sdMass = EvtPDL::getMass( getDaug( scalarDaughterIndex ) );
        const double ndMass = EvtPDL::getMass( getDaug( nonScalarDaughterIndex ) );
        const double pstar =
            sqrt( pow( parentMass, 2 ) - pow( ( sdMass + ndMass ), 2 ) ) *
            sqrt( pow( parentMass, 2 ) - pow( ( ndMass - sdMass ), 2 ) ) /
            ( 2 * parentMass );

        EvtVector4R p4_sd, p4_nd;

        const int nsteps = 16;

        double prob_max = 0;
        double theta_max = 0;

        for ( int i = 0; i <= nsteps; i++ ) {
            const double theta = i * EvtConst::pi / nsteps;

            p4_sd.set( sqrt( pow( pstar, 2 ) + pow( sdMass, 2 ) ), 0,
                       +pstar * sin( theta ), +pstar * cos( theta ) );

            p4_nd.set( sqrt( pow( pstar, 2 ) + pow( ndMass, 2 ) ), 0,
                       -pstar * sin( theta ), -pstar * cos( theta ) );

            scalarDaughter->init( getDaug( scalarDaughterIndex ), p4_sd );
            nonScalarDaughter->init( getDaug( nonScalarDaughterIndex ), p4_nd );

            calcAmp( parent, amp );

            const double i_prob = rho.normalizedProb( amp.getSpinDensity() );

            if ( i_prob > prob_max ) {
                prob_max = i_prob;
                theta_max = theta;
            }
        }

        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << "  - probability " << prob_max
            << " found at theta = " << theta_max << std::endl;
        theProbMax *= 1.01 * prob_max;
    }

    setProbMax( theProbMax );

    EvtGenReport( EVTGEN_INFO, "EvtGen" )
        << " EvtSSD_DirectCP - set up maximum probability to " << theProbMax
        << std::endl;
}

void EvtSSD_DirectCP::decay( EvtParticle* parent )
{
    bool flip = false;
    EvtId daugs[2];

    // decide it is B or Bbar:
    if ( EvtRandom::Flat( 0., 1. ) < ( ( 1. - m_acp ) / 2. ) ) {
        // it is a B
        if ( EvtPDL::getStdHep( getParentId() ) < 0 )
            flip = true;
    } else {
        // it is a Bbar
        if ( EvtPDL::getStdHep( getParentId() ) > 0 )
            flip = true;
    }

    if ( flip ) {
        if ( ( isB0Mixed( *parent ) ) || ( isBsMixed( *parent ) ) ) {
            parent->getParent()->setId(
                EvtPDL::chargeConj( parent->getParent()->getId() ) );
            parent->setId( EvtPDL::chargeConj( parent->getId() ) );
        } else {
            parent->setId( EvtPDL::chargeConj( parent->getId() ) );
        }
        daugs[0] = EvtPDL::chargeConj( getDaug( 0 ) );
        daugs[1] = EvtPDL::chargeConj( getDaug( 1 ) );
    } else {
        daugs[0] = getDaug( 0 );
        daugs[1] = getDaug( 1 );
    }

    parent->initializePhaseSpace( 2, daugs );

    calcAmp( *parent, m_amp2 );
}

bool EvtSSD_DirectCP::isB0Mixed( const EvtParticle& p )
{
    if ( !( p.getParent() ) )
        return false;

    static const EvtId B0 = EvtPDL::getId( "B0" );
    static const EvtId B0B = EvtPDL::getId( "anti-B0" );

    if ( ( p.getId() != B0 ) && ( p.getId() != B0B ) )
        return false;

    if ( ( p.getParent()->getId() == B0 ) || ( p.getParent()->getId() == B0B ) )
        return true;

    return false;
}

bool EvtSSD_DirectCP::isBsMixed( const EvtParticle& p )
{
    if ( !( p.getParent() ) )
        return false;

    static const EvtId BS0 = EvtPDL::getId( "B_s0" );
    static const EvtId BSB = EvtPDL::getId( "anti-B_s0" );

    if ( ( p.getId() != BS0 ) && ( p.getId() != BSB ) )
        return false;

    if ( ( p.getParent()->getId() == BS0 ) || ( p.getParent()->getId() == BSB ) )
        return true;

    return false;
}

std::string EvtSSD_DirectCP::getParamName( int i )
{
    switch ( i ) {
        case 0:
            return "ACP";
        default:
            return "";
    }
}

void EvtSSD_DirectCP::calcAmp( const EvtParticle& parent, EvtAmp& amp ) const
{
    const EvtSpinType::spintype d1type = EvtPDL::getSpinType(
        parent.getDaug( 0 )->getId() );
    const EvtSpinType::spintype d2type = EvtPDL::getSpinType(
        parent.getDaug( 1 )->getId() );

    if ( d1type == EvtSpinType::SCALAR && d2type == EvtSpinType::SCALAR ) {
        amp.vertex( 1. );
        return;
    }

    const EvtVector4R p4_parent = parent.getP4Restframe();

    const int nonScalarDaughterIndex = d1type == EvtSpinType::SCALAR ? 1 : 0;

    const EvtParticle& daughter = *parent.getDaug( nonScalarDaughterIndex );

    const EvtSpinType::spintype nonScalarType = EvtPDL::getSpinType(
        daughter.getId() );

    if ( nonScalarType == EvtSpinType::VECTOR ) {
        const EvtVector4R momv = daughter.getP4();
        const double norm = momv.mass() / ( momv.d3mag() * parent.mass() );

        amp.vertex( 0, norm * p4_parent * ( daughter.epsParent( 0 ) ) );
        amp.vertex( 1, norm * p4_parent * ( daughter.epsParent( 1 ) ) );
        amp.vertex( 2, norm * p4_parent * ( daughter.epsParent( 2 ) ) );

    }

    // This is for the EvtSpinType::TENSOR case.
    else {
        const double norm = daughter.mass() * daughter.mass() /
                            ( p4_parent.mass() * daughter.getP4().d3mag() *
                              daughter.getP4().d3mag() );

        amp.vertex( 0, norm * daughter.epsTensorParent( 0 ).cont1( p4_parent ) *
                           p4_parent );
        amp.vertex( 1, norm * daughter.epsTensorParent( 1 ).cont1( p4_parent ) *
                           p4_parent );
        amp.vertex( 2, norm * daughter.epsTensorParent( 2 ).cont1( p4_parent ) *
                           p4_parent );
        amp.vertex( 3, norm * daughter.epsTensorParent( 3 ).cont1( p4_parent ) *
                           p4_parent );
        amp.vertex( 4, norm * daughter.epsTensorParent( 4 ).cont1( p4_parent ) *
                           p4_parent );
    }
}
