
/***********************************************************************
* Copyright 1998-2022 CERN for the benefit of the EvtGen authors       *
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

#include "EvtGenModels/EvtRareLbToLllFF.hh"

#include "EvtGenBase/EvtIdSet.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtVector4R.hh"

//-----------------------------------------------------------------------------
// Implementation file for class : EvtRareLbToLllFF
//
// 2013-11-25 : Thomas Blake
//-----------------------------------------------------------------------------

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================

EvtRareLbToLllFF::FormFactorDependence::FormFactorDependence() :
    m_a0( 0 ), m_a2( 0 ), m_a4( 0 ), m_al( 0 ), m_ap( 0 )
{
}

EvtRareLbToLllFF::FormFactorDependence::FormFactorDependence( const double al,
                                                              const double ap ) :
    m_a0( 0 ), m_a2( 0 ), m_a4( 0 ), m_al( al ), m_ap( ap )
{
}

EvtRareLbToLllFF::FormFactorDependence::FormFactorDependence( const double a0,
                                                              const double a2,
                                                              const double a4,
                                                              const double al,
                                                              const double ap ) :
    m_a0( a0 ), m_a2( a2 ), m_a4( a4 ), m_al( al ), m_ap( ap )
{
}

EvtRareLbToLllFF::FormFactorDependence::FormFactorDependence(
    const EvtRareLbToLllFF::FormFactorDependence& other ) :
    m_a0( other.m_a0 ),
    m_a2( other.m_a2 ),
    m_a4( other.m_a4 ),
    m_al( other.m_al ),
    m_ap( other.m_ap )
{
}

EvtRareLbToLllFF::FormFactorDependence*
EvtRareLbToLllFF::FormFactorDependence::clone() const
{
    return new EvtRareLbToLllFF::FormFactorDependence( m_a0, m_a2, m_a4, m_al,
                                                       m_ap );
}

EvtRareLbToLllFF::FormFactorSet::FormFactorSet()
{
}

EvtRareLbToLllFF::FormFactorSet::FormFactorSet(
    const EvtRareLbToLllFF::FormFactorSet& other ) :
    m_F1( other.m_F1 ),
    m_F2( other.m_F2 ),
    m_F3( other.m_F3 ),
    m_F4( other.m_F4 ),
    m_G1( other.m_G1 ),
    m_G2( other.m_G2 ),
    m_G3( other.m_G3 ),
    m_G4( other.m_G4 ),
    m_H1( other.m_H1 ),
    m_H2( other.m_H2 ),
    m_H3( other.m_H3 ),
    m_H4( other.m_H4 ),
    m_H5( other.m_H5 ),
    m_H6( other.m_H6 )
{
}

void EvtRareLbToLllFF::FormFactorDependence::param( const double al,
                                                    const double ap )
{
    m_al = al;
    m_ap = ap;
}

void EvtRareLbToLllFF::FormFactorDependence::param( const double a0,
                                                    const double a2,
                                                    const double a4,
                                                    const double al,
                                                    const double ap )
{
    m_a0 = a0;
    m_a2 = a2;
    m_a4 = a4;
    m_al = al;
    m_ap = ap;
}

void EvtRareLbToLllFF::init()
{
    // Parameters for Lambda0
    auto L1115 = std::make_unique<EvtRareLbToLllFF::FormFactorSet>();
    L1115->m_F1.param( 1.21, 0.319, -0.0177, 0.387, 0.372 );
    L1115->m_F2.param( -0.202, -0.219, 0.0103, 0.387, 0.372 );
    L1115->m_F3.param( -0.0615, 0.00102, -0.00139, 0.387, 0.372 );
    L1115->m_F4.param( 0.387, 0.372 );
    L1115->m_G1.param( 0.927, 0.104, -0.00553, 0.387, 0.372 );
    L1115->m_G2.param( -0.236, -0.233, 0.0110, 0.387, 0.372 );
    L1115->m_G3.param( 0.0756, 0.0195, -0.00115, 0.387, 0.372 );
    L1115->m_G4.param( 0.387, 0.372 );
    L1115->m_H1.param( 0.936, 0.0722, -0.00643, 0.387, 0.372 );
    L1115->m_H2.param( 0.227, 0.265, -0.0101, 0.387, 0.372 );
    L1115->m_H3.param( -0.0757, -0.0195, 0.00116, 0.387, 0.372 );
    L1115->m_H4.param( -0.0174, -0.00986, -0.000524, 0.387, 0.372 );
    L1115->m_H5.param( 0.387, 0.372 );
    L1115->m_H6.param( 0.387, 0.372 );

    // Parameters for Lambda(Lambda(1520)0)
    auto L1520 = std::make_unique<EvtRareLbToLllFF::FormFactorSet>();
    L1520->m_F1.param( -1.66, -0.295, 0.00924, 0.333, 0.308 );
    L1520->m_F2.param( 0.544, 0.194, -0.00420, 0.333, 0.308 );
    L1520->m_F3.param( 0.126, 0.00799, -0.000635, 0.333, 0.308 );
    L1520->m_F4.param( -0.0330, -0.00977, 0.00211, 0.303, 0.308 );
    L1520->m_G1.param( -0.964, -0.100, 0.00264, 0.333, 0.308 );
    L1520->m_G2.param( 0.625, 0.219, -0.00508, 0.333, 0.308 );
    L1520->m_G3.param( -0.183, -0.0380, 0.00351, 0.333, 0.308 );
    L1520->m_G4.param( 0.0530, 0.0161, -0.00221, 0.333, 0.308 );
    L1520->m_H1.param( -1.08, -0.0732, 0.00464, 0.333, 0.308 );
    L1520->m_H2.param( -0.507, -0.246, 0.00309, 0.333, 0.308 );
    L1520->m_H3.param( 0.187, 0.0295, -0.00107, 0.333, 0.308 );
    L1520->m_H4.param( 0.0772, 0.0267, -0.00217, 0.333, 0.308 );
    L1520->m_H5.param( -0.0517, -0.0173, 0.00259, 0.333, 0.308 );
    L1520->m_H6.param( 0.0206, 0.00679, -0.000220, 0.333, 0.308 );

    m_FFMap[EvtPDL::getId( "Lambda0" ).getId()] = L1115.get();
    m_FFMap[EvtPDL::getId( "anti-Lambda0" ).getId()] = L1115.get();
    m_FFMap[EvtPDL::getId( "Lambda(1520)0" ).getId()] = L1520.get();
    m_FFMap[EvtPDL::getId( "anti-Lambda(1520)0" ).getId()] = L1520.get();

    m_FF = { std::move( L1115 ), std::move( L1520 ) };

    EvtGenReport( EVTGEN_INFO, "EvtGen" )
        << " EvtRareLbToLll is using form factors from arXiv:1108.6129 "
        << std::endl;
}

//=============================================================================

double EvtRareLbToLllFF::func(
    const double p, const EvtRareLbToLllFF::FormFactorDependence& dep ) const
{
    static const double mq = 0.2848;
    static const double mtilde = 1.122;

    const double asq = 0.5 * ( dep.m_al * dep.m_al + dep.m_ap * dep.m_ap );
    const double psq = p * p;

    return ( dep.m_a0 + dep.m_a2 * psq + dep.m_a4 * psq * psq ) *
           exp( -( 3. * mq * mq * psq ) / ( 2. * mtilde * mtilde * asq ) );
}

void EvtRareLbToLllFF::DiracFF( const EvtParticle& parent,
                                const EvtParticle& lambda,
                                const EvtRareLbToLllFF::FormFactorSet& dep,
                                EvtRareLbToLllFF::FormFactors& FF ) const
{
    const double M = lambda.mass();
    const double MB = parent.mass();

    const double vdotv = calculateVdotV( parent, lambda );
    const double p = lambda.getP4().d3mag();

    FF.m_F[0] = func( p, dep.m_F1 );
    FF.m_F[1] = func( p, dep.m_F2 );
    FF.m_F[2] = func( p, dep.m_F3 );

    FF.m_G[0] = func( p, dep.m_G1 );
    FF.m_G[1] = func( p, dep.m_G2 );
    FF.m_G[2] = func( p, dep.m_G3 );

    const double H1 = func( p, dep.m_H1 );
    const double H2 = func( p, dep.m_H2 );
    const double H3 = func( p, dep.m_H3 );
    const double H4 = func( p, dep.m_H4 );

    if ( isNatural( lambda ) ) {
        FF.m_FT[0] = -( MB + M ) * H1 - ( MB - M * vdotv ) * H2 -
                     ( MB * vdotv - M ) * H3;
        FF.m_FT[1] = MB * H1 + ( MB - M ) * H2 + ( MB * vdotv - M ) * H4;
        FF.m_FT[2] = M * H1 + ( MB - M ) * H3 - ( MB - M * vdotv ) * H4;

        FF.m_GT[0] = ( MB - M ) * H1 - M * ( 1. - vdotv ) * H2 -
                     MB * ( 1. - vdotv ) * H3;
        FF.m_GT[1] = MB * H1 - M * H2 - MB * H3;
        FF.m_GT[2] = M * H1 + M * H2 + MB * H3;
    } else {
        FF.m_FT[0] = ( MB - M ) * H1 - ( MB - M * vdotv ) * H2 -
                     ( MB * vdotv - M ) * H3;
        FF.m_FT[1] = MB * H1 - ( MB + M ) * H2 + ( MB * vdotv - M ) * H4;
        FF.m_FT[2] = M * H1 - ( MB + M ) * H3 - ( MB - M * vdotv ) * H4;

        FF.m_GT[0] = -( MB + M ) * H1 + M * ( 1. + vdotv ) * H2 +
                     MB * ( 1. + vdotv ) * H3;
        FF.m_GT[1] = MB * H1 - M * H2 - MB * H3;
        FF.m_GT[2] = M * H1 - M * H2 - MB * H3;
    }
}

void EvtRareLbToLllFF::RaritaSchwingerFF(
    const EvtParticle& parent, const EvtParticle& lambda,
    const EvtRareLbToLllFF::FormFactorSet& FFset,
    EvtRareLbToLllFF::FormFactors& FF ) const
{
    const double M = lambda.mass();
    const double MB = parent.mass();

    const double vdotv = calculateVdotV( parent, lambda );
    const double p = lambda.getP4().d3mag();

    FF.m_F[0] = func( p, FFset.m_F1 );
    FF.m_F[1] = func( p, FFset.m_F2 );
    FF.m_F[2] = func( p, FFset.m_F3 );
    FF.m_F[3] = func( p, FFset.m_F4 );

    FF.m_G[0] = func( p, FFset.m_G1 );
    FF.m_G[1] = func( p, FFset.m_G2 );
    FF.m_G[2] = func( p, FFset.m_G3 );
    FF.m_G[3] = func( p, FFset.m_G4 );

    const double H1 = func( p, FFset.m_H1 );
    const double H2 = func( p, FFset.m_H2 );
    const double H3 = func( p, FFset.m_H3 );
    const double H4 = func( p, FFset.m_H4 );
    const double H5 = func( p, FFset.m_H5 );
    const double H6 = func( p, FFset.m_H6 );

    if ( isNatural( lambda ) ) {
        FF.m_FT[0] = -( MB + M ) * H1 - ( MB - M * vdotv ) * H2 -
                     ( MB * vdotv - M ) * H3 - MB * H5;
        FF.m_FT[1] = MB * H1 + ( MB - M ) * H2 + ( MB * vdotv - M ) * H4 -
                     MB * H6;
        FF.m_FT[2] = M * H1 + ( MB - M ) * H3 - ( MB - M * vdotv ) * H4;
        FF.m_FT[3] = ( MB - M ) * H5 + ( MB - M * vdotv ) * H6;

        FF.m_GT[0] = ( MB - M ) * H1 - M * ( 1. - vdotv ) * H2 -
                     MB * ( 1. - vdotv ) * H3 + MB * H5 + M * H6;
        FF.m_GT[1] = MB * H1 - M * H2 - MB * H3;
        FF.m_GT[2] = M * H1 + M * H2 + MB * H3 - M * H6;
        FF.m_GT[3] = ( MB + M ) * H5 + M * ( 1. + vdotv ) * H6;
    } else {
        FF.m_FT[0] = ( MB - M ) * H1 - ( MB - M * vdotv ) * H2 -
                     ( MB * vdotv - M ) * H3 - MB * H5;
        FF.m_FT[1] = MB * H1 - ( MB + M ) * H2 + ( MB * vdotv - M ) * H4 -
                     MB * H6;
        FF.m_FT[2] = M * H1 - ( MB + M ) * H3 - ( MB - M * vdotv ) * H4;
        FF.m_FT[3] = -( MB + M ) * H5 + ( MB - M * vdotv ) * H6;

        FF.m_GT[0] = -( MB + M ) * H1 + M * ( 1. + vdotv ) * H2 +
                     MB * ( 1. + vdotv ) * H3 + MB * H5 + M * H6;
        FF.m_GT[1] = MB * H1 - M * H2 - MB * H3;
        FF.m_GT[2] = M * H1 - M * H2 - MB * H3 - M * H6;
        FF.m_GT[3] = -( MB - M ) * H5 - M * ( 1. - vdotv ) * H6;
    }
}

void EvtRareLbToLllFF::getFF( const EvtParticle& parent,
                              const EvtParticle& lambda,
                              EvtRareLbToLllFF::FormFactors& FF ) const
{
    // Find the FF information for this particle, start by setting all to zero
    FF.areZero();

    // Are the FF's for the particle known?
    auto it = m_FFMap.find( lambda.getId().getId() );

    if ( it == m_FFMap.end() ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << " EvtRareLbToLll does not contain FF for " << lambda.getId()
            << std::endl;
        return;
    }

    // Split by spin 1/2, spin 3/2
    const int spin = EvtPDL::getSpinType( lambda.getId() );

    if ( EvtSpinType::DIRAC == spin ) {
        DiracFF( parent, lambda, *( it->second ), FF );
    } else if ( spin == EvtSpinType::RARITASCHWINGER ) {
        RaritaSchwingerFF( parent, lambda, *( it->second ), FF );
    } else {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << " EvtRareLbToLll expects DIRAC or RARITASWINGER daughter "
            << std::endl;    // should add a warning here
        return;
    }

    return;
}
