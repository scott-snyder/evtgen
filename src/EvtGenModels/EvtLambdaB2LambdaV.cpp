
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

#include "EvtGenModels/EvtLambdaB2LambdaV.hh"

#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtReport.hh"

#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>

using std::fstream;
//************************************************************************
//*                                                                      *
//*                      Class EvtLambdaB2LambdaV                        *
//*                                                                      *
//************************************************************************
//DECLARE_COMPONENT( EvtLambdaB2LambdaV );

EvtLambdaB2LambdaV::EvtLambdaB2LambdaV()
{
    //set facility name
    fname = "EvtGen.EvtLambdaB2LambdaV";
}

//------------------------------------------------------------------------
// Method 'getName'
//------------------------------------------------------------------------
std::string EvtLambdaB2LambdaV::getName()
{
    return "LAMBDAB2LAMBDAV";
}

//------------------------------------------------------------------------
// Method 'clone'
//------------------------------------------------------------------------
EvtDecayBase* EvtLambdaB2LambdaV::clone()
{
    return new EvtLambdaB2LambdaV;
}

//------------------------------------------------------------------------
// Method 'initProbMax'
//------------------------------------------------------------------------
void EvtLambdaB2LambdaV::initProbMax()
{
    //maximum (case where C=0)
    double Max = 1.0 + fabs( A * B ) + 2.0 * fabs( A * abs( C ) );
    if ( verbose() ) {
        EvtGenReport( EVTGEN_DEBUG, fname.c_str() )
            << " PDF max value : " << Max << std::endl;
    }
    setProbMax( Max );
}

//------------------------------------------------------------------------
// Method 'init'
//------------------------------------------------------------------------
void EvtLambdaB2LambdaV::init()
{
    bool antiparticle = false;

    //introduction
    if ( verbose() ) {
        EvtGenReport( EVTGEN_DEBUG, fname.c_str() )
            << "*************************************************" << std::endl;
        EvtGenReport( EVTGEN_DEBUG, fname.c_str() )
            << "*    Event Model Class : EvtLambdaB2LambdaV     *" << std::endl;
        EvtGenReport( EVTGEN_DEBUG, fname.c_str() )
            << "*************************************************" << std::endl;
    }

    //check the number of arguments
    checkNArg( 2 );

    //check the number of daughters
    checkNDaug( 2 );

    const EvtId Id_mother = getParentId();
    const EvtId Id_daug1 = getDaug( 0 );
    const EvtId Id_daug2 = getDaug( 1 );

    //lambdab identification
    if ( Id_mother == EvtPDL::getId( "Lambda_b0" ) ) {
        antiparticle = false;
    } else if ( Id_mother == EvtPDL::getId( "anti-Lambda_b0" ) ) {
        antiparticle = true;
    } else {
        EvtGenReport( EVTGEN_ERROR, fname.c_str() )
            << " Mother is not a Lambda_b0 or an anti-Lambda_b0, but a "
            << EvtPDL::name( Id_mother ) << std::endl;
        abort();
    }

    //lambda
    if ( !( Id_daug1 == EvtPDL::getId( "Lambda0" ) && !antiparticle ) &&
         !( Id_daug1 == EvtPDL::getId( "anti-Lambda0" ) && antiparticle ) ) {
        if ( !antiparticle ) {
            EvtGenReport( EVTGEN_ERROR, fname.c_str() )
                << " Daughter1 is not a Lambda0, but a "
                << EvtPDL::name( Id_daug1 ) << std::endl;
        } else {
            EvtGenReport( EVTGEN_ERROR, fname.c_str() )
                << " Daughter1 is not an anti-Lambda0, but a "
                << EvtPDL::name( Id_daug1 ) << std::endl;
        }
        abort();
    }

    //identification meson V
    if ( getArg( 1 ) == 1 )
        Vtype = VID::JPSI;
    else if ( getArg( 1 ) == 2 )
        Vtype = VID::RHO;
    else if ( getArg( 1 ) == 3 )
        Vtype = VID::OMEGA;
    else if ( getArg( 1 ) == 4 )
        Vtype = VID::RHO_OMEGA_MIXING;
    else {
        EvtGenReport( EVTGEN_ERROR, fname.c_str() )
            << " Vtype " << getArg( 1 ) << " is unknown" << std::endl;
        abort();
    }

    //check vector meson V
    if ( Id_daug2 == EvtPDL::getId( "J/psi" ) && Vtype == VID::JPSI ) {
        if ( !antiparticle )
            EvtGenReport( EVTGEN_DEBUG, fname.c_str() )
                << " Decay mode successfully initialized : Lambda_b0 -> Lambda J/psi"
                << std::endl;
        else
            EvtGenReport( EVTGEN_DEBUG, fname.c_str() )
                << " Decay mode successfully initialized : anti-Lambda_b0 -> anti-Lambda J/psi"
                << std::endl;
    } else if ( Id_daug2 == EvtPDL::getId( "rho0" ) && Vtype == VID::RHO ) {
        if ( !antiparticle )
            EvtGenReport( EVTGEN_DEBUG, fname.c_str() )
                << " Decay mode successfully initialized : Lambda_b0 -> Lambda rho0"
                << std::endl;
        else
            EvtGenReport( EVTGEN_DEBUG, fname.c_str() )
                << " Decay mode successfully initialized : anti-Lambda_b0 -> anti-Lambda rho0"
                << std::endl;
    } else if ( Id_daug2 == EvtPDL::getId( "omega" ) && Vtype == VID::OMEGA ) {
        if ( !antiparticle )
            EvtGenReport( EVTGEN_DEBUG, fname.c_str() )
                << " Decay mode successfully initialized : Lambda_b0 -> Lambda omega"
                << std::endl;
        else
            EvtGenReport( EVTGEN_DEBUG, fname.c_str() )
                << " Decay mode successfully initialized : anti-Lambda_b0 -> anti-Lambda omega"
                << std::endl;
    } else if ( ( Id_daug2 == EvtPDL::getId( "omega" ) ||
                  Id_daug2 == EvtPDL::getId( "rho0" ) ) &&
                Vtype == VID::RHO_OMEGA_MIXING ) {
        if ( !antiparticle )
            EvtGenReport( EVTGEN_DEBUG, fname.c_str() )
                << " Decay mode successfully initialized : "
                << "Lambda_b0 -> Lambda rho-omega-mixing" << std::endl;
        else
            EvtGenReport( EVTGEN_DEBUG, fname.c_str() )
                << " Decay mode successfully initialized : "
                << "anti-Lambda_b0 -> anti-Lambda rho-omega-mixing" << std::endl;
    }

    else {
        EvtGenReport( EVTGEN_ERROR, fname.c_str() )
            << " Daughter2 is not a J/psi, phi or rho0 but a "
            << EvtPDL::name( Id_daug2 ) << std::endl;
        abort();
    }

    //fix dynamics parameters
    B = (double)getArg( 0 );
    C = EvtComplex( ( sqrt( 2. ) / 2. ), ( sqrt( 2. ) / 2. ) );
    switch ( Vtype ) {
        case VID::JPSI:
            A = 0.490;
            break;
        case VID::RHO:
        case VID::OMEGA:
        case VID::RHO_OMEGA_MIXING:
            A = 0.194;
            break;
        default:
            A = 0;
            break;
    }

    if ( verbose() ) {
        EvtGenReport( EVTGEN_DEBUG, fname.c_str() )
            << " LambdaB decay parameters : " << std::endl;
        EvtGenReport( EVTGEN_DEBUG, fname.c_str() )
            << "   - lambda asymmetry A = " << A << std::endl;
        EvtGenReport( EVTGEN_DEBUG, fname.c_str() )
            << "   - lambdab polarisation B = " << B << std::endl;
        EvtGenReport( EVTGEN_DEBUG, fname.c_str() )
            << "   - lambdab density matrix rho+- C = " << C << std::endl;
    }
}

//------------------------------------------------------------------------
// Method 'decay'
//------------------------------------------------------------------------
void EvtLambdaB2LambdaV::decay( EvtParticle* lambdab )
{
    // Using formulae from hep-ph/0602043

    lambdab->makeDaughters( getNDaug(), getDaugs() );

    //get lambda and V particles
    EvtParticle* lambda = lambdab->getDaug( 0 );
    EvtParticle* V = lambdab->getDaug( 1 );

    //get resonance masses
    // - LAMBDAB          -> mass given by the preceding class
    // - LAMBDA           -> nominal mass
    // - V                -> getVmass method
    double MASS_LAMBDAB = lambdab->mass();
    double MASS_LAMBDA = EvtPDL::getMeanMass( EvtPDL::getId( "Lambda0" ) );
    double MASS_V = getVMass( MASS_LAMBDAB, MASS_LAMBDA );

    //generate random angles
    double phi = EvtRandom::Flat( 0, 2 * EvtConst::pi );
    double theta = acos( EvtRandom::Flat( -1, +1 ) );

    // Get 4-vectors
    double E_lambda = ( MASS_LAMBDAB * MASS_LAMBDAB +
                        MASS_LAMBDA * MASS_LAMBDA - MASS_V * MASS_V ) /
                      ( 2 * MASS_LAMBDAB );
    double E_V = ( MASS_LAMBDAB * MASS_LAMBDAB + MASS_V * MASS_V -
                   MASS_LAMBDA * MASS_LAMBDA ) /
                 ( 2 * MASS_LAMBDAB );
    double P = sqrt( E_lambda * E_lambda - lambda->mass() * lambda->mass() );

    EvtVector4R P_lambdab = lambdab->getP4();

    double px = P_lambdab.get( 1 );
    double py = P_lambdab.get( 2 );
    double pz = P_lambdab.get( 3 );
    double E = P_lambdab.get( 0 );

    const double pxSq = px * px;
    const double pySq = py * py;
    const double pT = sqrt( pxSq + pySq );
    const double invPT = pT > 0.0 ? 1.0 / pT : 0.0;
    EvtVector4R q_lambdab2( E, pT, 0.0, pz );

    EvtVector4R q_lambdab3( E, q_lambdab2.get( 3 ), q_lambdab2.get( 1 ),
                            q_lambdab2.get( 2 ) );

    EvtVector4R q_lambda0( E_lambda, P * sin( theta ) * cos( phi ),
                           P * sin( theta ) * sin( phi ), P * cos( theta ) );

    EvtVector4R q_V0( E_V, -P * sin( theta ) * cos( phi ),
                      -P * sin( theta ) * sin( phi ), -P * cos( theta ) );

    EvtVector4R q_lambda1( E_lambda, q_lambda0.get( 2 ), q_lambda0.get( 3 ),
                           q_lambda0.get( 1 ) );

    EvtVector4R q_V1( E_V, q_V0.get( 2 ), q_V0.get( 3 ), q_V0.get( 1 ) );

    EvtVector4R q_lambda(
        E_lambda, invPT * ( px * q_lambda1.get( 1 ) - py * q_lambda1.get( 2 ) ),
        invPT * ( py * q_lambda1.get( 1 ) + px * q_lambda1.get( 2 ) ),
        q_lambda1.get( 3 ) );

    EvtVector4R q_V( E_V, invPT * ( px * q_V1.get( 1 ) - py * q_V1.get( 2 ) ),
                     invPT * ( py * q_V1.get( 1 ) + px * q_V1.get( 2 ) ),
                     q_V1.get( 3 ) );

    // Set particle 4-vectors
    lambda->init( getDaugs()[0], q_lambda );
    V->init( getDaugs()[1], q_V );

    // Find pdf (eq 11 in paper)
    double pdf = 1.0 + A * B * cos( theta ) +
                 2.0 * A * real( C * EvtComplex( cos( phi ), sin( phi ) ) ) *
                     sin( theta );

    //set probability
    setProb( pdf );

    if ( verbose() ) {
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Angular angles  : theta = " << theta << " ; phi = " << phi
            << std::endl;

        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << "E of lambdab:  " << P_lambdab.get( 0 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << "E of lambdab:  " << E << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " LambdaB  px: " << px << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " LambdaB  py: " << py << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " LambdaB  pz: " << pz << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " LambdaB  E: " << E << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Lambdab3  E:  " << q_lambdab3.get( 0 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Lambda 0 px:  " << q_lambda0.get( 1 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Lambda 0 py:  " << q_lambda0.get( 2 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Lambda 0 pz:  " << q_lambda0.get( 3 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Lambda 0 E:  " << q_lambda0.get( 0 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Lambda 1 px:  " << q_lambda1.get( 1 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Lambda 1 py:  " << q_lambda1.get( 2 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Lambda 1 pz:  " << q_lambda1.get( 3 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Lambda 1 E:  " << q_lambda1.get( 0 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Lambda  px:  " << q_lambda.get( 1 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Lambda  py:  " << q_lambda.get( 2 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Lambda  pz:  " << q_lambda.get( 3 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Lambda  E:  " << q_lambda0.get( 3 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " V 0 px:  " << q_V0.get( 1 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " V 0 py:  " << q_V0.get( 2 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " V 0 pz:  " << q_V0.get( 3 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " V 0 E:  " << q_V0.get( 0 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " V 1 px:  " << q_V1.get( 1 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " V 1 py:  " << q_V1.get( 2 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " V 1 pz:  " << q_V1.get( 3 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " V 1 E:  " << q_V1.get( 0 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " V  px:  " << q_V.get( 1 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " V  py:  " << q_V.get( 2 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " V  pz:  " << q_V.get( 3 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " V  E:  " << q_V0.get( 3 ) << std::endl;
        EvtGenReport( EVTGEN_DEBUG, fname.c_str() )
            << " LambdaB decay pdf value : " << pdf << std::endl;
    }
}

//------------------------------------------------------------------------
// Method 'BreitWignerRelPDF'
//------------------------------------------------------------------------
double EvtLambdaB2LambdaV::BreitWignerRelPDF( double m, double _m0, double _g0 )
{
    double a = ( _m0 * _g0 ) * ( _m0 * _g0 );
    double b = ( m * m - _m0 * _m0 ) * ( m * m - _m0 * _m0 );
    return a / ( b + a );
}

//------------------------------------------------------------------------
// Method 'RhoOmegaMixingPDF'
//------------------------------------------------------------------------
double EvtLambdaB2LambdaV::RhoOmegaMixingPDF( double m, double _mr, double _gr,
                                              double _mo, double _go )
{
    double a = m * m - _mr * _mr;
    double b = m * m - _mo * _mo;
    double c = _gr * _mr;
    double d = _go * _mo;
    double re_pi = -3500e-6;    //GeV^2
    double im_pi = -300e-6;     //GeV^2
    double va_pi = re_pi + im_pi;

    //compute pdf value
    double f = 1 / ( a * a + c * c ) *
               ( 1 + ( va_pi * va_pi + 2 * b * re_pi + 2 * d * im_pi ) /
                         ( b * b + d * d ) );

    //compute pdf max value
    a = 0;
    b = _mr * _mr - _mo * _mo;

    double maxi = 1 / ( a * a + c * c ) *
                  ( 1 + ( va_pi * va_pi + 2 * b * re_pi + 2 * d * im_pi ) /
                            ( b * b + d * d ) );

    return f / maxi;
}

//------------------------------------------------------------------------
// Method 'getVMass'
//------------------------------------------------------------------------
double EvtLambdaB2LambdaV::getVMass( double MASS_LAMBDAB, double MASS_LAMBDA )
{
    //JPSI case
    if ( Vtype == VID::JPSI ) {
        return EvtPDL::getMass( EvtPDL::getId( "J/psi" ) );
    }

    //RHO,OMEGA,MIXING case
    else {
        //parameters
        double MASS_RHO = EvtPDL::getMeanMass( EvtPDL::getId( "rho0" ) );
        double MASS_OMEGA = EvtPDL::getMeanMass( EvtPDL::getId( "omega" ) );
        double WIDTH_RHO = EvtPDL::getWidth( EvtPDL::getId( "rho0" ) );
        double WIDTH_OMEGA = EvtPDL::getWidth( EvtPDL::getId( "omega" ) );
        double MASS_PION = EvtPDL::getMeanMass( EvtPDL::getId( "pi-" ) );
        double _max = MASS_LAMBDAB - MASS_LAMBDA;
        double _min = 2 * MASS_PION;

        double mass = 0;
        bool test = false;
        int ntimes = 10000;
        do {
            double y = EvtRandom::Flat( 0, 1 );

            //generate mass
            mass = ( _max - _min ) * EvtRandom::Flat( 0, 1 ) + _min;

            //pdf calculate
            double f = 0;
            switch ( Vtype ) {
                case VID::RHO:
                    f = BreitWignerRelPDF( mass, MASS_RHO, WIDTH_RHO );
                    break;
                case VID::OMEGA:
                    f = BreitWignerRelPDF( mass, MASS_OMEGA, WIDTH_OMEGA );
                    break;
                case VID::RHO_OMEGA_MIXING:
                    f = RhoOmegaMixingPDF( mass, MASS_RHO, WIDTH_RHO,
                                           MASS_OMEGA, WIDTH_OMEGA );
                    break;
                default:
                    f = 1;
                    break;
            }

            ntimes--;
            if ( y < f )
                test = true;
        } while ( ntimes && !test );

        //looping 10000 times
        if ( ntimes == 0 ) {
            EvtGenReport( EVTGEN_INFO, fname.c_str() )
                << "Tried accept/reject:10000"
                << " times, and rejected all the times!" << std::endl;
            EvtGenReport( EVTGEN_INFO, fname.c_str() )
                << "Is therefore accepting the last event!" << std::endl;
        }

        //return V particle mass
        return mass;
    }
}

//************************************************************************
//*                                                                      *
//*                Class EvtLambda2PPiForLambdaB2LambdaV                 *
//*                                                                      *
//************************************************************************

//------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------
EvtLambda2PPiForLambdaB2LambdaV::EvtLambda2PPiForLambdaB2LambdaV()
{
    //set facility name
    fname = "EvtGen.EvtLambda2PPiForLambdaB2LambdaV";
}

//------------------------------------------------------------------------
// Method 'getName'
//------------------------------------------------------------------------
std::string EvtLambda2PPiForLambdaB2LambdaV::getName()
{
    return "LAMBDA2PPIFORLAMBDAB2LAMBDAV";
}

//------------------------------------------------------------------------
// Method 'clone'
//------------------------------------------------------------------------
EvtDecayBase* EvtLambda2PPiForLambdaB2LambdaV::clone()
{
    return new EvtLambda2PPiForLambdaB2LambdaV;
}

//------------------------------------------------------------------------
// Method 'initProbMax'
//------------------------------------------------------------------------
void EvtLambda2PPiForLambdaB2LambdaV::initProbMax()
{
    double Max = 1.0 + fabs( A * B ) + 2.0 * fabs( A * abs( D ) );
    if ( verbose() ) {
        EvtGenReport( EVTGEN_DEBUG, fname.c_str() )
            << " PDF max value : " << Max << std::endl;
    }
    setProbMax( Max );
}

//------------------------------------------------------------------------
// Method 'init'
//------------------------------------------------------------------------
void EvtLambda2PPiForLambdaB2LambdaV::init()
{
    bool antiparticle = false;

    if ( verbose() ) {
        //introduction
        EvtGenReport( EVTGEN_DEBUG, fname.c_str() )
            << " ***********************************************************"
            << std::endl;
        EvtGenReport( EVTGEN_DEBUG, fname.c_str() )
            << " *   Event Model Class : EvtLambda2PPiForLambdaB2LambdaV   *"
            << std::endl;
        EvtGenReport( EVTGEN_DEBUG, fname.c_str() )
            << " ***********************************************************"
            << std::endl;
    }

    //check the number of arguments
    checkNArg( 2 );

    //check the number of daughters
    checkNDaug( 2 );

    const EvtId Id_mother = getParentId();
    const EvtId Id_daug1 = getDaug( 0 );
    const EvtId Id_daug2 = getDaug( 1 );

    //lambda  identification
    if ( Id_mother == EvtPDL::getId( "Lambda0" ) ) {
        antiparticle = false;
    } else if ( Id_mother == EvtPDL::getId( "anti-Lambda0" ) ) {
        antiparticle = true;
    } else {
        EvtGenReport( EVTGEN_ERROR, fname.c_str() )
            << " Mother is not a Lambda0 or an anti-Lambda0, but a "
            << EvtPDL::name( Id_mother ) << std::endl;
        abort();
    }

    //proton identification
    if ( !( Id_daug1 == EvtPDL::getId( "p+" ) && !antiparticle ) &&
         !( Id_daug1 == EvtPDL::getId( "anti-p-" ) && antiparticle ) ) {
        if ( !antiparticle ) {
            EvtGenReport( EVTGEN_ERROR, fname.c_str() )
                << " Daughter1 is not a p+, but a " << EvtPDL::name( Id_daug1 )
                << std::endl;
        } else {
            EvtGenReport( EVTGEN_ERROR, fname.c_str() )
                << " Daughter1 is not an anti-p-, but a "
                << EvtPDL::name( Id_daug1 ) << std::endl;
        }
        abort();
    }

    //pion identification
    if ( !( Id_daug2 == EvtPDL::getId( "pi-" ) && !antiparticle ) &&
         !( Id_daug2 == EvtPDL::getId( "pi+" ) && antiparticle ) ) {
        if ( !antiparticle ) {
            EvtGenReport( EVTGEN_ERROR, fname.c_str() )
                << " Daughter2 is not a p-, but a " << EvtPDL::name( Id_daug1 )
                << std::endl;
        } else {
            EvtGenReport( EVTGEN_ERROR, fname.c_str() )
                << " Daughter2 is not an p+, but a " << EvtPDL::name( Id_daug1 )
                << std::endl;
        }
        abort();
    }
    if ( !antiparticle )
        EvtGenReport( EVTGEN_DEBUG, fname.c_str() )
            << " Decay mode successfully initialized : Lambda0 -> p+ pi-"
            << std::endl;
    else
        EvtGenReport( EVTGEN_DEBUG, fname.c_str() )
            << " Decay mode successfully initialized : Anti-Lambda0 -> anti-p- pi+"
            << std::endl;

    //identification meson V
    if ( getArg( 1 ) == 1 ) {
        Vtype = VID::JPSI;
        if ( !antiparticle )
            EvtGenReport( EVTGEN_DEBUG, fname.c_str() )
                << " From : Lambda_b0 -> Lambda J/psi" << std::endl;
        else
            EvtGenReport( EVTGEN_DEBUG, fname.c_str() )
                << " From : anti-Lambda_b0 -> anti-Lambda J/psi" << std::endl;
    } else if ( getArg( 1 ) == 2 ) {
        Vtype = VID::RHO;
        if ( !antiparticle )
            EvtGenReport( EVTGEN_DEBUG, fname.c_str() )
                << " From : Lambda_b0 -> Lambda rho0" << std::endl;
        else
            EvtGenReport( EVTGEN_DEBUG, fname.c_str() )
                << " From : anti-Lambda_b0 -> anti-Lambda rho0" << std::endl;
    } else if ( getArg( 1 ) == 3 ) {
        Vtype = VID::OMEGA;
        if ( !antiparticle )
            EvtGenReport( EVTGEN_DEBUG, fname.c_str() )
                << " From : Lambda_b0 -> Lambda omega" << std::endl;
        else
            EvtGenReport( EVTGEN_DEBUG, fname.c_str() )
                << " From : anti-Lambda_b0 -> anti-Lambda omega" << std::endl;
    } else if ( getArg( 1 ) == 4 ) {
        Vtype = VID::RHO_OMEGA_MIXING;
    } else {
        EvtGenReport( EVTGEN_ERROR, fname.c_str() )
            << " Vtype " << getArg( 1 ) << " is unknown" << std::endl;
        if ( !antiparticle ) {
            EvtGenReport( EVTGEN_DEBUG, fname.c_str() )
                << " From : Lambda_b0 -> Lambda rho-omega-mixing" << std::endl;
        } else {
            EvtGenReport( EVTGEN_DEBUG, fname.c_str() )
                << " From : anti-Lambda_b0 -> anti-Lambda rho-omega-mixing"
                << std::endl;
        }
        abort();
    }

    //constants
    A = 0.642;
    C = (double)getArg( 0 );
    switch ( Vtype ) {
        case VID::JPSI:
            B = -0.167;
            D = EvtComplex( 0.25, 0.0 );
            break;
        case VID::RHO:
        case VID::OMEGA:
        case VID::RHO_OMEGA_MIXING:
            B = -0.21;
            D = EvtComplex( 0.31, 0.0 );
            break;
        default:
            B = 0;
            D = EvtComplex( 0, 0 );
            break;
    }

    if ( verbose() ) {
        EvtGenReport( EVTGEN_DEBUG, fname.c_str() )
            << " Lambda decay parameters : " << std::endl;
        EvtGenReport( EVTGEN_DEBUG, fname.c_str() )
            << "   - proton asymmetry A = " << A << std::endl;
        EvtGenReport( EVTGEN_DEBUG, fname.c_str() )
            << "   - lambda polarisation B = " << B << std::endl;
        EvtGenReport( EVTGEN_DEBUG, fname.c_str() )
            << "   - lambdaB polarisation C = " << C << std::endl;
        EvtGenReport( EVTGEN_DEBUG, fname.c_str() )
            << "   - lambda density matrix rho+- D = " << D << std::endl;
    }

    // Update max prob
    initProbMax();
}

//------------------------------------------------------------------------
// Method 'decay'
//------------------------------------------------------------------------
void EvtLambda2PPiForLambdaB2LambdaV::decay( EvtParticle* lambda )
{
    lambda->makeDaughters( getNDaug(), getDaugs() );

    //get proton and pion particles
    EvtParticle* proton = lambda->getDaug( 0 );
    EvtParticle* pion = lambda->getDaug( 1 );

    //get resonance masses
    // - LAMBDA        -> mass given by EvtLambdaB2LambdaV class
    // - PROTON & PION -> nominal mass
    double MASS_LAMBDA = lambda->mass();
    double MASS_PROTON = EvtPDL::getMeanMass( EvtPDL::getId( "p+" ) );
    double MASS_PION = EvtPDL::getMeanMass( EvtPDL::getId( "pi-" ) );

    //generate random angles
    double phi = EvtRandom::Flat( 0.0, 2.0 * EvtConst::pi );
    double theta = acos( EvtRandom::Flat( -1.0, 1.0 ) );

    // Find 4-vectors
    double E_proton = ( MASS_LAMBDA * MASS_LAMBDA + MASS_PROTON * MASS_PROTON -
                        MASS_PION * MASS_PION ) /
                      ( 2 * MASS_LAMBDA );
    double E_pion = ( MASS_LAMBDA * MASS_LAMBDA + MASS_PION * MASS_PION -
                      MASS_PROTON * MASS_PROTON ) /
                    ( 2 * MASS_LAMBDA );
    double P = sqrt( E_proton * E_proton - proton->mass() * proton->mass() );

    EvtVector4R P_lambda = lambda->getP4();
    EvtParticle* Mother_lambda = lambda->getParent();
    EvtVector4R lambdab = Mother_lambda->getP4();

    double E_lambda = P_lambda.get( 0 );
    double px_M = lambdab.get( 1 );
    double py_M = lambdab.get( 2 );
    double pz_M = lambdab.get( 3 );
    double E_M = lambdab.get( 0 );

    const double px_MSq = px_M * px_M;
    const double py_MSq = py_M * py_M;
    const double pT_M = sqrt( px_MSq + py_MSq );
    const double invPT_M = pT_M > 0.0 ? 1.0 / pT_M : 0.0;

    EvtVector4R q_lambdab2( E_M, pT_M, 0.0, pz_M );

    EvtVector4R q_lambdab3( E_M, q_lambdab2.get( 3 ), q_lambdab2.get( 1 ),
                            q_lambdab2.get( 2 ) );

    EvtVector4R q_lambda1(
        E_lambda,
        invPT_M * ( px_M * P_lambda.get( 1 ) + py_M * P_lambda.get( 2 ) ),
        invPT_M * ( -py_M * P_lambda.get( 1 ) + px_M * P_lambda.get( 2 ) ),
        P_lambda.get( 3 ) );

    EvtVector4R q_lambda2( E_lambda, q_lambda1.get( 3 ), q_lambda1.get( 1 ),
                           q_lambda1.get( 2 ) );

    const double E = q_lambda2.get( 0 );
    const double px = q_lambda2.get( 1 );
    const double py = q_lambda2.get( 2 );
    const double pz = q_lambda2.get( 3 );
    const double pxSq = px * px;
    const double pySq = py * py;
    const double pzSq = pz * pz;
    const double pTSq = pxSq + pySq;
    const double p = sqrt( pTSq + pzSq );
    const double invP = p > 0.0 ? 1.0 / p : 0.0;
    const double invPSq = invP * invP;
    const double pT = sqrt( pTSq );

    //EvtVector4R q_lambda4( E, 0.0, 0.0, p );
    EvtVector4R q_lambda4( E, pT, 0.0, pz );

    EvtVector4R q_proton1( E_proton, P * sin( theta ) * cos( phi ),
                           P * sin( theta ) * sin( phi ), P * cos( theta ) );
    EvtVector4R q_pion1( E_pion, -P * sin( theta ) * cos( phi ),
                         -P * sin( theta ) * sin( phi ), -P * cos( theta ) );

    const double qx_proton1 = q_proton1.get( 1 );
    const double qy_proton1 = q_proton1.get( 2 );
    const double qz_proton1 = q_proton1.get( 3 );

    EvtVector4R q_proton3( E_proton,
                           invPSq * ( qx_proton1 * px * pz - qy_proton1 * py * p +
                                      qz_proton1 * pT * px ),
                           invPSq * ( qx_proton1 * py * pz + qy_proton1 * px * p +
                                      qz_proton1 * pT * py ),
                           invP * ( -qx_proton1 * pT + qz_proton1 * pz ) );

    const double qx_pion1 = q_pion1.get( 1 );
    const double qy_pion1 = q_pion1.get( 2 );
    const double qz_pion1 = q_pion1.get( 3 );

    EvtVector4R q_pion3(
        E_pion,
        invPSq * ( qx_pion1 * px * pz - qy_pion1 * py * p + qz_pion1 * pT * px ),
        invPSq * ( qx_pion1 * py * pz + qy_pion1 * px * p + qz_pion1 * pT * py ),
        invP * ( -qx_pion1 * pT + qz_pion1 * pz ) );

    EvtVector4R q_proton5( q_proton3.get( 0 ), q_proton3.get( 2 ),
                           q_proton3.get( 3 ), q_proton3.get( 1 ) );

    EvtVector4R q_pion5( q_pion3.get( 0 ), q_pion3.get( 2 ), q_pion3.get( 3 ),
                         q_pion3.get( 1 ) );

    EvtVector4R q_proton(
        q_proton5.get( 0 ),
        invPT_M * ( px_M * q_proton5.get( 1 ) - py_M * q_proton5.get( 2 ) ),
        invPT_M * ( py_M * q_proton5.get( 1 ) + px_M * q_proton5.get( 2 ) ),
        q_proton5.get( 3 ) );

    EvtVector4R q_pion(
        q_pion5.get( 0 ),
        invPT_M * ( px_M * q_pion5.get( 1 ) - py_M * q_pion5.get( 2 ) ),
        invPT_M * ( py_M * q_pion5.get( 1 ) + px_M * q_pion5.get( 2 ) ),
        q_pion5.get( 3 ) );

    // Set particle 4 vectors
    proton->init( getDaugs()[0], q_proton );
    pion->init( getDaugs()[1], q_pion );

    //compute pdf
    double pdf = 1.0 + A * B * cos( theta ) +
                 2.0 * A * real( D * EvtComplex( cos( phi ), sin( phi ) ) ) *
                     sin( theta );
    //set probability
    setProb( pdf );

    if ( verbose() ) {
        EvtGenReport( EVTGEN_DEBUG, fname.c_str() )
            << " Angular angles  : theta = " << theta << " ; phi = " << phi
            << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Lambdab  px: " << px_M << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Lambdab  py: " << py_M << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Lambdab  pz: " << pz_M << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Lambdab  E: " << E_M << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Lambdab2  px:  " << q_lambdab2.get( 1 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Lambdab2  py:  " << q_lambdab2.get( 2 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Lambdab2  pz:  " << q_lambdab2.get( 3 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Lambdab2  E:   " << q_lambdab2.get( 0 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Lambdab3  px:  " << q_lambdab3.get( 1 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Lambdab3  py:  " << q_lambdab3.get( 2 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Lambdab3  pz:  " << q_lambdab3.get( 3 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Lambdab3  E:   " << q_lambdab3.get( 0 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Lambda 0  px:  " << P_lambda.get( 1 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Lambda 0  py:  " << P_lambda.get( 2 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Lambda 0  pz:  " << P_lambda.get( 3 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Lambda 0  E:   " << P_lambda.get( 0 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Lambda 1  px:  " << q_lambda1.get( 1 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Lambda 1  py:  " << q_lambda1.get( 2 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Lambda 1  pz:  " << q_lambda1.get( 3 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Lambda 1  E:   " << q_lambda1.get( 0 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Lambda 2  px:  " << q_lambda2.get( 1 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Lambda 2  py:  " << q_lambda2.get( 2 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Lambda 2  pz:  " << q_lambda2.get( 3 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Lambda 2  E:   " << q_lambda2.get( 0 ) << std::endl;

        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Lambda  px: " << px << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Lambda  py: " << py << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Lambda  pz: " << pz << std::endl;

        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " pion 1 px:  " << q_pion1.get( 1 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " pion 1 py:  " << q_pion1.get( 2 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " pion 1 pz:  " << q_pion1.get( 3 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " pion 1 E:   " << q_pion1.get( 0 ) << std::endl;

        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " pion 3 px:  " << q_pion3.get( 1 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " pion 3 px:  " << q_pion3.get( 1 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " pion 3 py:  " << q_pion3.get( 2 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " pion 3 pz:  " << q_pion3.get( 3 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " pion 3 E:   " << q_pion3.get( 0 ) << std::endl;

        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " pion 5 px:  " << q_pion5.get( 1 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " pion 5 py:  " << q_pion5.get( 2 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " pion 5 pz:  " << q_pion5.get( 3 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " pion 5 E:   " << q_pion5.get( 0 ) << std::endl;

        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " proton 1  px:  " << q_proton1.get( 1 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " proton 1  py:  " << q_proton1.get( 2 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " proton 1  pz:  " << q_proton1.get( 3 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " proton 1  E:   " << q_proton1.get( 0 ) << std::endl;

        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " proton 3  px:  " << q_proton3.get( 1 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " proton 3  py:  " << q_proton3.get( 2 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " proton 3  pz:  " << q_proton3.get( 3 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " proton 3  E:   " << q_proton3.get( 0 ) << std::endl;

        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " proton 5  px:  " << q_proton5.get( 1 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " proton 5  py:  " << q_proton5.get( 2 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " proton 5  pz:  " << q_proton5.get( 3 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " proton 5  E:   " << q_proton5.get( 0 ) << std::endl;

        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " proton  px:  " << q_proton.get( 1 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " proton  py:  " << q_proton.get( 2 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << "proton  pz:  " << q_proton.get( 3 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " pion px:  " << q_pion.get( 1 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " pion py:  " << q_pion.get( 2 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " pion pz:  " << q_pion.get( 3 ) << std::endl;

        EvtGenReport( EVTGEN_DEBUG, fname.c_str() )
            << " Lambda decay pdf value : " << pdf << std::endl;
    }
}

//************************************************************************
//*                                                                      *
//*                Class EvtV2VpVmForLambdaB2LambdaV                     *
//*                                                                      *
//************************************************************************

//------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------
EvtV2VpVmForLambdaB2LambdaV::EvtV2VpVmForLambdaB2LambdaV()
{
    //set facility name
    fname = "EvtGen.EvtV2VpVmForLambdaB2LambdaV";
}

//------------------------------------------------------------------------
// Method 'getName'
//------------------------------------------------------------------------
std::string EvtV2VpVmForLambdaB2LambdaV::getName()
{
    return "V2VPVMFORLAMBDAB2LAMBDAV";
}

//------------------------------------------------------------------------
// Method 'clone'
//------------------------------------------------------------------------
EvtDecayBase* EvtV2VpVmForLambdaB2LambdaV::clone()
{
    return new EvtV2VpVmForLambdaB2LambdaV;
}

//------------------------------------------------------------------------
// Method 'initProbMax'
//------------------------------------------------------------------------
void EvtV2VpVmForLambdaB2LambdaV::initProbMax()
{
    //maximum
    double Max = 0.0;
    if ( Vtype == VID::JPSI ) {
        if ( ( 1.0 - 3.0 * A ) > 0.0 )
            Max = 2.0 * ( 1.0 - A );
        else
            Max = 1.0 + A;
    } else {
        if ( ( 3.0 * A - 1.0 ) >= 0.0 )
            Max = 2.0 * A;
        else
            Max = 1.0 - A;
    }

    if ( verbose() ) {
        EvtGenReport( EVTGEN_DEBUG, fname.c_str() )
            << " PDF max value : " << Max << std::endl;
    }
    setProbMax( Max );
}

//------------------------------------------------------------------------
// Method 'init'
//------------------------------------------------------------------------
void EvtV2VpVmForLambdaB2LambdaV::init()
{
    if ( verbose() ) {
        //introduction
        EvtGenReport( EVTGEN_DEBUG, fname.c_str() )
            << " ***********************************************************"
            << std::endl;
        EvtGenReport( EVTGEN_DEBUG, fname.c_str() )
            << " *     Event Model Class : EvtV2VpVmForLambdaB2LambdaV     *"
            << std::endl;
        EvtGenReport( EVTGEN_DEBUG, fname.c_str() )
            << " ***********************************************************"
            << std::endl;
    }

    //check the number of arguments
    checkNArg( 2 );

    //check the number of daughters
    checkNDaug( 2 );

    const EvtId Id_mother = getParentId();
    const EvtId Id_daug1 = getDaug( 0 );
    const EvtId Id_daug2 = getDaug( 1 );

    //identification meson V
    if ( getArg( 1 ) == 1 )
        Vtype = VID::JPSI;
    else if ( getArg( 1 ) == 2 )
        Vtype = VID::RHO;
    else if ( getArg( 1 ) == 3 )
        Vtype = VID::OMEGA;
    else if ( getArg( 1 ) == 4 )
        Vtype = VID::RHO_OMEGA_MIXING;
    else {
        EvtGenReport( EVTGEN_ERROR, fname.c_str() )
            << " Vtype " << getArg( 1 ) << " is unknown" << std::endl;
        abort();
    }

    //vector meson V
    if ( Id_mother == EvtPDL::getId( "J/psi" ) && Vtype == VID::JPSI ) {
    } else if ( Id_mother == EvtPDL::getId( "omega" ) && Vtype == VID::OMEGA ) {
    } else if ( Id_mother == EvtPDL::getId( "rho0" ) && Vtype == VID::RHO ) {
    } else if ( ( Id_mother == EvtPDL::getId( "rho0" ) ||
                  Id_mother == EvtPDL::getId( "omega" ) ) &&
                Vtype == VID::RHO_OMEGA_MIXING ) {
    } else {
        EvtGenReport( EVTGEN_ERROR, fname.c_str() )
            << " Mother is not a J/psi, phi or rho0 but a "
            << EvtPDL::name( Id_mother ) << std::endl;
        abort();
    }

    //daughters for each V possibility
    switch ( Vtype ) {
        case VID::JPSI:
            if ( Id_daug1 != EvtPDL::getId( "mu+" ) ) {
                EvtGenReport( EVTGEN_ERROR, fname.c_str() )
                    << " Daughter1 is not a mu+, but a "
                    << EvtPDL::name( Id_daug1 ) << std::endl;
                abort();
            }
            if ( Id_daug2 != EvtPDL::getId( "mu-" ) ) {
                EvtGenReport( EVTGEN_ERROR, fname.c_str() )
                    << " Daughter2 is not a mu-, but a "
                    << EvtPDL::name( Id_daug2 ) << std::endl;
                abort();
            }
            EvtGenReport( EVTGEN_DEBUG, fname.c_str() )
                << " Decay mode successfully initialized : J/psi -> mu+ mu-"
                << std::endl;
            break;

        case VID::RHO:
        case VID::OMEGA:
        case VID::RHO_OMEGA_MIXING:
            if ( Id_daug1 != EvtPDL::getId( "pi+" ) ) {
                EvtGenReport( EVTGEN_ERROR, fname.c_str() )
                    << " Daughter1 is not a pi+, but a "
                    << EvtPDL::name( Id_daug1 ) << std::endl;
                abort();
            }
            if ( Id_daug2 != EvtPDL::getId( "pi-" ) ) {
                EvtGenReport( EVTGEN_ERROR, fname.c_str() )
                    << " Daughter2 is not a pi-, but a "
                    << EvtPDL::name( Id_daug2 ) << std::endl;
                abort();
            }
            if ( Vtype == VID::RHO )
                EvtGenReport( EVTGEN_DEBUG, fname.c_str() )
                    << " Decay mode successfully initialized : rho0 -> pi+ pi-"
                    << std::endl;
            if ( Vtype == VID::OMEGA )
                EvtGenReport( EVTGEN_DEBUG, fname.c_str() )
                    << " Decay mode successfully initialized : omega -> pi+ pi-"
                    << std::endl;
            if ( Vtype == VID::RHO_OMEGA_MIXING )
                EvtGenReport( EVTGEN_DEBUG, fname.c_str() )
                    << " Decay mode successfully initialized : rho-omega mixing -> pi+ pi-"
                    << std::endl;
            break;

        default:
            EvtGenReport( EVTGEN_ERROR, fname.c_str() )
                << "No decay mode chosen ! " << std::endl;
            abort();
            break;
    }

    //fix dynamics parameters
    switch ( Vtype ) {
        case VID::JPSI:
            A = 0.66;
            break;
        case VID::RHO:
        case VID::OMEGA:
        case VID::RHO_OMEGA_MIXING:
            A = 0.79;
            break;
        default:
            A = 0;
            break;
    }

    if ( verbose() ) {
        EvtGenReport( EVTGEN_DEBUG, fname.c_str() )
            << " V decay parameters : " << std::endl;
        EvtGenReport( EVTGEN_DEBUG, fname.c_str() )
            << "   - V density matrix rho00 A = " << A << std::endl;
    }
}

//------------------------------------------------------------------------
// Method 'decay'
//------------------------------------------------------------------------
void EvtV2VpVmForLambdaB2LambdaV::decay( EvtParticle* V )
{
    V->makeDaughters( getNDaug(), getDaugs() );

    //get Vp and Vm particles
    EvtParticle* Vp = V->getDaug( 0 );
    EvtParticle* Vm = V->getDaug( 1 );

    //get resonance masses
    // - V         -> mass given by EvtLambdaB2LambdaV class
    // - Vp & Vm   -> nominal mass
    double MASS_V = V->mass();
    double MASS_VM = 0;
    switch ( Vtype ) {
        case VID::JPSI:
            MASS_VM = EvtPDL::getMeanMass( EvtPDL::getId( "mu-" ) );
            break;
        case VID::RHO:
        case VID::OMEGA:
        case VID::RHO_OMEGA_MIXING:
            MASS_VM = EvtPDL::getMeanMass( EvtPDL::getId( "pi-" ) );
            break;
        default:
            MASS_VM = 0;
            break;
    }
    double MASS_VP = MASS_VM;

    //generate random angles
    double phi = EvtRandom::Flat( 0, 2 * EvtConst::pi );
    double theta = acos( EvtRandom::Flat( -1, +1 ) );

    // Get resonance 4-vectors
    double E_Vp = ( MASS_V * MASS_V + MASS_VP * MASS_VP - MASS_VM * MASS_VM ) /
                  ( 2 * MASS_V );
    double E_Vm = ( MASS_V * MASS_V + MASS_VM * MASS_VM - MASS_VP * MASS_VP ) /
                  ( 2 * MASS_V );
    double P = sqrt( E_Vp * E_Vp - Vp->mass() * Vp->mass() );

    EvtVector4R P_V = V->getP4();
    EvtParticle* Mother_V = V->getParent();
    EvtVector4R lambdab = Mother_V->getP4();

    const double E_V = ( P_V.get( 0 ) );
    const double px_M = lambdab.get( 1 );
    const double py_M = lambdab.get( 2 );
    const double pz_M = lambdab.get( 3 );
    const double E_M = lambdab.get( 0 );
    const double px_MSq = px_M * px_M;
    const double py_MSq = py_M * py_M;
    const double pT_MSq = px_MSq + py_MSq;
    const double pT_M = sqrt( pT_MSq );
    const double invPT_M = pT_M > 0.0 ? 1.0 / pT_M : 0.0;

    EvtVector4R q_lambdab2( E_M, pT_M, 0.0, pz_M );

    EvtVector4R q_lambdab3( E_M, q_lambdab2.get( 3 ), q_lambdab2.get( 1 ),
                            q_lambdab2.get( 2 ) );

    EvtVector4R q_V1( E_V,
                      invPT_M * ( px_M * P_V.get( 1 ) + py_M * P_V.get( 2 ) ),
                      invPT_M * ( -py_M * P_V.get( 1 ) + px_M * P_V.get( 2 ) ),
                      P_V.get( 3 ) );

    EvtVector4R q_V2( E_V, q_V1.get( 3 ), q_V1.get( 1 ), q_V1.get( 2 ) );

    const double px = -q_V2.get( 1 );
    const double py = -q_V2.get( 2 );
    const double pz = -q_V2.get( 3 );
    const double pxSq = px * px;
    const double pySq = py * py;
    const double pzSq = pz * pz;
    const double pTSq = pxSq + pySq;
    const double pSq = pTSq + pzSq;
    const double p = sqrt( pSq );
    const double invP = p > 0.0 ? 1.0 / p : 0.0;
    const double pT = sqrt( pTSq );
    const double invPT = pT > 0.0 ? 1.0 / pT : 0.0;

    EvtVector4R q_V4( q_V2.get( 0 ), 0.0, pT, p );

    EvtVector4R q_Vp1( E_Vp, P * sin( theta ) * cos( phi ),
                       P * sin( theta ) * sin( phi ), P * cos( theta ) );
    EvtVector4R q_Vm1( E_Vm, -P * sin( theta ) * cos( phi ),
                       -P * sin( theta ) * sin( phi ), -P * cos( theta ) );

    EvtVector4R q_Vp3( q_Vp1.get( 0 ),
                       invP * invPT *
                           ( q_Vp1.get( 1 ) * px * pz + q_Vp1.get( 2 ) * py * p -
                             q_Vp1.get( 3 ) * pT * px ),
                       invP * invPT *
                           ( q_Vp1.get( 1 ) * py * pz - q_Vp1.get( 2 ) * px * p -
                             q_Vp1.get( 3 ) * pT * py ),
                       -invP * ( q_Vp1.get( 1 ) * pT + q_Vp1.get( 3 ) * pz ) );

    EvtVector4R q_Vm3( q_Vm1.get( 0 ),
                       invP * invPT *
                           ( q_Vm1.get( 1 ) * px * pz + q_Vm1.get( 2 ) * py * p -
                             q_Vm1.get( 3 ) * pT * px ),
                       invP * invPT *
                           ( q_Vm1.get( 1 ) * py * pz - q_Vm1.get( 2 ) * px * p -
                             q_Vm1.get( 3 ) * pT * py ),
                       -invP * ( q_Vm1.get( 1 ) * pT + q_Vm1.get( 3 ) * pz ) );

    EvtVector4R q_Vp5( q_Vp3.get( 0 ), q_Vp3.get( 2 ), q_Vp3.get( 3 ),
                       q_Vp3.get( 1 ) );

    EvtVector4R q_Vm5( q_Vm3.get( 0 ), q_Vm3.get( 2 ), q_Vm3.get( 3 ),
                       q_Vm3.get( 1 ) );

    EvtVector4R q_Vp( q_Vp5.get( 0 ),
                      invPT_M * ( px_M * q_Vp5.get( 1 ) - py_M * q_Vp5.get( 2 ) ),
                      invPT_M * ( py_M * q_Vp5.get( 1 ) + px_M * q_Vp5.get( 2 ) ),
                      q_Vp5.get( 3 ) );

    EvtVector4R q_Vm( q_Vm5.get( 0 ),
                      invPT_M * ( px_M * q_Vm5.get( 1 ) - py_M * q_Vm5.get( 2 ) ),
                      invPT_M * ( py_M * q_Vm5.get( 1 ) + px_M * q_Vm5.get( 2 ) ),
                      q_Vm5.get( 3 ) );

    // Set particle 4-momenta
    Vp->init( getDaugs()[0], q_Vp );
    Vm->init( getDaugs()[1], q_Vm );

    //compute pdf
    double pdf = 0.0;
    if ( Vtype == VID::JPSI ) {
        //leptonic case
        pdf = ( 1.0 - 3.0 * A ) * cos( theta ) * cos( theta ) + ( 1.0 + A );
    } else {
        //hadronic case
        pdf = ( 3.0 * A - 1.0 ) * cos( theta ) * cos( theta ) + ( 1.0 - A );
    }

    //set probability
    setProb( pdf );

    if ( verbose() ) {
        EvtGenReport( EVTGEN_DEBUG, fname.c_str() )
            << " Angular angles  : theta = " << theta << " ; phi = " << phi
            << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Lambdab  px: " << px_M << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Lambdab  py: " << py_M << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Lambdab  pz: " << pz_M << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Lambdab  E: " << E_M << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Lambdab2  px:  " << q_lambdab2.get( 1 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Lambdab2  py:  " << q_lambdab2.get( 2 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Lambdab2  pz:  " << q_lambdab2.get( 3 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Lambdab2  E:   " << q_lambdab2.get( 0 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Lambdab3  px:  " << q_lambdab3.get( 1 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Lambdab3  py:  " << q_lambdab3.get( 2 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Lambdab3  pz:  " << q_lambdab3.get( 3 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Lambdab3  E:   " << q_lambdab3.get( 0 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " V 0  px:  " << P_V.get( 1 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " V 0  py:  " << P_V.get( 2 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " V 0  pz:  " << P_V.get( 3 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " V 0  E:   " << P_V.get( 0 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " V 1  px:  " << q_V1.get( 1 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " V 1  py:  " << q_V1.get( 2 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " V 1  pz:  " << q_V1.get( 3 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " V 1  E:   " << q_V1.get( 0 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " V 2  px:  " << q_V2.get( 1 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " V 2  py:  " << q_V2.get( 2 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " V 2  pz:  " << q_V2.get( 3 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " V 2  E:   " << q_V2.get( 0 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " V  px: " << px << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " V  py: " << py << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " V  pz: " << pz << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Vm 1 px:  " << q_Vm1.get( 1 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Vm 1 py:  " << q_Vm1.get( 2 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Vm 1 pz:  " << q_Vm1.get( 3 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Vm 1 E:   " << q_Vm1.get( 0 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Vm 3 px:  " << q_Vm3.get( 1 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Vm 3 px:  " << q_Vm3.get( 1 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Vm 3 py:  " << q_Vm3.get( 2 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Vm 3 pz:  " << q_Vm3.get( 3 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Vm 3 E:   " << q_Vm3.get( 0 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Vm 5 px:  " << q_Vm5.get( 1 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Vm 5 py:  " << q_Vm5.get( 2 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Vm 5 pz:  " << q_Vm5.get( 3 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Vm 5 E:   " << q_Vm5.get( 0 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Vp 1  px:  " << q_Vp1.get( 1 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Vp 1  py:  " << q_Vp1.get( 2 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Vp 1  pz:  " << q_Vp1.get( 3 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Vp 1  E:   " << q_Vp1.get( 0 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Vp 3  px:  " << q_Vp3.get( 1 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Vp 3  py:  " << q_Vp3.get( 2 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Vp 3  pz:  " << q_Vp3.get( 3 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Vp 3  E:   " << q_Vp3.get( 0 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Vp 5  px:  " << q_Vp5.get( 1 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Vp 5  py:  " << q_Vp5.get( 2 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Vp 5  pz:  " << q_Vp5.get( 3 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Vp 5  E:   " << q_Vp5.get( 0 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Vp  px:  " << q_Vp.get( 1 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Vp  py:  " << q_Vp.get( 2 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << "Vp  pz:  " << q_Vp.get( 3 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Vm px:  " << q_Vm.get( 1 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Vm py:  " << q_Vm.get( 2 ) << std::endl;
        EvtGenReport( EVTGEN_INFO, fname.c_str() )
            << " Vm pz:  " << q_Vm.get( 3 ) << std::endl;
        EvtGenReport( EVTGEN_DEBUG, fname.c_str() )
            << " V decay pdf value : " << pdf << std::endl;
    }
}
