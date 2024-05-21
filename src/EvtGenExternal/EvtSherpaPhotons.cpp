
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

#ifdef EVTGEN_SHERPA

#include "EvtGenExternal/EvtSherpaPhotons.hh"

#include "EvtGenBase/EvtPhotonParticle.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtVector4R.hh"

#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Phys/Particle.H"
#include "SHERPA/Initialization/Initialization_Handler.H"
#include "SHERPA/Single_Events/Event_Handler.H"
#include "SHERPA/SoftPhysics/Soft_Photon_Handler.H"

#include <cstring>
#include <iostream>
#include <vector>

using std::endl;

EvtSherpaPhotons::EvtSherpaPhotons( const bool useEvtGenRandom,
                                    const double infraredCutOff, const int mode,
                                    const int useME ) :
    m_useEvtGenRandom{ useEvtGenRandom },
    m_infraredCutOff{ infraredCutOff },
    m_mode{ mode },
    m_useME{ useME }
{
}

void EvtSherpaPhotons::initialise()
{
    if ( m_initialised ) {
        return;
    }

    EvtGenReport( EVTGEN_INFO, "EvtGen" )
        << "Setting up Sherpa PHOTONS++ generator for FSR." << endl;

    /* Specify whether we are going to use EvtGen's random number generator 
     * (via EvtSherpaRandom interface) for Sherpa's PHOTONS++ simulation. */
    if ( m_useEvtGenRandom ) {
        m_configs.push_back( "SHERPA_LDADD=EvtGenExternal" );
        m_configs.push_back( "EXTERNAL_RNG=EvtSherpaRandom" );
    }

    /*
     * Internal PHOTONS++ settings.
     * Documentation at 
     * https://sherpa.hepforge.org/doc/Sherpa.html#QED-Corrections
     */

    /*
     * YFS_IR_CUTOFF sets the infrared cut-off which serves as a 
     * minimum photon energy in this frame for explicit photon generation.
     * The default is YFS_IR_CUTOFF = 1E-3 (GeV) but we set it to 1E-7 to equalize
     * with the cutoff used for PHOTOS.
     */
    m_configs.push_back( "YFS_IR_CUTOFF=" + ATOOLS::ToString( m_infraredCutOff ) );

    /*
     * The keyword YFS_MODE = [0,1,2] determines the mode of operation. 
     * YFS_MODE = 0 switches Photons off. 
     * YFS_MODE = 1 sets the mode to "soft only", meaning soft emissions will be 
     * treated correctly to all orders but no hard emission corrections will be included. 
     * YFS_MODE = 2 these hard emission corrections will also be included up to 
     * first order in alpha_QED. This is the default setting in Sherpa. 
     */

    m_configs.push_back( "YFS_MODE=" + ATOOLS::ToString( m_mode ) );

    /*
     * The switch YFS_USE_ME = [0,1] tells Photons how to correct hard emissions to 
     * first order in alpha_QED. If YFS_USE_ME = 0, then Photons will use collinearly 
     * approximated real emission matrix elements. Virtual emission matrix elements of 
     * order alpha_QED are ignored. If, however, YFS_USE_ME=1, then exact real and/or 
     * virtual emission matrix elements are used wherever possible. 
     * These are presently available for V->FF, V->SS, S->FF, S->SS, S->Slnu, 
     * S->Vlnu type decays, Z->FF decays and leptonic tau and W decays. For all other 
     * decay types general collinearly approximated matrix elements are used. 
     * In both approaches all hadrons are treated as point-like objects. 
     * The default setting is YFS_USE_ME = 1. This switch is only effective if YFS_MODE = 2.
     */
    m_configs.push_back( "YFS_USE_ME=" + ATOOLS::ToString( m_useME ) );

    // TODO: Virtual photon splitting into l+l- will be part of ME corrections from Sherpa 3.0.0 on,
    // Once released, this should be switched off by default or
    // taken into account while retrieving the event.

    // Set up run-time arguments for Sherpa.
    std::vector<char*> argv = this->addParameters();

    // Create instance and initialise Sherpa.
    m_sherpaGen = std::make_unique<SHERPA::Sherpa>();
    m_sherpaGen->InitializeTheRun( argv.size(), &argv[0] );
    m_sherpaGen->InitializeTheEventHandler();

    m_gammaId = EvtPDL::getId( m_photonType );
    m_gammaPDG = EvtPDL::getStdHep( m_gammaId );

    this->updateParticleLists();

    m_initialised = true;
}

std::vector<char*> EvtSherpaPhotons::addParameters()
{
    std::vector<char*> argv;
    argv.reserve( m_configs.size() );

    for ( auto& config : m_configs ) {
        if ( config != "Sherpa" ) {
            EvtGenReport( EVTGEN_INFO, "EvtGen" )
                << " EvtSherpaPhotons::initialise: Setting parameter '"
                << config << "'" << endl;
        }
        argv.push_back( config.data() );
    }

    return argv;
}

void EvtSherpaPhotons::updateParticleLists()
{
    /* Use the EvtGen decay and pdl tables (decay/user.dec and evt.pdl) 
     * to update Sherpa's KF_Table which contains all defined particles 
     * and their properties. 
     * Loop over all entries in the EvtPDL particle data table.
    */
    const int nPDL = EvtPDL::entries();

    for ( int iPDL = 0; iPDL < nPDL; iPDL++ ) {
        const EvtId particleId = EvtPDL::getEntry( iPDL );

        // Sherpa and EvtGen should use the same PDG codes for particles
        const int PDGCode = EvtPDL::getStdHep( particleId );
        const long unsigned int particleCode = std::abs( PDGCode );

        // If particle is found in Sherpa's KF_Table, then continue.
        // Updating the pole mass is not necessary because we provide the 4-momenta
        // and the final mass for final-state radiation generation.
        ATOOLS::KFCode_ParticleInfo_Map::const_iterator kf_it(
            ATOOLS::s_kftable.find( particleCode ) );
        if ( kf_it != ATOOLS::s_kftable.end() ) {
            continue;
        }

        // If the PDG is negative, use the charge conjugated ID because
        // the KF_Table is defined only for positive PDG numbers.
        const EvtId particleIdConj = EvtPDL::chargeConj( particleId );

        // Get mass, width, charge and spin
        const double mass = EvtPDL::getMeanMass( particleId );
        const double width = EvtPDL::getWidth( particleId );
        const int charge = PDGCode < 0 ? EvtPDL::chg3( particleIdConj )
                                       : EvtPDL::chg3( particleId );
        const int spin = EvtSpinType::getSpin2(
            EvtPDL::getSpinType( particleId ) );

        // Use as reference the name of the particle (not anti-particle).
        const std::string idName = PDGCode < 0 ? particleIdConj.getName()
                                               : particleId.getName();

        // Create new entry in the KF_Table. The 0 and 1 in the constructor below
        // stand for m_on = 0 and m_stable = 1. These properties are irrelevant
        // for final-state radiation and set to default values.
        ATOOLS::s_kftable[particleCode] = new ATOOLS::Particle_Info(
            particleCode, mass, width, charge, spin, 0, 1, idName, idName );
    }
}

void EvtSherpaPhotons::doRadCorr( EvtParticle* theParticle )
{
    if ( !theParticle ) {
        return;
    }

    /* Skip running FSR if the particle has no daughters, since we can't add FSR.
     * Also skip FSR if the particle has too many daughters (>= 10) as done for PHOTOS.
     */
    const int nDaughters( theParticle->getNDaug() );
    if ( nDaughters == 0 || nDaughters >= 10 ) {
        return;
    }

    ATOOLS::Blob_List* blobs = m_sherpaGen->GetEventHandler()->GetBlobs();

    // Use the hadron decay flag always.
    // Internally, Sherpa handles tau decays in the same way as hadron decays.
    ATOOLS::Blob* blob = blobs->AddBlob( ATOOLS::btp::Hadron_Decay );

    // Tell Sherpa that the blob needs FSR (that is extra QED)
    blob->SetStatus( ATOOLS::blob_status::needs_extraQED );

    // Create mother particle and add it to blob
    ATOOLS::Flavour mother_flav( EvtPDL::getStdHep( theParticle->getId() ) );
    mother_flav.SetStable( false );
    const double motherM0 = theParticle->mass();
    const ATOOLS::Vec4D mother_mom( motherM0, 0., 0., 0. );

    // Add mother to blob with m_number=-1. The m_number will start from 0 for daughters.
    ATOOLS::Particle* mother_part = new ATOOLS::Particle( -1, mother_flav,
                                                          mother_mom );
    mother_part->SetFinalMass( motherM0 );
    mother_part->SetStatus( ATOOLS::part_status::decayed );
    // Set m_info="I" (initial) for mother. The m_info will be "F" (final) for daughters.
    mother_part->SetInfo( 'I' );

    blob->AddToInParticles( mother_part );

    // Add daughters to the blob
    for ( int iDaug = 0; iDaug < nDaughters; iDaug++ ) {
        const EvtParticle* theDaughter = theParticle->getDaug( iDaug );
        ATOOLS::Flavour daughter_flav( EvtPDL::getStdHep( theDaughter->getId() ) );
        daughter_flav.SetStable( true );

        const double daugE = theDaughter->getP4().get( 0 );
        const double daugPx = theDaughter->getP4().get( 1 );
        const double daugPy = theDaughter->getP4().get( 2 );
        const double daugPz = theDaughter->getP4().get( 3 );
        const double daugM0 = theDaughter->mass();

        const ATOOLS::Vec4D daughter_mom( daugE, daugPx, daugPy, daugPz );

        ATOOLS::Particle* daughter_part =
            new ATOOLS::Particle( iDaug, daughter_flav, daughter_mom );
        daughter_part->SetFinalMass( daugM0 );

        daughter_part->SetStatus( ATOOLS::part_status::active );
        daughter_part->SetInfo( 'F' );

        blob->AddToOutParticles( daughter_part );
    }

    const SHERPA::Initialization_Handler* inithandler =
        m_sherpaGen->GetInitHandler();
    SHERPA::Soft_Photon_Handler* softphotonhandler =
        inithandler->GetSoftPhotonHandler();

    // Simulate the radiation
    softphotonhandler->AddRadiation( blob );

    // Get number of final-state particles
    const int nFinal( blob->NOutP() );

    // Retrieve the event from Sherpa blob if FSR photons were added
    if ( nFinal > nDaughters ) {
        for ( int iLoop = 0; iLoop < nFinal; iLoop++ ) {
            const ATOOLS::Particle* outParticle = blob->OutParticle( iLoop );
            const ATOOLS::Vec4D daughter_mom = outParticle->Momentum();
            const long int pdgId = outParticle->Flav();
            const EvtVector4R newP4( daughter_mom[0], daughter_mom[1],
                                     daughter_mom[2], daughter_mom[3] );
            const char daugInfo = outParticle->Info();

            if ( iLoop < nDaughters ) {
                // Update momenta of initial particles
                EvtParticle* daugParticle = theParticle->getDaug( iLoop );
                if ( daugParticle ) {
                    daugParticle->setP4WithFSR( newP4 );
                }
            } else if ( pdgId == m_gammaPDG && daugInfo == 'S' ) {
                // Add FSR photons. PHOTONS flags them with info 'S'.
                // Create a new photon particle and add it to the list of daughters
                EvtPhotonParticle* gamma = new EvtPhotonParticle();
                gamma->init( m_gammaId, newP4 );
                gamma->setFSRP4toZero();
                gamma->setAttribute( "FSR", 1 );
                gamma->addDaug( theParticle );
            }
        }
    }

    return;
}

#endif
