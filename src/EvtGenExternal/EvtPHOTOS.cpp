
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

#ifdef EVTGEN_PHOTOS

#include "EvtGenExternal/EvtPHOTOS.hh"

#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtPhotonParticle.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtVector4R.hh"

#include <iostream>
#include <sstream>
#include <vector>

using std::endl;

// Mutex PHOTOS as it is not thread safe.
bool EvtPHOTOS::m_initialised = false;
std::mutex EvtPHOTOS::m_photos_mutex;

EvtPHOTOS::EvtPHOTOS( const std::string& photonType, const bool useEvtGenRandom,
                      const double infraredCutOff,
                      const double maxWtInterference ) :
    m_useEvtGenRandom{ useEvtGenRandom },
    m_infraredCutOff{ infraredCutOff },
    m_maxWtInterference{ maxWtInterference },
    m_photonType{ photonType }
{
}

void EvtPHOTOS::initialise()
{
    const std::lock_guard<std::mutex> lock( m_photos_mutex );

    m_gammaId = EvtPDL::getId( m_photonType );
    if ( m_gammaId == EvtId( -1, -1 ) ) {
        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << "Error in EvtPHOTOS. Do not recognise the photon type "
            << m_photonType << ". Setting this to \"gamma\". " << endl;
        m_gammaId = EvtPDL::getId( "gamma" );
    }
    m_gammaPDG = EvtPDL::getStdHep( m_gammaId );

    if ( m_initialised ) {
        return;
    }

    EvtGenReport( EVTGEN_INFO, "EvtGen" ) << "Setting up PHOTOS." << endl;

    if ( m_useEvtGenRandom ) {
        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << "Using EvtGen random number engine also for Photos++" << endl;

        Photospp::Photos::setRandomGenerator( EvtRandom::Flat );
    }

    Photospp::Photos::initialize();

    Photospp::Photos::setExponentiation( true );

    /* Sets the infrared cutoff (default at 1e-7)
    * Reset the minimum photon energy, if required, in units of half of the decaying particle mass.
    * This must be done after exponentiation! Keep the cut at 1e-7, i.e. 0.1 keV at the 1 GeV scale,
    * which is appropriate for B decays
    */
    Photospp::Photos::setInfraredCutOff( m_infraredCutOff );

    Photospp::Photos::setInterference( true );

    /* Increase the maximum possible value of the interference weight
     * corresponding to 2^n, where n = number of charges (+,-).
     */
    Photospp::Photos::maxWtInterference( m_maxWtInterference );

#ifdef EVTGEN_PHOTOS_NEWLIBS
    // This turns off/on virtual photons splitting into l+l-
    // It is false by default, but just to keep track of it here.
    Photospp::Photos::setPairEmission( false );
#endif

    m_initialised = true;
}

void EvtPHOTOS::doRadCorr( EvtParticle* theParticle )
{
    if ( !theParticle ) {
        return;
    }

    /* Skip running Photos if the particle has no daughters, since we can't add FSR.
     * Also skip Photos if the particle has too many daughters (>= 10) to avoid a problem
     * with a hard coded upper limit in the PHOENE subroutine.
     */
    const int nDaughters( theParticle->getNDaug() );
    if ( nDaughters == 0 || nDaughters >= 10 ) {
        return;
    }

    /* Create a dummy HepMC GenEvent containing a single vertex, with the mother
     * assigned as the incoming particle and its daughters as outgoing particles.
     * We then pass this event to Photos for processing.
     * It will return a modified version of the event, updating the momentum of
     * the original particles and will contain any new photon particles.
     * We add these extra photons to the mother particle daughter list.
     * First, Create the dummy event.
     */
    auto theEvent = std::make_unique<GenEvent>( Units::GEV, Units::MM );

    // Create the decay "vertex".
    GenVertexPtr theVertex = newGenVertexPtr();
    theEvent->add_vertex( theVertex );

    // Add the mother particle as the incoming particle to the vertex.
    const GenParticlePtr hepMCMother = this->createGenParticle( *theParticle,
                                                                true );
    theVertex->add_particle_in( hepMCMother );

    /* Find all daughter particles and assign them as outgoing particles to the vertex.
     * Keep track of the number of photons already in the decay (e.g. we may have B -> K* gamma)
     */
    int iDaug( 0 ), nGamma( 0 );
    for ( iDaug = 0; iDaug < nDaughters; iDaug++ ) {
        const EvtParticle& theDaughter = *theParticle->getDaug( iDaug );
        const GenParticlePtr hepMCDaughter =
            this->createGenParticle( theDaughter, false );
        theVertex->add_particle_out( hepMCDaughter );

        if ( theDaughter.getPDGId() == m_gammaPDG ) {
            nGamma++;
        }
    }

    {
        const std::lock_guard<std::mutex> lock( m_photos_mutex );
        /* Now pass the event to Photos for processing
     * Create a Photos event object */
#ifdef EVTGEN_HEPMC3
        Photospp::PhotosHepMC3Event photosEvent( theEvent.get() );
#else
        Photospp::PhotosHepMCEvent photosEvent( theEvent.get() );
#endif

        // Run the Photos algorithm
        photosEvent.process();
    }

    // Find the number of (outgoing) photons in the event
    const int nPhotons = this->getNumberOfPhotons( theVertex );

    int iLoop( 0 );
    // See whether Photos has created additional photons. If not, do nothing extra.
    if ( nPhotons > nGamma ) {
        // We have extra particles from Photos; these would have been appended
        // to the outgoing particle list

        // Get the iterator of outgoing particles for this vertex
#ifdef EVTGEN_HEPMC3
        for ( const auto& outParticle : theVertex->particles_out() ) {
#else
        HepMC::GenVertex::particles_out_const_iterator outIter;
        for ( outIter = theVertex->particles_out_const_begin();
              outIter != theVertex->particles_out_const_end(); ++outIter ) {
            // Get the next HepMC GenParticle
            const HepMC::GenParticle* outParticle = *outIter;
#endif
            // Get the three-momentum Photos result for this particle, and the PDG id
            if ( outParticle ) {
                const FourVector HepMCP4 = outParticle->momentum();
                const double px = HepMCP4.px();
                const double py = HepMCP4.py();
                const double pz = HepMCP4.pz();
                const double energy = HepMCP4.e();
                const int pdgId = outParticle->pdg_id();
                const EvtVector4R newP4( energy, px, py, pz );

                if ( iLoop < nDaughters ) {
                    /* Original daughters: keep the original particle mass, 
		     * but set the three-momentum according to what Photos 
		     * has modified. 
		     */
                    EvtParticle* daugParticle = theParticle->getDaug( iLoop );
                    if ( daugParticle ) {
                        daugParticle->setP4WithFSR( newP4 );
                    }

                } else if ( pdgId == m_gammaPDG ) {
                    // Extra photon particle. Setup the four-momentum object

                    // Create a new photon particle and add it to the list of daughters
                    EvtPhotonParticle* gamma = new EvtPhotonParticle();
                    gamma->init( m_gammaId, newP4 );
                    // Set the pre-FSR photon momentum to zero
                    gamma->setFSRP4toZero();
                    // Set its particle attribute to specify it is a FSR photon
                    gamma->setAttribute( "FSR", 1 );    // it is a FSR photon
                    // Let the mother know about this new photon
                    gamma->addDaug( theParticle );
                }

                // Increment the loop counter for detecting additional photon particles
                iLoop++;
            }
        }
    }

    // Cleanup
    theEvent->clear();

    return;
}

GenParticlePtr EvtPHOTOS::createGenParticle( const EvtParticle& theParticle,
                                             bool incoming ) const
{
    /* Method to create an HepMC::GenParticle version of the given EvtParticle.
     */

    // Get the 4-momentum (E, px, py, pz) for the EvtParticle
    const EvtVector4R p4 = incoming ? theParticle.getP4Restframe()
                                    : theParticle.getP4();

    // Convert this to the HepMC 4-momentum
    const double E = p4.get( 0 );
    const double px = p4.get( 1 );
    const double py = p4.get( 2 );
    const double pz = p4.get( 3 );

    const FourVector hepMC_p4( px, py, pz, E );

    const int PDGInt = EvtPDL::getStdHep( theParticle.getId() );

    // Set the status flag for the particle. This is required, otherwise Photos++
    // will crash from out-of-bounds array index problems.
    const int status = incoming ? Photospp::PhotosParticle::HISTORY
                                : Photospp::PhotosParticle::STABLE;

    const GenParticlePtr genParticle = newGenParticlePtr( hepMC_p4, PDGInt,
                                                          status );

    return genParticle;
}

int EvtPHOTOS::getNumberOfPhotons( const GenVertexPtr theVertex ) const
{
    /* Find the number of photons from the outgoing particle list
     */

    if ( !theVertex ) {
        return 0;
    }

    int nPhotons( 0 );

    // Get the iterator of outgoing particles for this vertex
#ifdef EVTGEN_HEPMC3
    for ( const auto& outParticle : theVertex->particles_out() ) {
#else
    HepMC::GenVertex::particles_out_const_iterator outIter;
    for ( outIter = theVertex->particles_out_const_begin();
          outIter != theVertex->particles_out_const_end(); ++outIter ) {
        // Get the next HepMC GenParticle
        const HepMC::GenParticle* outParticle = *outIter;
#endif

        // Get the PDG id
        int pdgId( 0 );
        if ( outParticle != nullptr ) {
            pdgId = outParticle->pdg_id();
        }

        // Keep track of how many photons there are
        if ( pdgId == m_gammaPDG ) {
            nPhotons++;
        }
    }

    return nPhotons;
}

#endif
