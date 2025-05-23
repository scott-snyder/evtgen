
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

#ifdef EVTGEN_VINCIA
#include "EvtGenExternal/EvtVinciaQED.hh"

#include "EvtGenBase/EvtDecayTable.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticleFactory.hh"
#include "EvtGenBase/EvtPhotonParticle.hh"
#include "EvtGenBase/EvtReport.hh"

#include "Pythia8/Event.h"
#include "Pythia8/ParticleData.h"

using std::endl;

EvtVinciaQED::EvtVinciaQED( std::string xmlDir, bool useEvtGenRandom,
                            double infraredCutOff ) :
    m_useEvtGenRandom{ useEvtGenRandom }, m_infraredCutOff{ infraredCutOff }
{
    /*
     * Create a Pythia generator to call the Vincia shower model.
     */

    EvtGenReport( EVTGEN_INFO, "EvtGen" )
        << " Setting up Pythia generator to call Vincia generator for FSR."
        << endl;
    m_vinciaPythiaGen = std::make_unique<Pythia8::Pythia>( xmlDir );

    if ( m_useEvtGenRandom ) {
        m_evtgenRandom = std::make_shared<EvtPythiaRandom>();
    }
}

void EvtVinciaQED::initialise()
{
    if ( m_initialised ) {
        return;
    }

    /* Set the random number generator
     * Specify if we are going to use the random number generator (engine)
     * from EvtGen for Pythia 8.
     */
    if ( m_useEvtGenRandom ) {
        m_vinciaPythiaGen->setRndmEnginePtr( m_evtgenRandom );
    }

    /* Hadron-level processes only (hadronized, string fragmentation and secondary decays).
     * We do not want to generate the full pp or e+e- event structure etc..
     */
    m_vinciaPythiaGen->readString( "ProcessLevel:all = off" );

    // Turn off Pythia warnings, e.g. changes to particle properties
    m_vinciaPythiaGen->readString( "Print:quiet = on" );

    // Setup to use Vincia QED shower
    m_vinciaPythiaGen->readString( "PartonShowers:model = 2" );

    // This turns off virtual photons splitting into l+l-
    m_vinciaPythiaGen->readString( "Vincia:nGammaToLepton = 0" );

    /* Sets the infrared cutoff (default at 1e-7)
     * Note that the Vincia default is 1e-6, so we are setting it lower
     * to equalize with PHOTOS configuration.
     */
    std::stringstream QminChgL;
    QminChgL << "Vincia:QminChgL = " << m_infraredCutOff;
    m_vinciaPythiaGen->readString( QminChgL.str() );

    m_gammaId = EvtPDL::getId( "gamma" );
    m_gammaPDG = EvtPDL::getStdHep( m_gammaId );

    m_vinciaPythiaGen->init();

    m_initialised = true;
}

void EvtVinciaQED::doRadCorr( EvtParticle* theParticle )
{
    if ( !theParticle ) {
        return;
    }

    if ( !m_initialised ) {
        this->initialise();
    }

    const EvtId particleId = theParticle->getId();

    const Pythia8::TimeShowerPtr timeShowerPtr =
        m_vinciaPythiaGen->getShowerModelPtr()->getTimeShower();

    /* Need to use the reference to the Pythia8::Event object,
     * otherwise it will just return a new empty, default event object.
     */
    Pythia8::Event& theEvent = m_vinciaPythiaGen->event;
    theEvent.reset();

    /* Skip running FSR if the particle has no daughters, since we can't add FSR.
     * Also skip FSR if the particle has too many daughters (>= 10) as done for PHOTOS.
     */
    const int nDaughters( theParticle->getNDaug() );
    if ( nDaughters == 0 || nDaughters >= 10 ) {
        return;
    }

    const int PDGCode = EvtPDL::getStdHep( particleId );
    // Status = -1 Means that the particle has decayed and is not remaining
    const int status( -1 );
    const int grandMother1( 0 );
    const int grandMother2( 0 );
    const int daughter1( 2 );
    const int daughter2( nDaughters + 1 );

    const int colour( 0 ), anticolour( 0 );
    const double px( 0.0 ), py( 0.0 ), pz( 0.0 );
    const double m0 = theParticle->mass();
    const double E = m0;

    theEvent.append( PDGCode, status, grandMother1, grandMother2, daughter1,
                     daughter2, colour, anticolour, px, py, pz, E, m0 );

    // Status = 91 means normal decay product
    const int daugStatus( 91 );

    // Add all the daughters and keep track of number of photons
    for ( int iDaug = 0; iDaug < nDaughters; iDaug++ ) {
        const EvtParticle& theDaughter = *theParticle->getDaug( iDaug );
        const int daugPDGCode = theDaughter.getPDGId();

        // Particle has only one mother in EvtGen, so first mother points to parent Id and second to none.
        const int mother1( 1 ), mother2( 0 ), grandDaughter1( 0 ),
            grandDaughter2( 0 );

        const int daugColour( 0 ), daugAntiColour( 0 );
        const double daugPx = theDaughter.getP4().get( 1 );
        const double daugPy = theDaughter.getP4().get( 2 );
        const double daugPz = theDaughter.getP4().get( 3 );
        const double daugE = theDaughter.getP4().get( 0 );
        const double daugM0 = theDaughter.mass();

        theEvent.append( daugPDGCode, daugStatus, mother1, mother2,
                         grandDaughter1, grandDaughter2, daugColour,
                         daugAntiColour, daugPx, daugPy, daugPz, daugE, daugM0 );
    }

    /* Generate the Vincia QED shower
     * It updates the event internally
     */
    timeShowerPtr->showerQED( daughter1, daughter2, theEvent, m0 );

    // Get number of final-state particles
    const int nFinal( theEvent.nFinal() );

    // Retrieve the event from Pythia if FSR photons were added
    if ( nFinal > nDaughters ) {
        /* Store the daughters for this particle from the Pythia decay tree.
         * This is a recursive function that will continue looping through
         * all available daughters until the first set of final-state
         * particles are encountered in the Pythia Event structure.
         * It will add FSR photons to and update the daughters
         * of the (parent) EvtGen particle.
         */
        this->storeDaughterInfo( theEvent, theParticle, 1, 1 );
    }
}

void EvtVinciaQED::storeDaughterInfo( Pythia8::Event& theEvent,
                                      EvtParticle* theParticle, int startInt,
                                      int ancestorId )
{
    const std::vector<int> daugList = theEvent.daughterList( startInt );

    std::vector<int>::const_iterator daugIter;
    for ( daugIter = daugList.begin(); daugIter != daugList.end(); ++daugIter ) {
        const int daugInt = *daugIter;

        Pythia8::Particle& daughter_part = theEvent[daugInt];

        const int statusCode = daughter_part.status();

        /* Keep track of the initial final-state daughters
         * We need to do so because we have to update the momenta of the EvtGen particles.
         * In the Vincia event record, the status of the initial daughters
         * is 91 if they did not radiate, or -91 if they did.
         * Intermediate or final daughters in the shower record have different status.
         * When we loop recursively through the shower record, the (absolute) status of the particles
         * identifies whether they are the initial daughters.
         * If the current particle is not an initial daughter,
         * we do not use the current index, but the provided ancestor index.
         */
        const int daugId = std::abs( statusCode ) == 91 ? daugInt : ancestorId;

        if ( !daughter_part.isFinal() ) {
            // Recursively search for correct daughter type
            this->storeDaughterInfo( theEvent, theParticle, daugInt, daugId );

        } else {
            /* We have a daughter that is a final-state particle.
             * Make sure we are not double counting particles.
             * Set the status flag for any "new" particle to say that we have stored it.
             * Use status flag = 1000 (within the user allowed range for Pythia codes).
             */

            // Check that the status flag for the particle is not equal to 1000
            if ( statusCode == 1000 ) {
                continue;
            }

            const int daugPDGInt = daughter_part.id();

            const double px = daughter_part.px();
            const double py = daughter_part.py();
            const double pz = daughter_part.pz();
            const double E = daughter_part.e();
            const EvtVector4R newP4( E, px, py, pz );

            /* Status = 51 means that this is an outgoing particle produced by final-state shower
             * Actually, 51 means exactly outgoing produced by parton branching,
             * but it is the status used for particles produced by the Vincia QED shower.
             */
            const bool daugIsFSRphoton =
                ( daugPDGInt == m_gammaPDG && statusCode == 51 ) ? true : false;

            /* Id of ancestor to correct momenta of existing particles
             * Remove offset of 2 due to event header and parent in Pythia
             */
            const int daugIndex = daugId - 2;

            /* Set the status flag for the Pythia particle to let us know
             * that we have already considered it to avoid double counting.
             */
            daughter_part.status( 1000 );

            // Update daughters of EvtGen particle and attach new FSR photons
            if ( !daugIsFSRphoton ) {
                // Update momenta of initial particles
                EvtParticle* daugParticle = theParticle->getDaug( daugIndex );
                if ( daugParticle ) {
                    daugParticle->setP4WithFSR( newP4 );
                }
            } else {
                // Create a new photon particle and add it to the list of daughters
                EvtPhotonParticle* gamma = new EvtPhotonParticle();
                gamma->init( m_gammaId, newP4 );
                gamma->setFSRP4toZero();
                gamma->setAttribute( "FSR", 1 );
                gamma->addDaug( theParticle );
            }
        }
    }
}

#endif
