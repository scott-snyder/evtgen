
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

#include "EvtGenBase/EvtParticleDecayList.hh"

#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtStatus.hh"

#include <ctype.h>
#include <fstream>
#include <iostream>
#include <stdlib.h>
using std::endl;
using std::fstream;

EvtParticleDecayList::EvtParticleDecayList( const EvtParticleDecayList& o )
{
    m_nmode = o.m_nmode;
    m_rawbrfrsum = o.m_rawbrfrsum;
    m_decaylist = new EvtParticleDecayPtr[m_nmode];

    int i;
    for ( i = 0; i < m_nmode; i++ ) {
        m_decaylist[i] = new EvtParticleDecay;

        EvtDecayBase* tModel = o.m_decaylist[i]->getDecayModel();

        EvtDecayBase* tModelNew = tModel->clone();
        if ( tModel->getFSR() ) {
            tModelNew->setFSR();
        }
        if ( tModel->verbose() ) {
            tModelNew->setVerbose();
        }
        if ( tModel->summary() ) {
            tModelNew->setSummary();
        }
        std::vector<std::string> args;
        int j;
        for ( j = 0; j < tModel->getNArg(); j++ ) {
            args.push_back( tModel->getArgStr( j ) );
        }
        tModelNew->saveDecayInfo( tModel->getParentId(), tModel->getNDaug(),
                                  tModel->getDaugs(), tModel->getNArg(), args,
                                  tModel->getModelName(),
                                  tModel->getBranchingFraction() );
        m_decaylist[i]->setDecayModel( tModelNew );

        m_decaylist[i]->setBrfrSum( o.m_decaylist[i]->getBrfrSum() );
        m_decaylist[i]->setMassMin( o.m_decaylist[i]->getMassMin() );
    }
}

EvtParticleDecayList::~EvtParticleDecayList()
{
    int i;
    for ( i = 0; i < m_nmode; i++ ) {
        delete m_decaylist[i];
    }

    if ( m_decaylist != nullptr )
        delete[] m_decaylist;
}

void EvtParticleDecayList::printSummary()
{
    int i;
    for ( i = 0; i < m_nmode; i++ ) {
        m_decaylist[i]->printSummary();
    }
}

void EvtParticleDecayList::removeDecay()
{
    int i;
    for ( i = 0; i < m_nmode; i++ ) {
        delete m_decaylist[i];
    }

    delete[] m_decaylist;
    m_decaylist = nullptr;
    m_nmode = 0;
    m_rawbrfrsum = 0.0;
}

EvtDecayBase* EvtParticleDecayList::getDecayModel( int imode )
{
    EvtDecayBase* theModel( nullptr );
    if ( imode >= 0 && imode < m_nmode ) {
        EvtParticleDecay* theDecay = m_decaylist[imode];
        if ( theDecay != nullptr ) {
            theModel = theDecay->getDecayModel();
        }
    }

    return theModel;
}

EvtDecayBase* EvtParticleDecayList::getDecayModel( EvtParticle* p )
{
    if ( p->getNDaug() != 0 ) {
        assert( p->getChannel() >= 0 );
        return getDecay( p->getChannel() ).getDecayModel();
    }
    if ( p->getChannel() > ( -1 ) ) {
        return getDecay( p->getChannel() ).getDecayModel();
    }

    if ( getNMode() == 0 ) {
        return nullptr;
    }
    if ( getRawBrfrSum() < 0.00000001 ) {
        return nullptr;
    }

    if ( getNMode() == 1 ) {
        p->setChannel( 0 );
        return getDecay( 0 ).getDecayModel();
    }

    if ( p->getChannel() > ( -1 ) ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" ) << "Internal error!!!" << endl;
        ::abort();
    }

    int j;

    for ( j = 0; j < 10000000; j++ ) {
        double u = EvtRandom::Flat();

        int i;
        bool breakL = false;
        for ( i = 0; i < getNMode(); i++ ) {
            if ( breakL )
                continue;
            if ( u < getDecay( i ).getBrfrSum() ) {
                breakL = true;
                //special case for decay of on particel to another
                // e.g. K0->K0S

                if ( getDecay( i ).getDecayModel()->getNDaug() == 1 ) {
                    p->setChannel( i );
                    return getDecay( i ).getDecayModel();
                }

                if ( p->hasValidP4() ) {
                    if ( getDecay( i ).getMassMin() < p->mass() ) {
                        p->setChannel( i );
                        return getDecay( i ).getDecayModel();
                    }
                } else {
                    //Lange apr29-2002 - dont know the mass yet
                    p->setChannel( i );
                    return getDecay( i ).getDecayModel();
                }
            }
        }
    }

    //Ok, we tried 10000000 times above to pick a decay channel that is
    //kinematically allowed! Now we give up and search all channels!
    //if that fails, the particle will not be decayed!

    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "Tried 10000000 times to generate decay of "
        << EvtPDL::name( p->getId() ) << " with mass=" << p->mass() << endl;
    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "Will take first kinematically allowed decay in the decay table"
        << endl;

    int i;

    //Need to check that we don't use modes with 0 branching fractions.
    double previousBrSum = 0.0;
    for ( i = 0; i < getNMode(); i++ ) {
        if ( getDecay( i ).getBrfrSum() != previousBrSum ) {
            if ( getDecay( i ).getMassMin() < p->mass() ) {
                p->setChannel( i );
                return getDecay( i ).getDecayModel();
            }
        }
        previousBrSum = getDecay( i ).getBrfrSum();
    }

    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "Could not decay:" << EvtPDL::name( p->getId() ).c_str()
        << " with mass:" << p->mass() << " will throw event away! " << endl;

    EvtStatus::setRejectFlag();
    return nullptr;
}

void EvtParticleDecayList::setNMode( int nmode )
{
    EvtParticleDecayPtr* m_decaylist_new = new EvtParticleDecayPtr[nmode];

    if ( m_nmode != 0 ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Error m_nmode not equal to zero!!!" << endl;
        ::abort();
    }
    if ( m_decaylist != nullptr ) {
        delete[] m_decaylist;
    }
    m_decaylist = m_decaylist_new;
    m_nmode = nmode;
}

EvtParticleDecay& EvtParticleDecayList::getDecay( int nchannel ) const
{
    if ( nchannel >= m_nmode ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Error getting channel:" << nchannel << " with only " << m_nmode
            << " stored!" << endl;
        ::abort();
    }
    return *( m_decaylist[nchannel] );
}

void EvtParticleDecayList::makeChargeConj( EvtParticleDecayList* conjDecayList )
{
    m_rawbrfrsum = conjDecayList->m_rawbrfrsum;

    setNMode( conjDecayList->m_nmode );

    int i;

    for ( i = 0; i < m_nmode; i++ ) {
        m_decaylist[i] = new EvtParticleDecay;
        m_decaylist[i]->chargeConj( conjDecayList->m_decaylist[i] );
    }
}

void EvtParticleDecayList::addMode( EvtDecayBase* decay, double brfrsum,
                                    double massmin )
{
    EvtParticleDecayPtr* newlist = new EvtParticleDecayPtr[m_nmode + 1];

    int i;
    for ( i = 0; i < m_nmode; i++ ) {
        newlist[i] = m_decaylist[i];
    }

    m_rawbrfrsum = brfrsum;

    newlist[m_nmode] = new EvtParticleDecay;

    newlist[m_nmode]->setDecayModel( decay );
    newlist[m_nmode]->setBrfrSum( brfrsum );
    newlist[m_nmode]->setMassMin( massmin );

    EvtDecayBase* newDec = newlist[m_nmode]->getDecayModel();
    for ( i = 0; i < m_nmode; i++ ) {
        if ( newDec->matchingDecay( *( newlist[i]->getDecayModel() ) ) ) {
            //sometimes its ok..
            if ( newDec->getModelName() == "JETSET" ||
                 newDec->getModelName() == "PYTHIA" )
                continue;
            if ( newDec->getModelName() == "JSCONT" ||
                 newDec->getModelName() == "PYCONT" )
                continue;
            if ( newDec->getModelName() == "PYGAGA" )
                continue;
            if ( newDec->getModelName() == "LUNDAREALAW" )
                continue;
            if ( newDec->getModelName() == "TAUOLA" )
                continue;
            EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                << "Two matching decays with same parent in decay table\n";
            EvtGenReport( EVTGEN_ERROR, "EvtGen" ) << "Please fix that\n";
            EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                << "Parent " << EvtPDL::name( newDec->getParentId() ).c_str()
                << endl;
            for ( int j = 0; j < newDec->getNDaug(); j++ )
                EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                    << "Daughter "
                    << EvtPDL::name( newDec->getDaug( j ) ).c_str() << endl;
            assert( 0 );
        }
    }

    if ( m_nmode != 0 ) {
        delete[] m_decaylist;
    }

    if ( ( m_nmode == 0 ) && ( m_decaylist != nullptr ) )
        delete[] m_decaylist;

    m_nmode++;

    m_decaylist = newlist;
}

void EvtParticleDecayList::finalize()
{
    if ( m_nmode > 0 ) {
        if ( m_rawbrfrsum < 0.000001 ) {
            EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                << "Please give me a "
                << "branching fraction sum greater than 0\n";
            assert( 0 );
        }
        if ( fabs( m_rawbrfrsum - 1.0 ) > 0.0001 ) {
            EvtGenReport( EVTGEN_INFO, "EvtGen" )
                << "Warning, sum of branching fractions for "
                << EvtPDL::name( m_decaylist[0]->getDecayModel()->getParentId() )
                       .c_str()
                << " is " << m_rawbrfrsum << endl;
            EvtGenReport( EVTGEN_INFO, "EvtGen" ) << "rescaled to one! " << endl;
        }

        int i;

        for ( i = 0; i < m_nmode; i++ ) {
            double brfrsum = m_decaylist[i]->getBrfrSum() / m_rawbrfrsum;
            m_decaylist[i]->setBrfrSum( brfrsum );
        }
    }
}

EvtParticleDecayList& EvtParticleDecayList::operator=( const EvtParticleDecayList& o )
{
    if ( this != &o ) {
        removeDecay();
        m_nmode = o.m_nmode;
        m_rawbrfrsum = o.m_rawbrfrsum;
        m_decaylist = new EvtParticleDecayPtr[m_nmode];

        int i;
        for ( i = 0; i < m_nmode; i++ ) {
            m_decaylist[i] = new EvtParticleDecay;

            EvtDecayBase* tModel = o.m_decaylist[i]->getDecayModel();

            EvtDecayBase* tModelNew = tModel->clone();
            if ( tModel->getFSR() ) {
                tModelNew->setFSR();
            }
            if ( tModel->verbose() ) {
                tModelNew->setVerbose();
            }
            if ( tModel->summary() ) {
                tModelNew->setSummary();
            }
            std::vector<std::string> args;
            int j;
            for ( j = 0; j < tModel->getNArg(); j++ ) {
                args.push_back( tModel->getArgStr( j ) );
            }
            tModelNew->saveDecayInfo( tModel->getParentId(), tModel->getNDaug(),
                                      tModel->getDaugs(), tModel->getNArg(),
                                      args, tModel->getModelName(),
                                      tModel->getBranchingFraction() );
            m_decaylist[i]->setDecayModel( tModelNew );

            //m_decaylist[i]->setDecayModel(tModel);
            m_decaylist[i]->setBrfrSum( o.m_decaylist[i]->getBrfrSum() );
            m_decaylist[i]->setMassMin( o.m_decaylist[i]->getMassMin() );
        }
    }
    return *this;
}

void EvtParticleDecayList::removeMode( EvtDecayBase* decay )
{
    // here we will delete a decay with the same final state particles
    // and recalculate the branching fractions for the remaining modes
    int match = -1;
    int i;
    double match_bf;

    for ( i = 0; i < m_nmode; i++ ) {
        if ( decay->matchingDecay( *( m_decaylist[i]->getDecayModel() ) ) ) {
            match = i;
        }
    }

    if ( match < 0 ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << " Attempt to remove undefined mode for" << endl
            << "Parent " << EvtPDL::name( decay->getParentId() ).c_str() << endl
            << "Daughters: ";
        for ( int j = 0; j < decay->getNDaug(); j++ )
            EvtGenReport( EVTGEN_ERROR, "" )
                << EvtPDL::name( decay->getDaug( j ) ).c_str() << " ";
        EvtGenReport( EVTGEN_ERROR, "" ) << endl;
        ::abort();
    }

    if ( match == 0 ) {
        match_bf = m_decaylist[match]->getBrfrSum();
    } else {
        match_bf = ( m_decaylist[match]->getBrfrSum() -
                     m_decaylist[match - 1]->getBrfrSum() );
    }

    double divisor = 1 - match_bf;
    if ( divisor < 0.000001 && m_nmode > 1 ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Removing requested mode leaves "
            << EvtPDL::name( decay->getParentId() ).c_str()
            << " with zero sum branching fraction," << endl
            << "but more than one decay mode remains. Aborting." << endl;
        ::abort();
    }

    EvtParticleDecayPtr* newlist = new EvtParticleDecayPtr[m_nmode - 1];

    for ( i = 0; i < match; i++ ) {
        newlist[i] = m_decaylist[i];
        newlist[i]->setBrfrSum( newlist[i]->getBrfrSum() / divisor );
    }
    for ( i = match + 1; i < m_nmode; i++ ) {
        newlist[i - 1] = m_decaylist[i];
        newlist[i - 1]->setBrfrSum(
            ( newlist[i - 1]->getBrfrSum() - match_bf ) / divisor );
    }

    delete[] m_decaylist;

    m_nmode--;

    m_decaylist = newlist;

    if ( m_nmode == 0 ) {
        delete[] m_decaylist;
    }
}

bool EvtParticleDecayList::isJetSet() const
{
    int i;
    EvtDecayBase* decayer;

    for ( i = 0; i < getNMode(); i++ ) {
        decayer = getDecay( i ).getDecayModel();
        if ( decayer->getModelName() == "PYTHIA" )
            return true;
    }

    return false;
}
