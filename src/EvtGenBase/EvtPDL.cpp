
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

#include "EvtGenBase/EvtPDL.hh"

#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtPartProp.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtReport.hh"

#include <cstring>
#include <fstream>
#include <iostream>

using std::endl;

EvtPDL& EvtPDL::getInstance()
{
    static thread_local EvtPDL theInstance;
    return theInstance;
}

void EvtPDL::reset()
{
    m_firstAlias = std::numeric_limits<std::size_t>::max();
    m_partlist.clear();
    m_particleNameLookup.clear();
}

void EvtPDL::read( const std::string& fname )
{
    std::ifstream pdtIn( fname );
    if ( !pdtIn ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Could not open:" << fname << "EvtPDL" << endl;
        return;
    }
    readPDT( pdtIn );
    pdtIn.close();
}

void EvtPDL::readPDT( std::istream& indec )
{
    getInstance().reset();

    char cmnd[100];
    char xxxx[100];

    char pname[100];
    int stdhepid;
    double mass;
    double pwidth;
    double pmaxwidth;
    int chg3;
    int spin2;
    double ctau;
    int lundkc;
    EvtId i;
    std::size_t nentries = 0;

    indec.seekg( 0, std::ios::beg );

    do {
        char ch, ch1;

        do {
            indec.get( ch );
            if ( ch == '\n' )
                indec.get( ch );
            if ( ch != '*' ) {
                indec.putback( ch );
            } else {
                while ( indec.get( ch1 ), ch1 != '\n' )
                    ;
            }
        } while ( ch == '*' );

        indec >> cmnd;

        if ( strcmp( cmnd, "end" ) ) {
            if ( !strcmp( cmnd, "add" ) ) {
                indec >> xxxx;
                indec >> xxxx;
                indec >> pname;
                indec >> stdhepid;
                indec >> mass;
                indec >> pwidth;
                indec >> pmaxwidth;
                indec >> chg3;
                indec >> spin2;
                indec >> ctau;
                indec >> lundkc;

                i = EvtId( nentries, nentries );

                EvtPartProp tmp;

                tmp.setSpinType( EvtSpinType::SCALAR );

                if ( spin2 == 0 )
                    tmp.setSpinType( EvtSpinType::SCALAR );
                if ( spin2 == 1 )
                    tmp.setSpinType( EvtSpinType::DIRAC );
                if ( spin2 == 2 )
                    tmp.setSpinType( EvtSpinType::VECTOR );
                if ( spin2 == 3 )
                    tmp.setSpinType( EvtSpinType::RARITASCHWINGER );
                if ( spin2 == 4 )
                    tmp.setSpinType( EvtSpinType::TENSOR );
                if ( spin2 == 5 )
                    tmp.setSpinType( EvtSpinType::SPIN5HALF );
                if ( spin2 == 6 )
                    tmp.setSpinType( EvtSpinType::SPIN3 );
                if ( spin2 == 7 )
                    tmp.setSpinType( EvtSpinType::SPIN7HALF );
                if ( spin2 == 8 )
                    tmp.setSpinType( EvtSpinType::SPIN4 );
                if ( spin2 == 2 && mass < 0.0001 )
                    tmp.setSpinType( EvtSpinType::PHOTON );
                if ( spin2 == 1 && mass < 0.0001 )
                    tmp.setSpinType( EvtSpinType::NEUTRINO );

                if ( !strcmp( pname, "string" ) ) {
                    tmp.setSpinType( EvtSpinType::STRING );
                }

                if ( !strcmp( pname, "vpho" ) ) {
                    tmp.setSpinType( EvtSpinType::VECTOR );
                }

                tmp.setId( i );
                tmp.setIdChgConj( EvtId( -1, -1 ) );
                tmp.setStdHep( stdhepid );
                tmp.setLundKC( lundkc );
                tmp.setName( pname );
                if ( getInstance().m_particleNameLookup.find( std::string(
                         pname ) ) != getInstance().m_particleNameLookup.end() ) {
                    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                        << "The particle name:" << pname
                        << " is already defined." << endl;
                    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                        << "Will terminate execution.";
                    ::abort();
                }
                getInstance().m_particleNameLookup[std::string( pname )] = nentries;
                tmp.setctau( ctau );
                tmp.setChg3( chg3 );

                tmp.initLineShape( mass, pwidth, pmaxwidth );

                getInstance().m_partlist.push_back( tmp );
                nentries++;
            }

            // if find a set read information and discard it

            if ( !strcmp( cmnd, "set" ) ) {
                indec >> xxxx;
                indec >> xxxx;
                indec >> xxxx;
                indec >> xxxx;
            }
        }

    } while ( strcmp( cmnd, "end" ) );
}

void EvtPDL::aliasChgConj( EvtId a, EvtId abar )
{
    if ( EvtPDL::chargeConj( EvtId( a.getId(), a.getId() ) ) !=
         EvtId( abar.getId(), abar.getId() ) ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Can't charge conjugate the two aliases:"
            << EvtPDL::name( a ).c_str() << " and "
            << EvtPDL::name( abar ).c_str() << endl;

        ::abort();
    }

    getInstance().m_partlist[a.getAlias()].setIdChgConj( abar );
    getInstance().m_partlist[abar.getAlias()].setIdChgConj( a );
}

EvtId EvtPDL::chargeConj( EvtId id )
{
    int index = id.getAlias();
    EvtId idchg;
    if ( index > -1 ) {
        idchg = getInstance().m_partlist[id.getAlias()].getIdChgConj();
    }

    if ( idchg != EvtId( -1, -1 ) )
        return idchg;

    if ( id.getId() != id.getAlias() ) {
        if ( chargeConj( EvtId( id.getId(), id.getId() ) ) ==
             EvtId( id.getId(), id.getId() ) ) {
            getInstance().m_partlist[id.getAlias()].setIdChgConj( id );
            return id;
        }
    }

    if ( id.getAlias() != id.getId() ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Trying to charge conjugate alias particle:" << name( id ).c_str()
            << " without defining the alias!" << endl;

        ::abort();
    }

    for ( size_t i = 0; i < getInstance().m_partlist.size(); i++ ) {
        if ( getInstance().m_partlist[i].getStdHep() ==
             -getInstance().m_partlist[id.getId()].getStdHep() ) {
            getInstance().m_partlist[id.getId()].setIdChgConj(
                getInstance().m_partlist[i].getId() );
            return getInstance().m_partlist[i].getId();
        }
    }

    getInstance().m_partlist[id.getId()].setIdChgConj( id );
    return id;
}

EvtId EvtPDL::evtIdFromStdHep( int stdhep )
{
    for ( size_t i = 0; i < getInstance().m_partlist.size(); i++ ) {
        if ( getInstance().m_partlist[i].getStdHep() == stdhep )
            return getInstance().m_partlist[i].getId();
    }

    return EvtId( -1, -1 );
}

void EvtPDL::alias( EvtId num, const std::string& newname )
{
    if ( getInstance().m_firstAlias < getInstance().m_partlist.size() ) {
        for ( size_t i = getInstance().m_firstAlias;
              i < getInstance().m_partlist.size(); i-- ) {
            if ( newname == getInstance().m_partlist[i].getName() ) {
                EvtGenReport( EVTGEN_WARNING, "EvtGen" )
                    << "Redefining alias:" << newname.c_str()
                    << " will be ignored!" << endl;
                return;
            }
        }
    } else {
        getInstance().m_firstAlias = getInstance().m_partlist.size();
    }

    getInstance().m_partlist.push_back( getInstance().m_partlist[num.getId()] );
    int entry = getInstance().m_partlist.size() - 1;
    getInstance().m_partlist[entry].setName( newname );
    if ( getInstance().m_particleNameLookup.find( std::string( newname ) ) !=
         getInstance().m_particleNameLookup.end() ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "The particle name:" << newname << " is already defined." << endl;
        EvtGenReport( EVTGEN_ERROR, "EvtGen" ) << "Will terminate execution.";
        ::abort();
    }
    getInstance().m_particleNameLookup[std::string( newname )] = entry;
    getInstance().m_partlist[entry].setId( EvtId( num.getId(), entry ) );
    //Lange - Dec7, 2003. Unset the charge conjugate.
    getInstance().m_partlist[entry].setIdChgConj( EvtId( -1, -1 ) );
}

EvtId EvtPDL::getId( const std::string& name )
{
    std::map<std::string, int>::iterator it =
        getInstance().m_particleNameLookup.find( std::string( name ) );
    if ( it == getInstance().m_particleNameLookup.end() )
        return EvtId( -1, -1 );

    return getInstance().m_partlist[it->second].getId();
}

// Function to get EvtId from LundKC ( == Pythia Hep Code , KF )
EvtId EvtPDL::evtIdFromLundKC( int pythiaId )
{
    unsigned int i;

    for ( i = 0; i < getInstance().m_partlist.size(); i++ ) {
        if ( getInstance().m_partlist[i].getLundKC() == pythiaId )
            return getInstance().m_partlist[i].getId();
    }

    return EvtId( -1, -1 );
}

double EvtPDL::getMeanMass( EvtId i )
{
    return getInstance().m_partlist[i.getId()].getMass();
}

double EvtPDL::getMass( EvtId i )
{
    return getInstance().m_partlist[i.getId()].rollMass();
}

double EvtPDL::getRandMass( EvtId i, EvtId* parId, int nDaug, EvtId* dauId,
                            EvtId* othDaugId, double maxMass, double* dauMasses )
{
    return getInstance().m_partlist[i.getId()].getRandMass( parId, nDaug, dauId,
                                                            othDaugId, maxMass,
                                                            dauMasses );
}

double EvtPDL::getMassProb( EvtId i, double mass, double massPar, int nDaug,
                            double* massDau )
{
    return getInstance().m_partlist[i.getId()].getMassProb( mass, massPar,
                                                            nDaug, massDau );
}

double EvtPDL::getMaxMass( EvtId i )
{
    return getInstance().m_partlist[i.getId()].getMassMax();
}

double EvtPDL::getMinMass( EvtId i )
{
    return getInstance().m_partlist[i.getId()].getMassMin();
}

double EvtPDL::getMaxRange( EvtId i )
{
    return getInstance().m_partlist[i.getId()].getMaxRange();
}

double EvtPDL::getWidth( EvtId i )
{
    return getInstance().m_partlist[i.getId()].getWidth();
}

double EvtPDL::getctau( const EvtId i )
{
    return getInstance().m_partlist[i.getId()].getctau();
}

int EvtPDL::getStdHep( EvtId id )
{
    return getInstance().m_partlist[id.getId()].getStdHep();
}

int EvtPDL::getLundKC( EvtId id )
{
    return getInstance().m_partlist[id.getId()].getLundKC();
}

int EvtPDL::chg3( EvtId i )
{
    return getInstance().m_partlist[i.getId()].getChg3();
}

EvtSpinType::spintype EvtPDL::getSpinType( EvtId i )
{
    return getInstance().m_partlist[i.getId()].getSpinType();
}

std::string EvtPDL::name( EvtId i )
{
    return getInstance().m_partlist[i.getAlias()].getName();
}

size_t EvtPDL::entries()
{
    return getInstance().m_partlist.size();
}

EvtId EvtPDL::getEntry( int i )
{
    return getInstance().m_partlist[i].getId();
}

void EvtPDL::reSetMass( EvtId i, double mass )
{
    getInstance().m_partlist[i.getId()].reSetMass( mass );
}

void EvtPDL::reSetWidth( EvtId i, double width )
{
    getInstance().m_partlist[i.getId()].reSetWidth( width );
}

void EvtPDL::reSetMassMin( EvtId i, double mass )
{
    getInstance().m_partlist[i.getId()].reSetMassMin( mass );
}

void EvtPDL::reSetMassMax( EvtId i, double mass )
{
    getInstance().m_partlist[i.getId()].reSetMassMax( mass );
}

void EvtPDL::reSetBlatt( EvtId i, double blatt )
{
    getInstance().m_partlist[i.getId()].reSetBlatt( blatt );
}

void EvtPDL::reSetBlattBirth( EvtId i, double blatt )
{
    getInstance().m_partlist[i.getId()].reSetBlattBirth( blatt );
}

void EvtPDL::includeBirthFactor( EvtId i, bool yesno )
{
    getInstance().m_partlist[i.getId()].includeBirthFactor( yesno );
}

void EvtPDL::includeDecayFactor( EvtId i, bool yesno )
{
    getInstance().m_partlist[i.getId()].includeDecayFactor( yesno );
}

void EvtPDL::changeLS( EvtId i, std::string& newLS )
{
    getInstance().m_partlist[i.getId()].newLineShape( newLS );
}

void EvtPDL::setPWForDecay( EvtId i, int spin, EvtId d1, EvtId d2 )
{
    getInstance().m_partlist[i.getId()].setPWForDecay( spin, d1, d2 );
}

void EvtPDL::setPWForBirthL( EvtId i, int spin, EvtId par, EvtId othD )
{
    getInstance().m_partlist[i.getId()].setPWForBirthL( spin, par, othD );
}
