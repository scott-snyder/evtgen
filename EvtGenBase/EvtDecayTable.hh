
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

#ifndef EVTDECAYTABLE_HH
#define EVTDECAYTABLE_HH

#include "EvtGenBase/EvtDecayBase.hh"
#include "EvtGenBase/EvtParticleDecayList.hh"

#include <vector>

class EvtId;

// Description: Class to read in and handle the decays available
//              to EvtGen for each particle, and the model to be
//              used for each one.

class EvtDecayTable {
  public:
    static EvtDecayTable* getInstance();

    int getNMode( int ipar ) const;

    EvtDecayBase* getDecay( int ipar, int imode );

    void readDecayFile( const std::string dec_name, bool verbose = true );
    void readXMLDecayFile( const std::string dec_name, bool verbose = true );

    bool stringToBoolean( std::string valStr ) const;
    void checkParticle( std::string particle ) const;

    int findChannel( EvtId parent, std::string model, int ndaug, EvtId* daugs,
                     int narg, std::string* args ) const;

    int inChannelList( EvtId parent, int ndaug, EvtId* daugs ) const;

    EvtDecayBase* getDecayFunc( EvtParticle* p );

    void printSummary() const;

    void checkConj() const;

    std::vector<EvtParticleDecayList> getDecayTable() const
    {
        return m_decaytable;
    }

    EvtDecayBase* findDecayModel( int aliasInt, int modeInt );
    EvtDecayBase* findDecayModel( EvtId id, int modeInt );

    bool hasPythia( int aliasInt ) const;
    bool hasPythia( EvtId id ) const;

    int getNModes( int aliasInt ) const;
    int getNModes( EvtId id ) const;

    std::vector<std::string> splitString( std::string& theString,
                                          std::string& splitter ) const;

  protected:
    EvtDecayTable();
    ~EvtDecayTable();

  private:
    std::vector<EvtParticleDecayList> m_decaytable;

    EvtDecayTable( const EvtDecayTable& ){};
    //EvtDecayTable& operator=(const EvtDecayTable&) {};
};

#endif
