
/***********************************************************************
* Copyright 1998-2021 CERN for the benefit of the EvtGen authors       *
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

#ifndef EVTIDSET_HH
#define EVTIDSET_HH

#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtPatches.hh"

#include <initializer_list>
#include <set>
#include <string>

class EvtIdSet {
  public:
    explicit EvtIdSet( std::initializer_list<EvtId> idList );
    explicit EvtIdSet( std::initializer_list<std::string> nameList );
    EvtIdSet( const EvtIdSet& set1, const EvtIdSet& set2 );

    std::size_t sizeOfSet() const { return m_list.size(); };
    bool contains( const EvtId& id ) const;
    bool contains( const std::string& name ) const;

    void append( const EvtIdSet& set1 );

  private:
    std::set<EvtId> m_list;
};

#endif
