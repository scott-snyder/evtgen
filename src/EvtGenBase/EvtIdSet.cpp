
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

#include "EvtGenBase/EvtIdSet.hh"

#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtPatches.hh"

#include <initializer_list>
#include <set>
#include <string>

EvtIdSet::EvtIdSet( std::initializer_list<EvtId> idList ) : m_list{ idList }
{
}

EvtIdSet::EvtIdSet( std::initializer_list<std::string> nameList )
{
    for ( const std::string& name : nameList ) {
        m_list.insert( EvtPDL::getId( name ) );
    }
}

EvtIdSet::EvtIdSet( const EvtIdSet& set1, const EvtIdSet& set2 ) :
    m_list{ set1.m_list }
{
    this->append( set2 );
}

bool EvtIdSet::contains( const EvtId& id ) const
{
    const auto iter = m_list.find( id );
    return iter != m_list.end();
}

bool EvtIdSet::contains( const std::string& name ) const
{
    const EvtId id{ EvtPDL::getId( name ) };
    return this->contains( id );
}

void EvtIdSet::append( const EvtIdSet& set1 )
{
    m_list.insert( set1.m_list.begin(), set1.m_list.end() );
}
