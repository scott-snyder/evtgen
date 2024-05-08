
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

#include "EvtGenBase/EvtModelAlias.hh"

#include "EvtGenBase/EvtPatches.hh"

EvtModelAlias::EvtModelAlias( std::string alias, std::string model,
                              std::vector<std::string> args ) :

    m_aliasName( alias ), m_model( model ), m_modelArgs( args )

{
}

EvtModelAlias::EvtModelAlias( const EvtModelAlias& copyMe ) :

    m_aliasName( copyMe.m_aliasName ),
    m_model( copyMe.m_model ),
    m_modelArgs( copyMe.m_modelArgs )

{
}

EvtModelAlias EvtModelAlias::operator=( const EvtModelAlias& copyMe )
{
    m_aliasName = copyMe.m_aliasName;
    m_model = copyMe.m_model;
    m_modelArgs = copyMe.m_modelArgs;

    return *this;
}

std::vector<std::string> EvtModelAlias::getArgList()
{
    return m_modelArgs;
}
