
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

#include "EvtGenBase/EvtParserXml.hh"

#include "EvtGenBase/EvtReport.hh"

#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <string.h>
#include <vector>
using namespace std;

bool EvtParserXml::open( std::string filename )
{
    if ( !expandEnvVars( filename ) ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Error while expanding environment variables in file name '"
            << filename.c_str() << "'" << endl;
        return false;
    }

    m_fin.open( filename.c_str() );
    if ( !m_fin ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Could not open file '" << filename.c_str() << "'" << endl;
        return false;
    }

    return true;
}

bool EvtParserXml::close()
{
    m_fin.close();
    return true;
}

bool EvtParserXml::readNextTag()
{
    if ( !processTagTree() ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Unexpected end tag " << m_tagTitle << " found near line "
            << m_lineNo << endl;
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Will terminate execution!" << endl;
        return false;
    }    //first process the previous tag to find out where we are in the tag tree

    while ( m_line.find( "<" ) ==
            std::string::npos ) {    //add lines until we find start of a tag
        std::string addLine;
        if ( !std::getline( m_fin, addLine ) )
            return false;
        m_lineNo++;
        m_line += " ";
        m_line += addLine;
    }

    unsigned int startTag;
    unsigned int endTag;
    unsigned int endTagTitle;

    startTag = m_line.find( "<" );

    if ( m_line[startTag + 1] ==
         '?' ) {    //XML header tag - ignore then read the next tag
        while ( m_line.find( "?>", startTag ) == std::string::npos ) {
            std::string addLine;
            if ( !std::getline( m_fin, addLine ) )
                return false;
            m_lineNo++;
            m_line += " ";
            m_line += addLine;
        }
        endTag = m_line.find( "?>", startTag );
        m_line = m_line.substr( endTag + 2 );
        return readNextTag();
    } else if ( m_line[startTag + 1] ==
                '!' ) {    //XML comment tag - ignore then read the next tag
        while ( m_line.find( "-->", startTag ) == std::string::npos ) {
            std::string addLine;
            if ( !std::getline( m_fin, addLine ) )
                return false;
            m_lineNo++;
            m_line += " ";
            m_line += addLine;
        }
        endTag = m_line.find( "-->", startTag );
        m_line = m_line.substr( endTag + 3 );
        m_tagTitle = "";
        m_tag = "";
        return readNextTag();
    } else {    //parsable

        while ( m_line.find( ">", startTag ) ==
                std::string::npos ) {    //find end of a tag
            std::string addLine;
            if ( !std::getline( m_fin, addLine ) )
                return false;
            m_lineNo++;
            m_line += " ";
            m_line += addLine;
        }
        endTag = m_line.find( ">", startTag );
        m_inLineTag = false;
        if ( m_line.find( "/>", startTag ) < endTag ) {
            endTag--;
            m_inLineTag = true;
        }

        if ( m_line.find( " ", startTag ) != std::string::npos &&
             m_line.find( " ", startTag ) <
                 endTag ) {    //find end of the first word in the tag
            endTagTitle = m_line.find( " ", startTag );
        } else {
            endTagTitle = endTag;
        }

        m_tagTitle = m_line.substr( startTag + 1, endTagTitle - startTag - 1 );
        m_tag = m_line.substr( startTag + 1, endTag - startTag - 1 );

        //now we have the tag lets remove it from the line
        if ( m_inLineTag ) {
            m_line = m_line.substr( endTag + 2 );
        } else {
            m_line = m_line.substr( endTag + 1 );
        }
        return true;
    }
}

std::string EvtParserXml::getParentTagTitle()
{
    if ( m_tagTree.empty() )
        return "";
    else
        return m_tagTree.back();
}

std::string EvtParserXml::readAttribute( std::string attribute,
                                         std::string defaultValue )
{
    std::string whitespace = " \t\n\v\f\r";
    for ( unsigned int i = 0; i < whitespace.size(); i++ ) {
        //find any whitespace followed by the attribute name followed by an '='
        std::string attName = whitespace[i] + attribute + "=";
        if ( m_tag.find( attName ) != std::string::npos ) {
            int startAttri = m_tag.find( attName );
            int startQuote = m_tag.find( "\"", startAttri + 1 );
            int endQuote = m_tag.find( "\"", startQuote + 1 );
            return m_tag.substr( startQuote + 1, endQuote - startQuote - 1 );
        }
    }
    return defaultValue;
}

bool EvtParserXml::readAttributeBool( std::string attribute, bool defaultValue )
{
    std::string valStr = readAttribute( attribute );
    if ( !defaultValue )
        return ( valStr == "true" || valStr == "1" || valStr == "on" ||
                 valStr == "yes" );
    else
        return ( valStr != "false" && valStr != "0" && valStr != "off" &&
                 valStr != "no" );
}

int EvtParserXml::readAttributeInt( std::string attribute, int defaultValue )
{
    std::string valStr = readAttribute( attribute );
    if ( valStr == "" )
        return defaultValue;
    std::istringstream valStream( valStr );
    int retVal;
    valStream >> retVal;
    return retVal;
}

double EvtParserXml::readAttributeDouble( std::string attribute,
                                          double defaultValue )
{
    std::string valStr = readAttribute( attribute );
    if ( valStr == "" )
        return defaultValue;
    std::istringstream valStream( valStr );
    double retVal;
    valStream >> retVal;
    return retVal;
}

bool EvtParserXml::processTagTree()
{
    if ( m_tagTitle == "" )
        return true;
    if ( m_tagTitle[0] == '/' ) {
        if ( m_tagTitle.substr( 1 ) == m_tagTree.back() ) {
            m_tagTree.pop_back();
        } else {
            return false;
        }
    } else if ( !m_inLineTag ) {
        m_tagTree.push_back( m_tagTitle );
    }
    return true;
}

bool EvtParserXml::expandEnvVars( std::string& str )
{
    while ( str.find( '$' ) != std::string::npos ) {
        size_t varStart = str.find( '$' );
        size_t varNameLength;
        std::string varName;

        //if this is the last character then just remove the $
        if ( varStart == str.length() - 1 ) {
            str.erase( varStart );
            return true;
        }

        if ( str[varStart + 1] == '{' ) {
            //deal with environment variables in {}s
            size_t braceStart = varStart + 1;
            size_t braceEnd = str.find( '}', braceStart );

            if ( braceEnd == std::string::npos ) {
                EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                    << "Incomplete environment variable found in text: " << str
                    << endl;
                EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                    << "Will terminate execution!" << endl;
                return false;
            }

            varName = str.substr( braceStart + 1, braceEnd - braceStart - 1 );
            varNameLength = braceEnd - varStart;

        } else {
            //deal with everything else
            varNameLength = 0;

            while ( varNameLength + varStart + 1 < str.length() &&
                    isAlphaNum( str[varStart + varNameLength + 1] ) ) {
                ++varNameLength;
            }

            varName = str.substr( varStart + 1, varNameLength );
        }

        char* envVar = getenv( varName.c_str() );

        if ( envVar )
            str.replace( varStart, varNameLength + 1, envVar );
        else {
            EvtGenReport( EVTGEN_WARNING, "EvtGen" )
                << "Undefined environment variable found in text: " << varName
                << endl;
            str.replace( varStart, varNameLength + 1, "" );
        }
    }
    return true;
}

bool EvtParserXml::isAlphaNum( char c )
{
    if ( c >= '0' && c <= '9' )
        return true;
    if ( c >= 'A' && c <= 'Z' )
        return true;
    if ( c >= 'a' && c <= 'z' )
        return true;
    if ( c == '_' )
        return true;
    return false;
}
