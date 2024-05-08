
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

#include "EvtGenBase/EvtParser.hh"

#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"

#include <fstream>
#include <sstream>
#include <string.h>
using namespace std;

#define MAXBUF 1024

EvtParser::EvtParser()
{
    m_ntoken = 0;
    m_lengthoftokenlist = 0;
    m_tokenlist = nullptr;
    m_linelist = nullptr;
}

EvtParser::~EvtParser()
{
    delete[] m_tokenlist;
    delete[] m_linelist;
}

int EvtParser::getNToken()
{
    return m_ntoken;
}

const std::string& EvtParser::getToken( int i )
{
    return m_tokenlist[i];
}

int EvtParser::getLineofToken( int i )
{
    return m_linelist[i];
}

int EvtParser::read( const std::string filename )
{
    ifstream fin;

    fin.open( filename.c_str() );
    if ( !fin ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Could not open file '" << filename.c_str() << "'" << endl;
        return -1;
    }

    char buf[MAXBUF];
    char buf2[MAXBUF];
    char c;

    int line = 0;
    int i;

    while ( fin.peek() != EOF ) {
        line++;

        i = 0;
        while ( ( c = fin.get() ) != '\n' && c != EOF && i < MAXBUF ) {
            buf[i] = c;
            i++;
        }
        if ( i == MAXBUF ) {
            EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                << "Error in EvtParser: line:" << line << " to long" << endl;
        } else {
            buf[i] = '\0';
        }

        //search for '#' which indicates comment for rest of line!
        i = 0;
        do {
            if ( buf[i] == '#' )
                buf[i] = 0;
            i++;
        } while ( buf[i - 1] != 0 );

        string tmp( buf, strlen( buf ) );

        //read each token
        istringstream ist( tmp );
        while ( ist >> buf2 ) {
            i = 0;
            int semicolon = 0;
            do {
                if ( buf2[i] == ';' ) {
                    buf2[i] = 0;
                    semicolon = 1;
                }
            } while ( buf2[i++] != 0 );
            if ( buf2[0] != 0 ) {
                addToken( line, buf2 );
            }
            if ( semicolon )
                addToken( line, ";" );
        }
    }

    fin.close();

    return 0;
}

void EvtParser::addToken( int line, const std::string& string )
{
    //EvtGenReport(EVTGEN_INFO,"EvtGen") <<m_ntoken<<" "<<line<<" "<<string<<endl;

    if ( m_ntoken == m_lengthoftokenlist ) {
        int new_length = 1000 + 4 * m_lengthoftokenlist;

        int* newlinelist = new int[new_length];
        std::string* newtokenlist = new std::string[new_length];

        int i;

        for ( i = 0; i < m_ntoken; i++ ) {
            newlinelist[i] = m_linelist[i];
            newtokenlist[i] = m_tokenlist[i];
        }

        delete[] m_tokenlist;
        delete[] m_linelist;

        m_tokenlist = newtokenlist;
        m_linelist = newlinelist;

        m_lengthoftokenlist = new_length;
    }

    m_tokenlist[m_ntoken] = string;

    m_linelist[m_ntoken] = line;

    m_ntoken++;

    //EvtGenReport(EVTGEN_INFO,"EvtGen") << "First:"<<m_tokenlist[0]<<" last:"<<m_tokenlist[m_ntoken-1]<<endl;
}
