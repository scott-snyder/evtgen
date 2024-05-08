
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

#ifndef EVTSTRINGHASH_HH
#define EVTSTRINGHASH_HH

#include <string>

template <class T>
class EvtStringHash {
  public:
    inline EvtStringHash( int size );
    inline void add( const std::string& str, T* data );
    inline T* get( const std::string& str );
    inline ~EvtStringHash();

  private:
    EvtStringHash();
    int m_size;
    inline int hash( const std::string& str );
    std::string*** m_strings;
    T*** m_data;
    int* m_entries;
};

template <class T>
EvtStringHash<T>::EvtStringHash( int size )
{
    m_size = size;

    typedef std::string** EvtStringPtrPtr;
    typedef T** TPtrPtr;

    m_strings = new EvtStringPtrPtr[m_size];
    m_data = new TPtrPtr[m_size];
    m_entries = new int[m_size];

    int i;

    for ( i = 0; i < m_size; i++ ) {
        m_entries[i] = 0;
    }
}

template <class T>
EvtStringHash<T>::~EvtStringHash()
{
    int i;
    for ( i = 0; i < m_size; i++ ) {
        int j;
        for ( j = 0; j < m_entries[i]; j++ ) {
            delete m_strings[i][j];
        }
        if ( m_entries[i] > 0 ) {
            delete[] m_strings[i];
            delete[] m_data[i];
        }
    }

    delete[] m_strings;
    delete[] m_data;
    delete[] m_entries;
}

template <class T>
void EvtStringHash<T>::add( const std::string& str, T* data )
{
    int ihash = hash( str );

    typedef std::string* EvtStringPtr;
    typedef T* TPtr;

    std::string** newstrings = new EvtStringPtr[m_entries[ihash] + 1];
    T** newdata = new TPtr[m_entries[ihash] + 1];

    int i;

    for ( i = 0; i < m_entries[ihash]; i++ ) {
        newstrings[i] = m_strings[ihash][i];
        newdata[i] = m_data[ihash][i];
    }

    newstrings[m_entries[ihash]] = new std::string;
    *( newstrings[m_entries[ihash]] ) = str;
    newdata[m_entries[ihash]] = data;

    if ( m_entries[ihash] != 0 ) {
        delete[] m_strings[ihash];
        delete[] m_data[ihash];
    }

    m_entries[ihash]++;

    m_strings[ihash] = newstrings;
    m_data[ihash] = newdata;
}

template <class T>
T* EvtStringHash<T>::get( const std::string& str )
{
    int ihash = hash( str );

    int i;

    for ( i = 0; i < m_entries[ihash]; i++ ) {
        if ( *( m_strings[ihash][i] ) == str )
            return m_data[ihash][i];
    }

    return 0;
}

template <class T>
int EvtStringHash<T>::hash( const std::string& str )
{
    const char* cstr = str.c_str();

    int i = 0;

    int value = 0;

    while ( cstr[i] != 0 ) {
        value += (int)cstr[i];
        i++;
    }

    return value % m_size;
}

#endif
