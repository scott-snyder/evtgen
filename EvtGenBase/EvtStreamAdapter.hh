
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

#ifndef EVT_STREAM_ADAPTER_HH
#define EVT_STREAM_ADAPTER_HH

// Stream adapters are used to convert a stream-like input (for example,
// a file containing N entries) to an STL like iterator interface. There
// must be a way to get point from the stream, and also an indicator of the
// end of the stream.

template <class Point>
class EvtStreamAdapter {
  public:
    EvtStreamAdapter() {}
    virtual ~EvtStreamAdapter() {}
    virtual EvtStreamAdapter* clone() const = 0;
    virtual Point currentValue() = 0;
    virtual void advance() = 0;
    virtual bool pastEnd() = 0;
};

// N points are read from a generated stream.

template <class Point, class Generator>
class EvtGenStreamAdapter : public EvtStreamAdapter<Point> {
  public:
    EvtGenStreamAdapter( Generator gen, int count ) :
        m_gen( gen ), m_count( count )
    {
    }

    virtual ~EvtGenStreamAdapter() {}

    EvtStreamAdapter<Point>* clone() const override
    {
        return new EvtGenStreamAdapter( *this );
    }
    Point currentValue() override { return m_gen(); }
    bool pastEnd() override { return ( m_count <= 0 ); }
    void advance() override { m_count--; }

  private:
    Generator m_gen;
    int m_count;    // also serves as past the end indicator
};

// Only points satisfying a predicate are read from the stream.

template <class Point, class Iterator, class Predicate>
class EvtPredStreamAdapter : public EvtStreamAdapter<Point> {
  public:
    EvtPredStreamAdapter( Predicate pred, Iterator it, Iterator end ) :
        m_pred( pred ), m_it( it ), m_end( end )
    {
    }
    virtual ~EvtPredStreamAdapter() {}

    virtual EvtStreamAdapter<Point>* clone() const
    {
        return new EvtPredStreamAdapter( *this );
    }
    virtual Point currentValue()
    {
        Point value;
        while ( !pastEnd() ) {
            value = *m_it;
            if ( m_pred( value ) )
                break;
            m_it++;
        }
        return value;
    }

    virtual bool pastEnd() { return m_it == m_end; }
    virtual void advance() { m_it++; }

  private:
    Predicate m_pred;
    Iterator m_it;
    Iterator m_end;
};

#endif
