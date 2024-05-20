
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

#ifdef EVTGEN_SHERPA

#ifndef EVTSHERPARANDOM_HH
#define EVTSHERPARANDOM_HH

#include "EvtGenBase/EvtRandom.hh"

#include "ATOOLS/Math/Random.H"

/* Description: Class to specify the chosen EvtGen random number (engine)
 * to be used for Sherpa. 
 */

class EvtSherpaRandom : public ATOOLS::External_RNG {
  public:
    double Get() { return EvtRandom::Flat(); }

  private:
};

// The following code makes EvtSherpaRandom loadable in Sherpa

//! \cond DOXYGEN_IGNORE

// In Sherpa versions v2.2.15 and earlier the DECLARE_GETTER macro produces code that doesn't compile with gcc versions > 12 and C++20
// So, in those cases, instead use explicit code that has been modified to work correctly
#ifdef EVTGEN_SHERPA_USEMACRO
DECLARE_GETTER( EvtSherpaRandom, "EvtSherpaRandom", ATOOLS::External_RNG,
                ATOOLS::RNG_Key );
#else
namespace ATOOLS {
    template <>
    class Getter<External_RNG, RNG_Key, EvtSherpaRandom> :
        public Getter_Function<External_RNG, RNG_Key> {
      public:
        Getter( const bool& d = true ) :
            Getter_Function<External_RNG, RNG_Key>{ "EvtSherpaRandom" }
        {
            SetDisplay( d );
        }

      protected:
        void PrintInfo( std::ostream& str, const size_t width ) const;
        Object_Type* operator()( const Parameter_Type& parameters ) const;

      private:
        static Getter<External_RNG, RNG_Key, EvtSherpaRandom> s_initializer;
    };
}    // namespace ATOOLS

ATOOLS::Getter<ATOOLS::External_RNG, ATOOLS::RNG_Key, EvtSherpaRandom>
    ATOOLS::Getter<ATOOLS::External_RNG, ATOOLS::RNG_Key, EvtSherpaRandom>::s_initializer{
        true };
#endif

ATOOLS::External_RNG*
ATOOLS::Getter<ATOOLS::External_RNG, ATOOLS::RNG_Key, EvtSherpaRandom>::operator()(
    const ATOOLS::RNG_Key& ) const
{
    return new EvtSherpaRandom();
}

// This eventually prints a help message
void ATOOLS::Getter<ATOOLS::External_RNG, ATOOLS::RNG_Key, EvtSherpaRandom>::PrintInfo(
    std::ostream& str, const size_t ) const
{
    str << " EvtGen-Sherpa RNG interface";
}

//! \endcond

#endif

#endif
