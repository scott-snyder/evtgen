
/***********************************************************************
* Copyright 1998-2022 CERN for the benefit of the EvtGen authors       *
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

#ifndef EVTRARELBTOLLLFFLQCD_HH
#define EVTRARELBTOLLLFFLQCD_HH 1

// Include files

/** @class EvtRareLbToLllFF EvtRareLbToLllFF.hh EvtGenModels/EvtRareLbToLllFF.hh
 *
 *
 *  @author Michal Kreps
 *  @date   2016-04-19
 *  @date   2014-10-23
 */

#include "EvtGenBase/EvtIdSet.hh"
#include "EvtGenBase/EvtParticle.hh"

#include "EvtGenModels/EvtRareLbToLllFFBase.hh"

#include <map>
#include <string>

class EvtRareLbToLllFFlQCD : public EvtRareLbToLllFFBase {
  public:
    /// Standard constructor
    EvtRareLbToLllFFlQCD() = default;

    void init() override;

    void getFF( const EvtParticle& parent, const EvtParticle& lambda,
                EvtRareLbToLllFFBase::FormFactors& FF ) const override;

  private:
    double formFactorParametrization( const double q2, const double a0,
                                      const double a1, const double pole ) const;
    double zvar( const double q2 ) const;

    double m_fconsts[3][3];
    double m_gconsts[3][3];
    double m_hconsts[3][3];
    double m_htildaconsts[3][3];

    double m_t0;
    double m_tplus;
};

#endif    // EVTRARELBTOLLLFF_HH
