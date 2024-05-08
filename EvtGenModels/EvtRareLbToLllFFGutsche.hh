
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

#ifndef EVTRARELBTOLLLFFGUTSCHE_HH
#define EVTRARELBTOLLLFFGUTSCHE_HH 1

// Include files

/** @class EvtRareLbToLllFF EvtRareLbToLllFF.hh EvtGenModels/EvtRareLbToLllFF.hh
 *
 *
 *  @author Michal Kreps
 *  @date   2014-10-21
 */

#include "EvtGenBase/EvtIdSet.hh"
#include "EvtGenBase/EvtParticle.hh"

#include "EvtGenModels/EvtRareLbToLllFFBase.hh"

#include <map>
#include <string>

class EvtRareLbToLllFFGutsche : public EvtRareLbToLllFFBase {
  public:
    void init() override;

    void getFF( const EvtParticle& parent, const EvtParticle& lambda,
                EvtRareLbToLllFFBase::FormFactors& FF ) const override;

  private:
    double formFactorParametrization( const double s, const double f0,
                                      const double a, const double b ) const;

    double m_Vconsts[3][3];
    double m_Aconsts[3][3];
    double m_TVconsts[3][3];
    double m_TAconsts[3][3];

    const EvtIdSet m_Parents{ "Lambda_b0", "anti-Lambda_b0" };
    const EvtIdSet m_Daughters{ "Lambda0", "anti-Lambda0" };
};

#endif    // EVTRARELBTOLLLFF_HH
