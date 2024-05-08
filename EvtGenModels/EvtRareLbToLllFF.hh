
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

#ifndef EVTRARELBTOLLLFF_HH
#define EVTRARELBTOLLLFF_HH 1

// Include files

/** @class EvtRareLbToLllFF EvtRareLbToLllFF.hh EvtGenModels/EvtRareLbToLllFF.hh
 *
 *
 *  @author Thomas Blake
 *  @date   2013-11-26
 */

#include "EvtGenBase/EvtIdSet.hh"
#include "EvtGenBase/EvtParticle.hh"

#include "EvtGenModels/EvtRareLbToLllFFBase.hh"

#include <array>
#include <map>
#include <memory>
#include <string>

class EvtRareLbToLllFF : public EvtRareLbToLllFFBase {
  public:
    class FormFactorDependence final {
      public:
        FormFactorDependence();

        FormFactorDependence( const double al, const double ap );

        FormFactorDependence( const double a0, const double a2, const double a4,
                              const double al, const double ap );

        FormFactorDependence( const FormFactorDependence& other );

        FormFactorDependence* clone() const;

        void param( const double al, const double ap );

        void param( const double a0, const double a2, const double a4,
                    const double al, const double ap );

        double m_a0;
        double m_a2;
        double m_a4;
        double m_al;
        double m_ap;
    };

    class FormFactorSet final {
      public:
        FormFactorSet();

        FormFactorSet( const FormFactorSet& other );

        EvtRareLbToLllFF::FormFactorDependence m_F1;
        EvtRareLbToLllFF::FormFactorDependence m_F2;
        EvtRareLbToLllFF::FormFactorDependence m_F3;
        EvtRareLbToLllFF::FormFactorDependence m_F4;

        EvtRareLbToLllFF::FormFactorDependence m_G1;
        EvtRareLbToLllFF::FormFactorDependence m_G2;
        EvtRareLbToLllFF::FormFactorDependence m_G3;
        EvtRareLbToLllFF::FormFactorDependence m_G4;

        EvtRareLbToLllFF::FormFactorDependence m_H1;
        EvtRareLbToLllFF::FormFactorDependence m_H2;
        EvtRareLbToLllFF::FormFactorDependence m_H3;
        EvtRareLbToLllFF::FormFactorDependence m_H4;
        EvtRareLbToLllFF::FormFactorDependence m_H5;
        EvtRareLbToLllFF::FormFactorDependence m_H6;
    };

    void init() override;

    void getFF( const EvtParticle& parent, const EvtParticle& lambda,
                EvtRareLbToLllFFBase::FormFactors& FF ) const override;

  private:
    double func( const double p,
                 const EvtRareLbToLllFF::FormFactorDependence& dep ) const;

    std::array<std::unique_ptr<EvtRareLbToLllFF::FormFactorSet>, 2> m_FF;
    std::map<int, EvtRareLbToLllFF::FormFactorSet*> m_FFMap;

    void DiracFF( const EvtParticle& parent, const EvtParticle& lambda,
                  const EvtRareLbToLllFF::FormFactorSet& FFset,
                  EvtRareLbToLllFF::FormFactors& FF ) const;

    void RaritaSchwingerFF( const EvtParticle& parent, const EvtParticle& lambda,
                            const EvtRareLbToLllFF::FormFactorSet& FFset,
                            EvtRareLbToLllFF::FormFactors& FF ) const;
};

#endif    // EVTRARELBTOLLLFF_HH
