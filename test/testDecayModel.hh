
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

#ifndef TEST_DECAY_MODEL_HH
#define TEST_DECAY_MODEL_HH

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"

#include "nlohmann/json.hpp"

#include <cmath>
#include <string>
#include <utility>
#include <vector>

class EvtGen;
class EvtParticle;

class TestInfo {
  public:
    TestInfo( const std::string& name, const int d1, const int d2,
              const std::string& nameY = "", const int d1Y = 0,
              const int d2Y = 0 ) :
        m_name{ name }, m_nameY{ nameY }, m_d1{ d1 }, m_d2{ d2 }, m_d1Y{ d1Y }, m_d2Y{ d2Y }
    {
    }

    const std::string getName( const int var = 1 ) const
    {
        return var == 2 ? m_nameY : m_name;
    }
    int getd1( const int var = 1 ) const { return var == 2 ? m_d1Y : m_d1; }
    int getd2( const int var = 1 ) const { return var == 2 ? m_d2Y : m_d2; }

  private:
    std::string m_name;
    std::string m_nameY;
    int m_d1;
    int m_d2;
    int m_d1Y;
    int m_d2Y;
};

class TestDecayModel {
  public:
    TestDecayModel( const nlohmann::json& config );

    bool run();

  private:
    bool checkMandatoryFields();

    std::string createDecFile(
        const std::string& parent, const std::vector<std::string>& daughterNames,
        const std::vector<std::vector<std::string>>& grandDaughterNames,
        const std::vector<std::string>& models,
        const std::vector<std::vector<std::string>>& parameters,
        const std::vector<bool>& doConjDecay,
        const std::vector<std::string>& extras,
        const std::string decFileName ) const;

    void defineHistos( TFile* theFile );

    void generateEvents( EvtGen& theGen, const std::string& decFile,
                         const std::string& parentName, const bool doConjDecay,
                         const int nEvents, const bool debugFlag );

    int findChargedDaugtherWithMaxE( const EvtParticle* parent ) const;

    double getValue( const EvtParticle* rootPart, const std::string& varName,
                     const int d1, const int d2 ) const;

    void compareHistos( const std::string& refFileName ) const;

    double getCosAcoplanarityAngle( const EvtParticle* selectedParent,
                                    const int sel_NDaugMax, const int d1,
                                    const int d2 ) const;

    const nlohmann::json& m_config;
    std::vector<std::pair<TestInfo, TH1*>> m_1DhistVect;
    std::vector<std::pair<TestInfo, TH2*>> m_2DhistVect;
    TH1* m_mixedHist{ nullptr };
};

#endif
