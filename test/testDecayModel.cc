
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

#include "testDecayModel.hh"

#include "EvtGen/EvtGen.hh"

#include "EvtGenBase/EvtAbsRadCorr.hh"
#include "EvtGenBase/EvtConst.hh"
#include "EvtGenBase/EvtDecayBase.hh"
#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtKine.hh"
#include "EvtGenBase/EvtMTRandomEngine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtParticleFactory.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtVector4R.hh"

#ifdef EVTGEN_EXTERNAL
#include "EvtGenExternal/EvtExternalGenList.hh"
#include "EvtGenExternal/EvtSherpaPhotons.hh"
#endif

#include "TROOT.h"

#include "tbb/blocked_range.h"
#include "tbb/global_control.h"
#include "tbb/parallel_reduce.h"

#include <chrono>
#include <fstream>
#include <future>
#include <iostream>
#include <list>
#include <memory>

using nlohmann::json;

std::once_flag TestDecayModel::m_createDecFile_threadlock;

TestHistos::TestHistos( const std::string& parentName, const json& config )
{
    // Histogram information
    const std::size_t nHistos{ config.size() };

    m_1DhistVect.reserve( nHistos );
    m_2DhistVect.reserve( nHistos );

    for ( const auto& hInfo : config ) {
        const auto varTitle{ hInfo.at( "title" ).get<std::string>() };

        const auto varName{ hInfo.at( "variable" ).get<std::string>() };
        // Integer values that define what particles need to be used
        // for invariant mass combinations or helicity angles etc
        const auto d1{ hInfo.at( "d1" ).get<int>() };
        const auto d2{ hInfo.at( "d2" ).get<int>() };

        const auto nBins{ hInfo.at( "nbins" ).get<int>() };
        const auto xmin{ hInfo.at( "xmin" ).get<double>() };
        const auto xmax{ hInfo.at( "xmax" ).get<double>() };

        std::string histName( varName.c_str() );
        if ( d1 != 0 ) {
            histName += "_";
            histName += std::to_string( d1 );
        }
        if ( d2 != 0 ) {
            histName += "_";
            histName += std::to_string( d2 );
        }

        if ( !hInfo.contains( "variableY" ) ) {
            auto hist = std::make_unique<TH1D>( histName.c_str(),
                                                varTitle.c_str(), nBins, xmin,
                                                xmax );
            m_1DhistVect.emplace_back(
                std::make_pair( HistInfo{ varName, d1, d2 }, std::move( hist ) ) );
        } else {
            const auto varNameY{ hInfo.at( "variableY" ).get<std::string>() };
            const auto d1Y{ hInfo.at( "d1Y" ).get<int>() };
            const auto d2Y{ hInfo.at( "d2Y" ).get<int>() };

            const auto nBinsY{ hInfo.at( "nbinsY" ).get<int>() };
            const auto ymin{ hInfo.at( "ymin" ).get<double>() };
            const auto ymax{ hInfo.at( "ymax" ).get<double>() };

            histName += "_";
            histName += varNameY;
            if ( d1Y != 0 ) {
                histName += "_";
                histName += std::to_string( d1Y );
            }
            if ( d2Y != 0 ) {
                histName += "_";
                histName += std::to_string( d2Y );
            }
            auto hist = std::make_unique<TH2D>( histName.c_str(),
                                                varTitle.c_str(), nBins, xmin,
                                                xmax, nBinsY, ymin, ymax );
            m_2DhistVect.emplace_back(
                std::make_pair( HistInfo{ varName, d1, d2, varNameY, d1Y, d2Y },
                                std::move( hist ) ) );
        }
    }

    // Add a mixed/unmixed histogram
    // Useful for the case where the parent is a neutral K, D or B
    const std::array<std::string, 14> parentsThatMix{
        "B_s0", "anti-B_s0", "B_s0L",   "B_s0H", "B0",      "anti-B0", "B0L",
        "B0H",  "D0",        "anti-D0", "K0",    "anti-K0", "K_S0",    "K_L0" };
    if ( std::find( parentsThatMix.begin(), parentsThatMix.end(), parentName ) !=
         parentsThatMix.end() ) {
        const std::string varTitle{ parentName + " mixed" };
        m_mixedHist = std::make_unique<TH1D>( "mixed", varTitle.c_str(), 2, 0.0,
                                              2.0 );
        // TODO maybe set bin labels?
    }
}

TestHistos::TestHistos( const TestHistos& rhs )
{
    m_1DhistVect.reserve( rhs.m_1DhistVect.size() );
    for ( auto& [info, hist] : rhs.m_1DhistVect ) {
        auto newHist = std::unique_ptr<TH1>{ static_cast<TH1*>( hist->Clone() ) };
        m_1DhistVect.push_back( std::make_pair( info, std::move( newHist ) ) );
    }

    m_2DhistVect.reserve( rhs.m_2DhistVect.size() );
    for ( auto& [info, hist] : rhs.m_2DhistVect ) {
        auto newHist = std::unique_ptr<TH2>{ static_cast<TH2*>( hist->Clone() ) };
        m_2DhistVect.push_back( std::make_pair( info, std::move( newHist ) ) );
    }

    if ( rhs.m_mixedHist ) {
        m_mixedHist.reset( static_cast<TH1*>( rhs.m_mixedHist->Clone() ) );
    }
}

TestHistos::TestHistos( TestHistos&& rhs ) noexcept
{
    this->swap( rhs );
}

TestHistos& TestHistos::operator=( const TestHistos& rhs )
{
    TestHistos tmp{ rhs };
    this->swap( tmp );
    return *this;
}

TestHistos& TestHistos::operator=( TestHistos&& rhs ) noexcept
{
    this->swap( rhs );
    return *this;
}

void TestHistos::swap( TestHistos& rhs ) noexcept
{
    m_1DhistVect.swap( rhs.m_1DhistVect );
    m_2DhistVect.swap( rhs.m_2DhistVect );
    std::swap( m_mixedHist, rhs.m_mixedHist );
}

void TestHistos::add( const TestHistos& rhs )
{
    // handle the special case where we have been default constructed and the rhs has not
    if ( m_1DhistVect.empty() && m_2DhistVect.empty() && !m_mixedHist ) {
        ( *this ) = rhs;
        return;
    }

    // TODO - should really check that the sets of histograms are the same between left and right

    const std::size_t n1DHists{ rhs.m_1DhistVect.size() };
    for ( std::size_t i{ 0 }; i < n1DHists; ++i ) {
        m_1DhistVect[i].second->Add( rhs.m_1DhistVect[i].second.get() );
    }

    const std::size_t n2DHists{ rhs.m_2DhistVect.size() };
    for ( std::size_t i{ 0 }; i < n2DHists; ++i ) {
        m_2DhistVect[i].second->Add( rhs.m_2DhistVect[i].second.get() );
    }

    if ( m_mixedHist && rhs.m_mixedHist ) {
        m_mixedHist->Add( rhs.m_mixedHist.get() );
    }
}

void TestHistos::normalise()
{
    for ( auto& [_, hist] : m_1DhistVect ) {
        const double area{ hist->Integral() };
        if ( area > 0.0 ) {
            hist->Scale( 1.0 / area );
        }
    }
    for ( auto& [_, hist] : m_2DhistVect ) {
        const double area{ hist->Integral() };
        if ( area > 0.0 ) {
            hist->Scale( 1.0 / area );
        }
    }
    if ( m_mixedHist ) {
        const double area{ m_mixedHist->Integral() };
        if ( area > 0.0 ) {
            m_mixedHist->Scale( 1.0 / area );
        }
    }
}

void TestHistos::save( TFile* outputFile )
{
    outputFile->cd();

    for ( auto& [info, hist] : m_1DhistVect ) {
        hist->SetDirectory( outputFile );
        hist->Write();
        hist.release();
    }

    for ( auto& [info, hist] : m_2DhistVect ) {
        hist->SetDirectory( outputFile );
        hist->Write();
        hist.release();
    }

    if ( m_mixedHist ) {
        m_mixedHist->SetDirectory( outputFile );
        m_mixedHist->Write();
        m_mixedHist.release();
    }
}

TestDecayModel::TestDecayModel( const json& config ) :
    m_config{ checkMandatoryFields( config )
                  ? readConfig( config )
                  : throw std::runtime_error{
                        "ERROR : json does not contain all required fields" } }
{
}

bool TestDecayModel::checkMandatoryFields( const json& config )
{
    const std::array<std::string, 7> mandatoryFields{ "parent",    "daughters",
                                                      "models",    "parameters",
                                                      "outfile",   "events",
                                                      "histograms" };
    const std::array<std::string, 7> mandatoryHistoFields{
        "title", "variable", "d1", "d2", "nbins", "xmin", "xmax" };
    const std::array<std::string, 6> extra2DHistoFields{ "variableY", "d1Y",
                                                         "d2Y",       "nbinsY",
                                                         "ymin",      "ymax" };

    bool allMandatoryFields{ true };

    for ( const auto& field : mandatoryFields ) {
        if ( !config.contains( field ) ) {
            std::cerr << "ERROR : json does not contain required field: " << field
                      << std::endl;
            allMandatoryFields = false;
            continue;
        }
        if ( field == "histograms" ) {
            const json& jHistos{ config.at( "histograms" ) };
            for ( const auto& hInfo : jHistos ) {
                for ( const auto& hField : mandatoryHistoFields ) {
                    if ( !hInfo.contains( hField ) ) {
                        std::cerr
                            << "ERROR : json does not contain required field for histogram definition: "
                            << hField << std::endl;
                        allMandatoryFields = false;
                    }
                }
                if ( hInfo.contains( extra2DHistoFields[0] ) ) {
                    for ( const auto& hField : extra2DHistoFields ) {
                        if ( !hInfo.contains( hField ) ) {
                            std::cerr
                                << "ERROR : json does not contain required field for 2D histogram definition: "
                                << hField << std::endl;
                            allMandatoryFields = false;
                        }
                    }
                }
            }
        }
    }

    return allMandatoryFields;
}

TestConfig TestDecayModel::readConfig( const json& config )
{
    TestConfig cfg;

    // Get all the mandatory fields first
    cfg.parentName = config.at( "parent" ).get<std::string>();
    cfg.daughterNames = config.at( "daughters" ).get<std::vector<std::string>>();
    cfg.modelNames = config.at( "models" ).get<std::vector<std::string>>();
    cfg.modelParameters =
        config.at( "parameters" ).get<std::vector<std::vector<std::string>>>();
    cfg.nEvents = config.at( "events" ).get<std::size_t>();

    // Histogram information
    cfg.testHistograms = TestHistos{ cfg.parentName, config.at( "histograms" ) };

    // Then check for optional fields, setting default values if not present

    if ( config.contains( "grand_daughters" ) &&
         config.at( "grand_daughters" ).is_array() ) {
        cfg.grandDaughterNames = config.at( "grand_daughters" )
                                     .get<std::vector<std::vector<std::string>>>();
    }

    if ( config.contains( "extras" ) && config.at( "extras" ).is_array() ) {
        cfg.extraCommands = config.at( "extras" ).get<std::vector<std::string>>();
    }

    // Set the number of threads to use, 1 by default
    cfg.nThreads = 1;
    if ( config.contains( "threads" ) ) {
        cfg.nThreads = config.at( "threads" ).get<std::size_t>();
    }

    // Set the FSR generator. Use PHOTOS by default.
    cfg.fsrGenerator = FSRGenerator::PHOTOS;
    if ( config.contains( "fsr_generator" ) ) {
        cfg.fsrGenerator = config.at( "fsr_generator" ).get<FSRGenerator>();
    }

    // Set the type of threading to use, stdlib by default
    cfg.threadModel = ThreadModel::StdLib;
    if ( config.contains( "thread_model" ) ) {
        cfg.threadModel = config.at( "thread_model" ).get<ThreadModel>();
    }

    // Set the RNG seed base (defaults to zero), to which the event number is added
    cfg.rngSeed = 0;
    if ( config.contains( "rng_seed" ) ) {
        cfg.rngSeed = config.at( "rng_seed" ).get<std::size_t>();
    }

    // Set reference and output file names
    // Insert fsrGenerator name if FSR simulation is not deactivated
    const bool noFSR = std::find( cfg.extraCommands.begin(),
                                  cfg.extraCommands.end(),
                                  "noFSR" ) != cfg.extraCommands.end();

    const std::string fileNameEnd =
        noFSR ? ".root" : "_" + to_string( cfg.fsrGenerator ) + ".root";

    const auto outFileStrSize = config.at( "outfile" ).get<std::string>().size() -
                                5;

    cfg.outFileName =
        config.at( "outfile" ).get<std::string>().substr( 0, outFileStrSize ) +
        fileNameEnd;
    cfg.refFileName = "Ref/" + cfg.outFileName;
    cfg.decFileName = cfg.outFileName.substr( 0, cfg.outFileName.size() - 5 ) +
                      ".dec";

    cfg.debugFlag = ( config.contains( "debug_flag" ) &&
                      config.at( "debug_flag" ).is_boolean() )
                      ? config.at( "debug_flag" ).get<bool>()
                      : false;

    if ( config.contains( "do_conjugate_decay" ) &&
         config.at( "do_conjugate_decay" ).is_array() ) {
        cfg.doConjDecay =
            config.at( "do_conjugate_decay" ).get<std::vector<bool>>();
    }
    if ( cfg.doConjDecay.size() != cfg.modelNames.size() ) {
        cfg.doConjDecay.resize( cfg.modelNames.size(), false );
    }

    return cfg;
}

void TestDecayModel::run()
{
    TestHistos theHistos;

    const auto start{ std::chrono::steady_clock::now() };

    if ( m_config.nThreads > 1 ) {
        // Run multi-threaded using the specified thread model
        switch ( m_config.threadModel ) {
            case ThreadModel::StdLib:
                theHistos = runStdThreads();
                break;
            case ThreadModel::TBB:
                theHistos = runTBBThreads();
                break;
        }
    } else {
        // Run in the main thread
        theHistos = runDecayBody( 0, m_config.nEvents );
    }

#ifdef EVTGEN_SHERPA
    // The Sherpa and PHOTONS++ instances are static global instances inside EvtGen
    // as they are not thread safe (they are mutexed). From Sherpa 3.0.0 on,
    // the internal Sherpa settings are also a static global instance.
    // Because of the static destruction order fiasco,
    // we need to control the order in which the static instances are destroyed.
    // This is done by hand inside the function below. If the function below is not called,
    // a SegFault occurs at the end of the program execution when using Sherpa for FSR.
    EvtSherpaPhotons::finalise();
#endif

    const auto end{ std::chrono::steady_clock::now() };
    const std::chrono::duration<double, std::milli> elapsed_ms{ ( end - start ) };
    const std::chrono::duration<double, std::milli> elapsed_ms_per_event{
        elapsed_ms / m_config.nEvents };

    std::cout << "Took " << elapsed_ms.count() << " ms to generate "
              << m_config.nEvents << " events using " << m_config.nThreads
              << " " << to_string( m_config.threadModel ) << " threads ("
              << elapsed_ms_per_event.count() << " ms per event)" << std::endl;

    // Normalise histograms.
    theHistos.normalise();

    // Compare with reference histograms
    compareHistos( theHistos, m_config.refFileName );

    // Create the root output file and write the histograms to it
    // Only save the mixed/unmixed histogram for neutral K, D, B
    std::unique_ptr<TFile> outFile{
        TFile::Open( m_config.outFileName.c_str(), "recreate" ) };
    theHistos.save( outFile.get() );
    outFile->Close();
    std::cout << "Created output file: " << m_config.outFileName.c_str()
              << std::endl;
}

TestHistos TestDecayModel::runStdThreads() const
{
    // Determine the number of threads and the number of events per thread
    const std::size_t nThreads{ std::min( m_config.nThreads, m_config.nEvents ) };
    const std::size_t nEventsPerThread{ ( m_config.nEvents % nThreads )
                                            ? m_config.nEvents / nThreads + 1
                                            : m_config.nEvents / nThreads };

    // Create the store for the results from each thread
    std::vector<std::future<TestHistos>> allHistos;
    allHistos.reserve( nThreads );

    // Launch the threads
    std::size_t firstEvent{ 0 };
    std::size_t nEvents{ nEventsPerThread };

    for ( std::size_t iThread{ 0 }; iThread < nThreads; ++iThread ) {
        // The last thread may need to generate slightly fewer events to get the required total
        if ( ( firstEvent + nEvents ) >= m_config.nEvents ) {
            nEvents = m_config.nEvents - firstEvent;
        }

        std::cout << "Thread " << iThread << " will generate " << nEvents
                  << " events" << std::endl;

        allHistos.emplace_back(
            std::async( std::launch::async, [this, firstEvent, nEvents]()
                        { return runDecayBody( firstEvent, nEvents ); } ) );

        firstEvent += nEvents;
    }

    // Now wait for each thread to finish
    bool complete{ false };
    do {
        // Set the flag to completed and set it back if we find incomplete threads
        complete = true;
        for ( auto& future : allHistos ) {
            auto status = future.wait_for( std::chrono::seconds( 10 ) );
            if ( status != std::future_status::ready ) {
                complete = false;
            }
        }
    } while ( !complete );

    // Accumulate the histograms from all threads
    TestHistos theHistos{ allHistos[0].get() };
    for ( std::size_t iThread{ 1 }; iThread < nThreads; ++iThread ) {
        theHistos.add( allHistos[iThread].get() );
    }

    return theHistos;
}

TestHistos TestDecayModel::runTBBThreads() const
{
    tbb::global_control gc{ tbb::global_control::parameter::max_allowed_parallelism,
                            m_config.nThreads };

    TestHistos null_hist;

    return tbb::parallel_reduce(
        tbb::blocked_range<std::size_t>( 0, m_config.nEvents ), null_hist,
        [this]( const tbb::blocked_range<std::size_t>& range,
                const TestHistos& init ) -> TestHistos
        {
            std::cout << "Thread "
                      << tbb::this_task_arena::current_thread_index()
                      << " will generate " << range.size()
                      << " events: " << range.begin() << " - " << range.end()
                      << std::endl;
            TestHistos tmp{ init };
            tmp.add( runDecayBody( range.begin(), range.size() ) );
            return tmp;
        },
        []( const TestHistos& lhs, const TestHistos& rhs ) -> TestHistos
        {
            TestHistos tmp{ lhs };
            tmp.add( rhs );
            return tmp;
        } );
}

TestHistos TestDecayModel::runDecayBody( const std::size_t firstEvent,
                                         const std::size_t nEvents ) const
{
    // Initialise the EvtGen object and hence the EvtPDL tables
    // The EvtGen object is used by generateEvents, while the
    // latter are also used within createDecFile

    // Define the random number generator
    static thread_local auto randomEngine{ std::make_unique<EvtMTRandomEngine>() };

    // TODO - need to streamline the extra models stuff
    static thread_local EvtAbsRadCorr* radCorrEngine{ nullptr };
    std::list<EvtDecayBase*> extraModels;

    static thread_local bool initialised{ false };

#ifdef EVTGEN_EXTERNAL
    if ( !initialised ) {
        bool convertPythiaCodes( false );
        bool useEvtGenRandom( true );
        bool seedTauolaFortran( true );
        EvtExternalGenList genList( convertPythiaCodes, "", "gamma",
                                    useEvtGenRandom, seedTauolaFortran );
        switch ( m_config.fsrGenerator ) {
            case FSRGenerator::PHOTOS:
                radCorrEngine = genList.getPhotosModel();
                break;
            case FSRGenerator::SherpaPhotons1:
                radCorrEngine = genList.getSherpaPhotonsModel( 1e-7, 1, 0 );
                break;
            case FSRGenerator::SherpaPhotons20:
                radCorrEngine = genList.getSherpaPhotonsModel( 1e-7, 2, 0 );
                break;
            case FSRGenerator::SherpaPhotons21:
                radCorrEngine = genList.getSherpaPhotonsModel( 1e-7, 2, 1 );
                break;
            case FSRGenerator::VinciaQED:
                radCorrEngine = genList.getVinciaQEDModel( 1.0e-7 );
        }
        extraModels = genList.getListOfModels();
    }
#endif

    static thread_local EvtGen theGen{ "../DECAY.DEC", "../evt.pdl",
                                       randomEngine.get(), radCorrEngine,
                                       &extraModels };

    // Creates a decay file based on json file input
    // We don't want this to be called by every thread!
    std::call_once( m_createDecFile_threadlock, [this]() { createDecFile(); } );

    // Read the decay file
    if ( !initialised ) {
        theGen.readUDecay( m_config.decFileName.c_str() );
        initialised = true;
    }

    // Define the histograms to be saved
    TestHistos theHistos{ m_config.testHistograms };

    // Generate events and fill histograms
    generateEvents( theGen, theHistos, firstEvent, nEvents );

    return theHistos;
}

void TestDecayModel::createDecFile() const
{
    // Create (or overwrite) the decay file
    std::ofstream decFile{ m_config.decFileName };

    // Create daughter aliases if needed
    std::vector<std::string> aliasPrefix;
    for ( std::size_t daughter_index{ 0 };
          daughter_index < m_config.daughterNames.size(); daughter_index++ ) {
        if ( !m_config.grandDaughterNames.empty() &&
             !m_config.grandDaughterNames[daughter_index].empty() ) {
            decFile << "Alias My" << m_config.daughterNames[daughter_index] << " "
                    << m_config.daughterNames[daughter_index] << std::endl;
            if ( m_config.doConjDecay[daughter_index + 1] ) {
                const EvtId daugID{
                    EvtPDL::getId( m_config.daughterNames[daughter_index] ) };
                const EvtId daugConjID{ EvtPDL::chargeConj( daugID ) };
                const std::string conjName{ daugConjID.getName() };
                std::string conjName_Alias{ daugConjID.getName() };
                if ( std::find( std::begin( m_config.daughterNames ),
                                std::end( m_config.daughterNames ),
                                daugConjID.getName() ) !=
                     std::end( m_config.daughterNames ) ) {
                    conjName_Alias = conjName_Alias + "_" + daughter_index;
                }
                decFile << "Alias My" << conjName_Alias << " " << conjName
                        << std::endl;
                decFile << "ChargeConj My"
                        << m_config.daughterNames[daughter_index] << "  My"
                        << conjName_Alias << std::endl;
            } else if ( m_config.doConjDecay[0] ) {
                decFile << "ChargeConj My"
                        << m_config.daughterNames[daughter_index] << "  My"
                        << m_config.daughterNames[daughter_index] << std::endl;
            }
            aliasPrefix.push_back( "My" );
        } else {
            aliasPrefix.push_back( "" );
        }
    }

    for ( const auto& iExtra : m_config.extraCommands ) {
        decFile << iExtra << std::endl;
    }

    // Parent decay
    decFile << "Decay " << m_config.parentName << std::endl;
    decFile << "1.0";

    for ( std::size_t daughter_index{ 0 };
          daughter_index < m_config.daughterNames.size(); daughter_index++ ) {
        decFile << " " << aliasPrefix[daughter_index]
                << m_config.daughterNames[daughter_index];
    }

    decFile << " " << m_config.modelNames[0];

    for ( const auto& par : m_config.modelParameters[0] ) {
        decFile << " " << par;
    }

    decFile << ";" << std::endl;
    decFile << "Enddecay" << std::endl;
    if ( m_config.doConjDecay[0] ) {
        EvtId parID{ EvtPDL::getId( m_config.parentName ) };
        EvtId parConjID{ EvtPDL::chargeConj( parID ) };
        decFile << "CDecay " << parConjID.getName() << std::endl;
    }

    // Daughter decays into granddaughters
    for ( std::size_t daughter_index{ 0 };
          daughter_index < m_config.grandDaughterNames.size(); daughter_index++ ) {
        if ( m_config.grandDaughterNames[daughter_index].empty() )
            continue;
        decFile << "Decay " << aliasPrefix[daughter_index]
                << m_config.daughterNames[daughter_index] << std::endl;
        decFile << "1.0";
        for ( std::size_t grandDaughter_index{ 0 };
              grandDaughter_index <
              m_config.grandDaughterNames[daughter_index].size();
              grandDaughter_index++ ) {
            decFile << " "
                    << m_config.grandDaughterNames[daughter_index][grandDaughter_index];
        }
        decFile << " " << m_config.modelNames[daughter_index + 1];
        for ( const auto& par : m_config.modelParameters[daughter_index + 1] ) {
            decFile << " " << par;
        }
        decFile << ";" << std::endl;
        decFile << "Enddecay" << std::endl;
        if ( m_config.doConjDecay[daughter_index + 1] ) {
            EvtId daugID{
                EvtPDL::getId( m_config.daughterNames[daughter_index] ) };
            EvtId daugConjID{ EvtPDL::chargeConj( daugID ) };
            std::string conjName_Alias{ daugConjID.getName() };
            if ( std::find( std::begin( m_config.daughterNames ),
                            std::end( m_config.daughterNames ),
                            daugConjID.getName() ) !=
                 std::end( m_config.daughterNames ) ) {
                conjName_Alias = conjName_Alias + "_" + daughter_index;
            }
            decFile << "CDecay " << aliasPrefix[daughter_index]
                    << conjName_Alias << std::endl;
        }
    }

    decFile << "End" << std::endl;

    decFile.close();
}

void TestDecayModel::generateEvents( EvtGen& theGen, TestHistos& theHistos,
                                     const std::size_t firstEvent,
                                     const std::size_t nEvents ) const
{
    // Generate the decays
    const EvtId parId{ EvtPDL::getId( m_config.parentName.c_str() ) };
    const EvtId conjId{ m_config.doConjDecay[0] ? EvtPDL::chargeConj( parId )
                                                : parId };

    for ( std::size_t i{ firstEvent }; i < ( firstEvent + nEvents ); i++ ) {
        if ( i % 1000 == 0 ) {
            std::cout << "Event " << firstEvent + nEvents - i << std::endl;
        }

        // seed the RNG based on the event number
        EvtRandom::setSeed( m_config.rngSeed + i );

        // Initial 4-momentum and particle
        EvtVector4R pInit( EvtPDL::getMass( parId ), 0.0, 0.0, 0.0 );
        EvtParticle* parent{ nullptr };
        if ( EvtRandom::Flat() < 0.5 ) {
            parent = EvtParticleFactory::particleFactory( parId, pInit );
        } else {
            parent = EvtParticleFactory::particleFactory( conjId, pInit );
        }

        theGen.generateDecay( parent );

        // Check for mixing (and fill histogram)
        EvtParticle* prodParent{ nullptr };
        TH1* mixedHist{ theHistos.getMixedHist() };
        if ( mixedHist ) {
            if ( parent->getNDaug() == 1 ) {
                prodParent = parent;
                parent = prodParent->getDaug( 0 );

                mixedHist->Fill( 1 );
            } else {
                mixedHist->Fill( 0 );
            }
        }

        // To debug

        if ( m_config.debugFlag ) {
            std::cout << "Parent PDG code: " << parent->getPDGId()
                      << " has daughters " << parent->getNDaug() << std::endl;
            for ( std::size_t iDaughter{ 0 }; iDaughter < parent->getNDaug();
                  iDaughter++ ) {
                std::cout << "Parent PDG code of daughter " << iDaughter
                          << " : " << parent->getDaug( iDaughter )->getPDGId()
                          << " has daughters "
                          << parent->getDaug( iDaughter )->getNDaug()
                          << std::endl;

                for ( std::size_t iGrandDaughter{ 0 };
                      iGrandDaughter < parent->getDaug( iDaughter )->getNDaug();
                      iGrandDaughter++ ) {
                    std::cout << "Parent PDG code of grand daughter "
                              << iGrandDaughter << " : "
                              << parent->getDaug( iDaughter )
                                     ->getDaug( iGrandDaughter )
                                     ->getPDGId()
                              << " has daughters "
                              << parent->getDaug( iDaughter )
                                     ->getDaug( iGrandDaughter )
                                     ->getNDaug()
                              << std::endl;
                }
            }
        }

        int leadingChargedDaughter = -1;
        const std::string fsrStr = "_FSRPhotons";
        const std::string ewStr = "_EnergyWeight";

        const std::string perDaughter = "_perDaughter";

        // Store information
        for ( auto& [info, hist] : theHistos.get1DHistos() ) {
            if ( !hist ) {
                continue;
            }

            const std::string varName = info.getName();
            std::string reducedVarName = varName;

            const std::string::size_type findFSRstr = varName.find( fsrStr );

            const std::string::size_type findPerDaugStr = varName.find(
                perDaughter );

            // If the variable name does not have the substrings '_perDaughter' and '_FSRPhotons', then add an entry per event and continue
            if ( findFSRstr == std::string::npos &&
                 findPerDaugStr == std::string::npos ) {
                const double value{
                    getValue( parent, varName, info.getd1(), info.getd2() ) };

                hist->Fill( value );
                continue;
            }
            // If the variable name has the substring '_perDaughter' and does not have the substring '_FSRPhotons', then add an entry per daughter and continue
            else if ( findFSRstr == std::string::npos ) {
                reducedVarName.erase( findPerDaugStr, perDaughter.length() );

                // Figure out whether the variable indicates to be saved for a requested particle type
                const std::string::size_type findIsStr = reducedVarName.find(
                    "_is" );
                const std::string requestedType = findIsStr == std::string::npos
                                                    ? ""
                                                    : reducedVarName.substr(
                                                          findIsStr + 3 );

                if ( !requestedType.empty() ) {
                    reducedVarName.erase( findIsStr,
                                          reducedVarName.length() - findIsStr );
                }

                for ( int iDaug{ 0 }; iDaug < (int)parent->getNDaug(); iDaug++ ) {
                    const double value{
                        getValue( parent, reducedVarName, iDaug + 1, 0 ) };

                    const int partGroup = getPartGroup(
                        parent->getDaug( iDaug )->getPDGId() );

                    if ( requestedType.empty() ||
                         isPartType( partGroup, requestedType ) ) {
                        hist->Fill( value );
                    }
                }
                continue;
            }

            // If otherwise the variable contains the substring '_FSRPhotons', then add an entry per photon
            if ( leadingChargedDaughter == -1 ) {
                leadingChargedDaughter = findChargedDaughterWithMaxE( parent );
            }

            reducedVarName.erase( findFSRstr, fsrStr.length() );

            const std::string::size_type findEwStr = reducedVarName.find( ewStr );

            // If variable name has substring '_EnergyWeight' then add an entry weighted by the photon energy
            bool energyWeight = false;
            if ( findEwStr != std::string::npos ) {
                reducedVarName.erase( findEwStr, ewStr.length() );
                energyWeight = true;
            }

            for ( int iDaug{ 0 }; iDaug < (int)parent->getNDaug(); iDaug++ ) {
                const EvtParticle* iDaughter = parent->getDaug( iDaug );

                if ( iDaughter->getAttribute( "FSR" ) == 1 ) {
                    const double value{ getValue( parent, reducedVarName,
                                                  iDaug + 1,
                                                  leadingChargedDaughter + 1 ) };
                    if ( energyWeight ) {
                        const double photonEnergy = iDaughter->getP4().get( 0 );
                        hist->Fill( value, photonEnergy );
                    } else {
                        hist->Fill( value );
                    }
                }
            }
        }

        for ( auto& [info, hist] : theHistos.get2DHistos() ) {
            if ( !hist ) {
                continue;
            }

            const double valueX{ getValue( parent,
                                           info.getName( HistInfo::Axis::X ),
                                           info.getd1( HistInfo::Axis::X ),
                                           info.getd2( HistInfo::Axis::X ) ) };
            const double valueY{ getValue( parent,
                                           info.getName( HistInfo::Axis::Y ),
                                           info.getd1( HistInfo::Axis::Y ),
                                           info.getd2( HistInfo::Axis::Y ) ) };
            hist->Fill( valueX, valueY );
        }

        if ( m_config.debugFlag ) {
            if ( prodParent ) {
                prodParent->printTree();
            } else {
                parent->printTree();
            }
        }

        // Cleanup
        if ( prodParent ) {
            prodParent->deleteTree();
        } else {
            parent->deleteTree();
        }
    }
}

int TestDecayModel::findChargedDaughterWithMaxE( const EvtParticle* parent ) const
{
    /* This function returns the index of the charged daughter with the highest energy 
     * following the sign convention below. */

    double max_E = 0;
    double max_index = 0;

    const int parentCh3 = EvtPDL::chg3( parent->getId() );

    for ( int iDaug{ 0 }; iDaug < (int)parent->getNDaug(); iDaug++ ) {
        const EvtParticle* iDaughter = parent->getDaug( iDaug );
        const int dauCh3 = EvtPDL::chg3( iDaughter->getId() );
        const double daugE = iDaughter->getP4LabBeforeFSR().get( 0 );
        // Sign convention: take negative daughter if mother is neutral,
        // otherwise the one with the same sign as the mother.
        if ( ( ( parentCh3 == 0 && dauCh3 < 0 ) || ( parentCh3 * dauCh3 > 0 ) ) &&
             daugE > max_E ) {
            max_E = daugE;
            max_index = iDaug;
        }
    }
    return max_index;
}

double TestDecayModel::getValue( const EvtParticle* parent,
                                 const std::string& varName, const int d1,
                                 const int d2 ) const
{
    double value{ std::numeric_limits<double>::quiet_NaN() };
    if ( !parent ) {
        return value;
    }

    const int NDaugMax( parent->getNDaug() );
    // If variable name contains "_daugX", we are interested in daughters of daughter X
    // Else we are interested in daughters
    const EvtParticle* selectedParent{ parent };
    std::string selectedVarName{ varName };
    if ( varName.find( "_daug" ) != std::string::npos ) {
        // Get daughter index from last character in string
        const int iDaughter{ varName.back() - '0' };
        selectedVarName = varName.substr( 0, varName.size() - 6 );
        if ( iDaughter > 0 && iDaughter <= NDaugMax ) {
            selectedParent = parent->getDaug( iDaughter - 1 );
        } else {
            return value;
        }
    }

    const int sel_NDaugMax( selectedParent->getNDaug() );
    const EvtParticle* par1{ d1 > 0 && d1 <= sel_NDaugMax
                                 ? selectedParent->getDaug( d1 - 1 )
                                 : nullptr };
    const EvtParticle* par2{ d2 > 0 && d2 <= sel_NDaugMax
                                 ? selectedParent->getDaug( d2 - 1 )
                                 : nullptr };
    const EvtParticle* par3{ sel_NDaugMax > 0
                                 ? selectedParent->getDaug( sel_NDaugMax - 1 )
                                 : nullptr };

    // 4-momenta in parent rest frame
    const EvtVector4R p1{ par1 != nullptr ? par1->getP4() : EvtVector4R() };
    const EvtVector4R p2{ par2 != nullptr ? par2->getP4() : EvtVector4R() };
    const EvtVector4R p3{ par3 != nullptr ? par3->getP4() : EvtVector4R() };

    // 4-momenta in lab frame (1st parent in decay tree)
    const EvtVector4R p1_lab{ par1 != nullptr ? par1->getP4Lab() : EvtVector4R() };
    const EvtVector4R p2_lab{ par2 != nullptr ? par2->getP4Lab() : EvtVector4R() };
    const EvtVector4R p3_lab{ par3 != nullptr ? par3->getP4Lab() : EvtVector4R() };

    if ( !selectedVarName.compare( "id" ) ) {
        // StdHep ID of one of the daughters (controlled by d1) or the parent
        if ( par1 ) {
            value = par1->getPDGId();
        } else {
            value = selectedParent->getPDGId();
        }
    } else if ( !selectedVarName.compare( "particleType" ) ) {
        if ( par1 ) {
            value = getPartGroup( par1->getPDGId() );
        } else {
            value = getPartGroup( selectedParent->getPDGId() );
        }
    } else if ( !selectedVarName.compare( "nDaug" ) ) {
        // Number of daughters
        value = sel_NDaugMax;

    } else if ( !selectedVarName.compare( "parMass" ) ) {
        // Parent invariant mass
        value = selectedParent->mass();

    } else if ( !selectedVarName.compare( "mass" ) ) {
        // Invariant mass
        if ( d2 != 0 ) {
            // Invariant 4-mass combination of particles d1 and d2
            value = ( p1 + p2 ).mass();
        } else {
            // Invariant mass of particle d1 only
            value = p1.mass();
        }

    } else if ( !selectedVarName.compare( "massSq" ) ) {
        // Invariant mass
        if ( d2 != 0 ) {
            // Invariant 4-mass combination of particles d1 and d2
            value = ( p1 + p2 ).mass2();
        } else {
            // Invariant mass of particle d1 only
            value = p1.mass2();
        }
    } else if ( !selectedVarName.compare( "mPrime" ) ) {
        // pick up first unused particle rather than last one
        if ( sel_NDaugMax != 3 ) {
            return -1;
        }
        int unused{ 0 };
        for ( int ii{ 1 }; ii <= sel_NDaugMax; ++ii ) {
            if ( ii != d1 && ii != d2 ) {
                unused = ii;
                break;
            }
        }
        if ( unused == 0 ) {
            unused = sel_NDaugMax;
        }
        const auto parL{ selectedParent->getDaug( unused - 1 ) };
        //            const auto& pL = parL->getP4();

        const double mB{ selectedParent->mass() };
        const double m1{ par1->mass() };
        const double m2{ par2->mass() };
        const double m3{ parL->mass() };
        const double m12{ ( p1 + p2 ).mass() };
        const double m12norm{
            2 * ( ( m12 - ( m1 + m2 ) ) / ( mB - ( m1 + m2 + m3 ) ) ) - 1 };
        value = acos( m12norm ) / EvtConst::pi;

    } else if ( !selectedVarName.compare( "thetaPrime" ) ) {
        // pick up first unused particle rather than last one
        if ( sel_NDaugMax != 3 ) {
            return -1;
        }
        int unused{ 0 };
        for ( int ii{ 1 }; ii <= sel_NDaugMax; ++ii ) {
            if ( ii != d1 && ii != d2 ) {
                unused = ii;
                break;
            }
        }
        if ( unused == 0 ) {
            unused = sel_NDaugMax;
        }
        const auto parL{ selectedParent->getDaug( unused - 1 ) };
        const auto& pL{ parL->getP4() };

        const double mB{ selectedParent->mass() };
        const double m1{ p1.mass() };
        const double m2{ p2.mass() };
        const double m3{ pL.mass() };
        double mBSq{ mB * mB };
        double m1Sq{ m1 * m1 };
        double m2Sq{ m2 * m2 };
        double m3Sq{ m3 * m3 };
        const double m12{ ( p1 + p2 ).mass() };
        const double m13{ ( p1 + p3 ).mass() };
        const double m12Sq{ m12 * m12 };
        const double m13Sq{ m13 * m13 };
        double en1{ ( m12Sq - m2Sq + m1Sq ) / ( 2.0 * m12 ) };
        double en3{ ( mBSq - m12Sq - m3Sq ) / ( 2.0 * m12 ) };
        double p1_12{ std::sqrt( en1 * en1 - m1Sq ) };
        double p3_12{ std::sqrt( en3 * en3 - m3Sq ) };
        double cosTheta{ ( -m13Sq + m1Sq + m3Sq + 2. * en1 * en3 ) /
                         ( 2. * p1_12 * p3_12 ) };
        value = acos( cosTheta ) / EvtConst::pi;

    } else if ( !selectedVarName.compare( "pSumSq" ) ) {
        // Invariant momentum sum squared
        value = ( p1 + p2 ).mass2();

    } else if ( !selectedVarName.compare( "pDiffSq" ) ) {
        // Invariant momentum difference squared
        value = ( p1 - p2 ).mass2();

    } else if ( !selectedVarName.compare( "mass3" ) ) {
        // Invariant mass of 3 daughters
        value = ( p1 + p2 + p3 ).mass();

    } else if ( !selectedVarName.compare( "mass3_specified" ) ) {
        // Invariant mass of 3 daughters, d1 is first daughter
        // second daughter is d1 + 1, d2 is last daughter
        const EvtParticle* daug2{ selectedParent->getDaug( d1 ) };
        const EvtParticle* daug3{ selectedParent->getDaug( d2 - 1 ) };

        const EvtVector4R p2_specified{ daug2 != nullptr ? daug2->getP4Lab()
                                                         : EvtVector4R() };
        const EvtVector4R p3_specified{ daug3 != nullptr ? daug3->getP4Lab()
                                                         : EvtVector4R() };

        value = ( p1 + p2_specified + p3_specified ).mass();

    } else if ( !selectedVarName.compare( "cosTheta3" ) ) {
        // Cosine of the polar angle of the momentum of d1 + d2 + d3
        const EvtVector4R p123{ p1 + p2 + p3 };
        if ( p123.d3mag() > 0.0 ) {
            value = p123.get( 3 ) / p123.d3mag();
        }

    } else if ( !selectedVarName.compare( "pLab" ) ) {
        // Momentum of particle d1 in lab frame
        value = p1_lab.d3mag();

    } else if ( !selectedVarName.compare( "p" ) ) {
        // Momentum of particle d1 in parent rest frame
        value = p1.d3mag();

    } else if ( !selectedVarName.compare( "pLabSq" ) ) {
        // Momentum squared of particle d1 (in lab frame)
        const double p1_lab_x{ p1_lab.get( 1 ) };
        const double p1_lab_y{ p1_lab.get( 2 ) };
        const double p1_lab_z{ p1_lab.get( 3 ) };
        value = p1_lab_x * p1_lab_x + p1_lab_y * p1_lab_y + p1_lab_z * p1_lab_z;

    } else if ( !selectedVarName.compare( "pSq" ) ) {
        // Momentum squared of particle d1 (in lab frame)
        const double p1_x{ p1.get( 1 ) };
        const double p1_y{ p1.get( 2 ) };
        const double p1_z{ p1.get( 3 ) };
        value = p1_x * p1_x + p1_y * p1_y + p1_z * p1_z;

    } else if ( !selectedVarName.compare( "pxLab" ) ) {
        // x momentum of particle d1 in lab frame
        value = p1_lab.get( 1 );

    } else if ( !selectedVarName.compare( "px" ) ) {
        // x momentum of particle d1 in parent frame
        value = p1.get( 1 );

    } else if ( !selectedVarName.compare( "pyLab" ) ) {
        // y momentum of particle d1 in lab frame
        value = p1_lab.get( 2 );

    } else if ( !selectedVarName.compare( "py" ) ) {
        // y momentum of particle d1 in parent frame
        value = p1.get( 2 );

    } else if ( !selectedVarName.compare( "pzLab" ) ) {
        // z momentum of particle d1 in lab frame
        value = p1_lab.get( 3 );

    } else if ( !selectedVarName.compare( "pz" ) ) {
        // z momentum of particle d1 in parent frame
        value = p1.get( 3 );

    } else if ( !selectedVarName.compare( "cosHel" ) ||
                !selectedVarName.compare( "absCosHel" ) ||
                !selectedVarName.compare( "cosHelParent" ) ) {
        // Cosine of helicity angle
        EvtVector4R p12;
        EvtVector4R p1Res;

        if ( !selectedVarName.compare( "cosHelParent" ) ) {
            p12 = selectedParent->getP4Lab();
            const EvtVector4R boost{ p12.get( 0 ), -p12.get( 1 ), -p12.get( 2 ),
                                     -p12.get( 3 ) };
            // Momentum of particle d1 in resonance frame, p1Res
            p1Res = boostTo( p1_lab, boost );
        } else if ( d2 != 0 ) {
            // Resonance center-of-mass system (d1 and d2)
            p12 = p1_lab + p2_lab;
            // Boost vector
            const EvtVector4R boost{ p12.get( 0 ), -p12.get( 1 ), -p12.get( 2 ),
                                     -p12.get( 3 ) };
            // Momentum of particle d1 in resonance frame, p1Res
            p1Res = boostTo( p1_lab, boost );
        } else {
            // The resonance is d1
            p12 = p1;

            // Find its first daughter
            const EvtParticle* gpar{ par1 != nullptr ? par1->getDaug( 0 )
                                                     : nullptr };

            p1Res = gpar != nullptr ? gpar->getP4() : EvtVector4R();
        }

        // Cosine of angle between p1Res and momentum of resonance in parent frame
        const double p1ResMag{ p1Res.d3mag() };
        const double p12Mag{ p12.d3mag() };
        if ( p1ResMag > 0.0 && p12Mag > 0.0 ) {
            value = -p1Res.dot( p12 ) / ( p1ResMag * p12Mag );
        }

        if ( !selectedVarName.compare( "absCosHel" ) ) {
            value = std::abs( value );
        }

    } else if ( !selectedVarName.compare( "cosHelTau" ) ) {
        // Works only for B -> X nu_1 tau -> pi nu_2.
        // Cosine of helicity angle between pi momentum and opposite W momentum -(nu_1 + tau)
        // in the tau rest frame. Index d1 must match with tau.
        // p3 (momentum of last daughter) is taken as the neutrino momentum

        // W momentum
        const EvtVector4R p_W{ p1_lab + p3_lab };

        // Index d2 must match the index of the pion (daughter of tau)
        const EvtParticle* pion{
            selectedParent->getDaug( d1 - 1 )->getDaug( d2 - 1 ) };
        const EvtVector4R p_pion{ pion != nullptr ? pion->getP4Lab()
                                                  : EvtVector4R() };

        // Boost vector to tau frame
        const EvtVector4R boost{ p1_lab.get( 0 ), -p1_lab.get( 1 ),
                                 -p1_lab.get( 2 ), -p1_lab.get( 3 ) };

        // Boost both momenta to tau frame
        const EvtVector4R p_W_boosted{ boostTo( p_W, boost ) };
        const EvtVector4R p_pion_boosted{ boostTo( p_pion, boost ) };
        // Cosine of angle between opposite W momentum and pion momentum in tau frame
        const double p_W_boostedMag{ p_W_boosted.d3mag() };
        const double p_pion_boostedMag{ p_pion_boosted.d3mag() };

        if ( p_W_boostedMag > 0.0 && p_pion_boostedMag > 0.0 ) {
            value = -p_W_boosted.dot( p_pion_boosted ) /
                    ( p_W_boostedMag * p_pion_boostedMag );
        }

    } else if ( !selectedVarName.compare( "cosHelDiTau" ) ||
                !selectedVarName.compare( "cosHelDiTau_over05" ) ||
                !selectedVarName.compare( "cosHelDiTau_under05" ) ) {
        // Works for B -> tau (pi nu) tau -> (pi nu), B -> phi (-> KK) l l decays,
        // B -> 4 pions, or similar decays..
        // Cosine of helicity angle between pi momentum and opposite B momentum in tau rest frame.
        // Index d1 must match with tau

        if ( sel_NDaugMax < 2 || sel_NDaugMax > 4 ) {
            return value;
        }

        // B momentum
        const EvtVector4R p_B{ selectedParent->getP4Lab() };

        // Index d2 must match the index of the pion (daughter of tau)
        const EvtParticle* pion_1{
            sel_NDaugMax <= 3
                ? selectedParent->getDaug( d1 - 1 )->getDaug( d2 - 1 )
                : selectedParent->getDaug( d1 - 1 ) };
        const EvtVector4R p_pion_1{ pion_1 != nullptr ? pion_1->getP4Lab()
                                                      : EvtVector4R() };

        const EvtVector4R p_first_daughter{
            sel_NDaugMax <= 3 ? selectedParent->getDaug( d1 - 1 )->getP4Lab()
                              : selectedParent->getDaug( d1 - 1 )->getP4Lab() +
                                    selectedParent->getDaug( d1 )->getP4Lab() };

        // Boost vector to tau frame
        const EvtVector4R boost_1{ p_first_daughter.get( 0 ),
                                   -p_first_daughter.get( 1 ),
                                   -p_first_daughter.get( 2 ),
                                   -p_first_daughter.get( 3 ) };

        // Boost both momenta to tau frame
        const EvtVector4R p_B_boosted_1{ boostTo( p_B, boost_1 ) };
        const EvtVector4R p_pion_boosted_1{ boostTo( p_pion_1, boost_1 ) };
        // Cosine of angle between opposite W momentum and pion momentum in tau frame
        const double p_B_boosted_1_Mag{ p_B_boosted_1.d3mag() };
        const double p_pion_boosted_1_Mag{ p_pion_boosted_1.d3mag() };

        double hel1{ std::numeric_limits<double>::quiet_NaN() };
        double hel{ std::numeric_limits<double>::quiet_NaN() };

        if ( p_B_boosted_1_Mag > 0.0 && p_pion_boosted_1_Mag > 0.0 ) {
            hel1 = -p_B_boosted_1.dot( p_pion_boosted_1 ) /
                   ( p_B_boosted_1_Mag * p_pion_boosted_1_Mag );
        }

        if ( ( !selectedVarName.compare( "cosHelDiTau_over05" ) ) ||
             ( !selectedVarName.compare( "cosHelDiTau_under05" ) ) ) {
            // Works for B -> tau (pi nu) tau -> (pi nu) or similar decays.
            // Cosine of helicity angle between pi momentum and opposite B momentum in tau rest frame
            // Index d1 must match with tau; cosHelicity above +0.5 or below -0.5

            // Index d2 must match the index of the pion (daughter of tau)
            const EvtParticle* pion{
                sel_NDaugMax <= 3 ? selectedParent->getDaug( 0 )->getDaug( d2 - 1 )
                                  : selectedParent->getDaug( 0 ) };
            const EvtVector4R p_pion{ pion != nullptr ? pion->getP4Lab()
                                                      : EvtVector4R() };
            // Boost vector to tau frame
            const EvtVector4R p_second_daughter{
                sel_NDaugMax == 2
                    ? selectedParent->getDaug( 0 )->getP4Lab()
                    : selectedParent->getDaug( 0 )->getP4Lab() +
                          selectedParent->getDaug( 1 )->getP4Lab() };

            const EvtVector4R boost{ p_second_daughter.get( 0 ),
                                     -p_second_daughter.get( 1 ),
                                     -p_second_daughter.get( 2 ),
                                     -p_second_daughter.get( 3 ) };

            // Boost both momenta to tau frame
            const EvtVector4R p_B_boosted{ boostTo( p_B, boost ) };
            const EvtVector4R p_pion_boosted{ boostTo( p_pion, boost ) };
            // Cosine of angle between opposite W momentum and pion momentum in tau frame
            const double p_B_boostedMag{ p_B_boosted.d3mag() };
            const double p_pion_boostedMag{ p_pion_boosted.d3mag() };

            if ( p_B_boostedMag > 0.0 && p_pion_boostedMag > 0.0 ) {
                hel = -p_B_boosted.dot( p_pion_boosted ) /
                      ( p_B_boostedMag * p_pion_boostedMag );
            }
        }

        if ( !selectedVarName.compare( "cosHelDiTau" ) ||
             ( hel > 0.5 && !selectedVarName.compare( "cosHelDiTau_over05" ) ) ||
             ( hel < -0.5 && !selectedVarName.compare( "cosHelDiTau_under05" ) ) ) {
            value = hel1;
        }

    } else if ( !selectedVarName.compare( "cosAcoplanarityAngle" ) ||
                !selectedVarName.compare( "acoplanarityAngle" ) ) {
        // Acoplanarity angle or cosine for B -> tau (pi nu) tau -> (pi nu),
        // B -> phi (-> KK) l l decays, B -> 4 pions, or similar decays

        value = getCosAcoplanarityAngle( selectedParent, sel_NDaugMax, d1, d2 );
        if ( !selectedVarName.compare( "acoplanarityAngle" ) ) {
            value = std::acos( value ) * 180.0 / EvtConst::pi;    // degrees
        }

    } else if ( !selectedVarName.compare( "cosTheta" ) ) {
        // Cosine of polar angle of first daughter in lab frame
        const double p1_lab_mag{ p1_lab.d3mag() };
        if ( p1_lab_mag > 0.0 ) {
            value = p1_lab.get( 3 ) / p1_lab_mag;
        }

    } else if ( !selectedVarName.compare( "phi" ) ) {
        // Azimuthal angle of first daughter in lab frame (degrees)
        const double p1_lab_mag{ p1_lab.d3mag() };
        if ( p1_lab_mag > 0.0 ) {
            value = atan2( p1_lab.get( 1 ), p1_lab.get( 2 ) ) * 180.0 /
                    EvtConst::pi;
        }

    } else if ( !selectedVarName.compare( "openingAngle" ) ) {
        // Polar angle between first and second daughters in parent's frame

        const double cost{ p1.dot( p2 ) / ( p1.d3mag() * p2.d3mag() ) };

        value = acos( cost ) * 180.0 / EvtConst::pi;

    } else if ( !selectedVarName.compare( "decayangle" ) ) {
        // Polar angle between first and second daughters in lab frame

        const EvtVector4R p{ selectedParent->getP4() };
        const EvtVector4R q{ p1 + p2 };
        const EvtVector4R d{ p1 };

        const double cost{ EvtDecayAngle( p, q, d ) };

        value = acos( cost ) * 180.0 / EvtConst::pi;

    } else if ( !selectedVarName.compare( "decayangle3" ) ) {
        // Polar angle between combined first and second daughters
        // with the third daughter in lab frame

        const EvtVector4R p{ selectedParent->getP4() };
        const EvtVector4R q{ p1 + p2 + p3 };
        const EvtVector4R d{ p1 + p2 };

        const double cost{ EvtDecayAngle( p, q, d ) };

        value = acos( cost ) * 180.0 / EvtConst::pi;

    } else if ( !selectedVarName.compare( "decayangle_BTO4PI" ) ) {
        // Polar angle between first and second daughters in lab frame.
        // Used in PIPIPI (BTO4PI_CP) model

        const EvtVector4R p{ p1 + p2 + p3 };
        const EvtVector4R q{ p1 + p2 };
        const EvtVector4R d{ p1 };

        const double cost{ EvtDecayAngle( p, q, d ) };

        value = acos( cost ) * 180.0 / EvtConst::pi;

    } else if ( !selectedVarName.compare( "chi" ) ) {
        // Chi angle (degrees) using all 4 daughters and parent 4 mom

        const EvtParticle* daug1{ selectedParent->getDaug( d1 - 1 ) };
        const EvtParticle* daug2{ selectedParent->getDaug( d1 ) };
        const EvtParticle* daug3{ selectedParent->getDaug( d1 + 1 ) };
        const EvtParticle* daug4{ selectedParent->getDaug( d1 + 2 ) };

        const EvtVector4R p_parent{ selectedParent->getP4() };
        const EvtVector4R p_daug1{ daug1 != nullptr ? daug1->getP4()
                                                    : EvtVector4R() };
        const EvtVector4R p_daug2{ daug2 != nullptr ? daug2->getP4()
                                                    : EvtVector4R() };
        const EvtVector4R p_daug3{ daug3 != nullptr ? daug3->getP4()
                                                    : EvtVector4R() };
        const EvtVector4R p_daug4{ daug4 != nullptr ? daug4->getP4()
                                                    : EvtVector4R() };

        const double chi{
            EvtDecayAngleChi( p_parent, p_daug1, p_daug2, p_daug3, p_daug4 ) };
        value = chi * 180.0 / EvtConst::pi;

    } else if ( !selectedVarName.compare( "cosThetaResNorm" ) ) {
        // P -> R p4 -> (p1 p2 p3) p4, where the resonance R decays to p1 p2 p3.
        // Theta is the angle between the normal of the plane containing p1, p2 & p3
        // and the bachelor particle p4 in the rest frame of resonance R.
        // The normal vector is given by the cross product p3 x p1

        if ( sel_NDaugMax > 1 ) {
            const EvtParticle* res{ selectedParent->getDaug( 0 ) };
            const EvtParticle* bac{ selectedParent->getDaug( 1 ) };

            // Check resonance has 3 daughters
            if ( res != nullptr && res->getNDaug() == 3 ) {
                const EvtParticle* daug1 = res->getDaug( 0 );
                const EvtParticle* daug3 = res->getDaug( 2 );

                // 4-momenta in base parent P lab frame
                const EvtVector4R p4_Res{ res->getP4Lab() };
                const EvtVector4R p4_p4{ bac != nullptr ? bac->getP4Lab()
                                                        : EvtVector4R() };
                const EvtVector4R p4_p1{ daug1 != nullptr ? daug1->getP4Lab()
                                                          : EvtVector4R() };
                const EvtVector4R p4_p3{ daug3 != nullptr ? daug3->getP4Lab()
                                                          : EvtVector4R() };

                // Boost 4-vector for resonance frame
                const EvtVector4R boost{ p4_Res.get( 0 ), -p4_Res.get( 1 ),
                                         -p4_Res.get( 2 ), -p4_Res.get( 3 ) };

                // Momentum of p1 and p3 in resonance frame
                const EvtVector4R p1Res{ boostTo( p4_p1, boost ) };
                const EvtVector4R p3Res{ boostTo( p4_p3, boost ) };

                // Plane normal vector (just uses 3-momentum components)
                const EvtVector4R norm{ p3Res.cross( p1Res ) };

                // Momentum of p4 in resonance frame
                const EvtVector4R p4Res{ boostTo( p4_p4, boost ) };

                // Cosine of the angle between the normal and p4 in the resonance frame
                const double normMag{ norm.d3mag() };
                const double p4ResMag{ p4Res.d3mag() };
                if ( normMag > 0.0 && p4ResMag > 0.0 ) {
                    value = norm.dot( p4Res ) / ( normMag * p4ResMag );
                }
            }
        }

    } else if ( !selectedVarName.compare( "cosBetaRes" ) ) {
        // For resonance P (parent) -> p1 p2 p3, beta is the
        // angle between p1 & p3 in the (p1 + p2) rest frame
        if ( sel_NDaugMax > 2 ) {
            const EvtParticle* daug1 = selectedParent->getDaug( 0 );
            const EvtParticle* daug2 = selectedParent->getDaug( 1 );
            const EvtParticle* daug3 = selectedParent->getDaug( 2 );

            // 4-momenta in base parent frame
            const EvtVector4R p4_p1{ daug1 != nullptr ? daug1->getP4Lab()
                                                      : EvtVector4R() };
            const EvtVector4R p4_p2{ daug2 != nullptr ? daug2->getP4Lab()
                                                      : EvtVector4R() };
            const EvtVector4R p4_p3{ daug3 != nullptr ? daug3->getP4Lab()
                                                      : EvtVector4R() };

            // p1 + p2
            const EvtVector4R p12{ p4_p1 + p4_p2 };

            // Boost 4-vector for p12 frame
            const EvtVector4R boost{ p12.get( 0 ), -p12.get( 1 ), -p12.get( 2 ),
                                     -p12.get( 3 ) };

            // Momentum of p1 & p3 in p12 frame
            const EvtVector4R p1_12{ boostTo( p4_p1, boost ) };
            const EvtVector4R p3_12{ boostTo( p4_p3, boost ) };

            // Cosine of angle between p1 & p3 in p12 frame
            const double p1_12Mag{ p1_12.d3mag() };
            const double p3_12Mag{ p3_12.d3mag() };
            if ( p1_12Mag > 0.0 && p3_12Mag > 0.0 ) {
                value = p1_12.dot( p3_12 ) / ( p1_12Mag * p3_12Mag );
            }
        }

    } else if ( !selectedVarName.compare( "E" ) ) {
        // Energy of first daughter in lab frame
        value = p1_lab.get( 0 );

    } else if ( !selectedVarName.compare( "E_over_Eparent" ) ) {
        // Energy of first daughter w.r.t parent energy (lab frame)
        value = p1_lab.get( 0 ) / selectedParent->getP4Lab().get( 0 );

    } else if ( !selectedVarName.compare( "E_over_Eparent_over05" ) ) {
        // First daughter E_over_Eparent (lab frame) if d2 granddaughter E ratio > 0.5
        const double E_over_Eparent_2{
            parent->getDaug( 0 )->getDaug( d2 - 1 )->getP4Lab().get( 0 ) /
            parent->getDaug( 0 )->getP4Lab().get( 0 ) };

        if ( E_over_Eparent_2 > 0.5 ) {
            value = p1_lab.get( 0 ) / selectedParent->getP4Lab().get( 0 );
        }

    } else if ( !selectedVarName.compare( "totE_over_Mparent" ) ) {
        // Energy of given daughters w.r.t parent mass (lab frame)
        value = ( p1_lab.get( 0 ) + p2_lab.get( 0 ) ) / selectedParent->mass();

    } else if ( !selectedVarName.compare( "prob" ) ) {
        // Decay probability
        const double* dProb{ selectedParent->decayProb() };
        if ( dProb ) {
            value = *dProb;
        }

    } else if ( !selectedVarName.compare( "lifetime" ) ) {
        // Lifetime of particle d1 (ps)
        if ( d1 ) {
            value = par1 != nullptr ? par1->getLifetime() * 1e12 / EvtConst::c
                                    : 0.0;
        } else {
            value = parent != nullptr
                      ? parent->getLifetime() * 1e12 / EvtConst::c
                      : 0.0;
        }

    } else if ( !selectedVarName.compare( "deltaT" ) ) {
        // Lifetime difference between particles d1 and d2
        const double t1{
            par1 != nullptr ? par1->getLifetime() * 1e12 / EvtConst::c : 0.0 };
        const double t2{
            par2 != nullptr ? par2->getLifetime() * 1e12 / EvtConst::c : 0.0 };
        value = t1 - t2;

    } else if ( !selectedVarName.compare( "decTime" ) ) {
        // Decay flight time of particle d1 in picoseconds.
        // Decay vertex = position of (1st) decay particle of d1
        const EvtParticle* gpar{ par1 != nullptr ? par1->getDaug( 0 ) : nullptr };
        const EvtVector4R vtxPos{ gpar != nullptr ? gpar->get4Pos()
                                                  : EvtVector4R() };
        const double p{ p1_lab.d3mag() };
        value = p > 0.0
                  ? 1e12 * vtxPos.d3mag() * p1_lab.mass() / ( p * EvtConst::c )
                  : 0.0;
    } else if ( !selectedVarName.compare( "nFSRPhotons" ) ) {
        // Loop over all daughters and get number of FSR photons

        double nFSRPhotons{ 0 };

        for ( std::size_t iDaughter{ 0 };
              iDaughter < selectedParent->getNDaug(); iDaughter++ ) {
            const EvtParticle* iDaug = selectedParent->getDaug( iDaughter );

            if ( iDaug->getAttribute( "FSR" ) == 1 )
                nFSRPhotons += 1.0;
        }

        value = nFSRPhotons;

    } else if ( !selectedVarName.compare( "totalFSREnergy" ) ) {
        // Loop over all daughters and get number of FSR photons

        double totalFSREnergy{ 0 };

        for ( std::size_t iDaughter{ 0 };
              iDaughter < selectedParent->getNDaug(); iDaughter++ ) {
            const EvtParticle* iDaug = selectedParent->getDaug( iDaughter );

            if ( iDaug->getAttribute( "FSR" ) == 1 )
                totalFSREnergy += iDaug->getP4Lab().get( 0 );
        }

        value = totalFSREnergy;
    } else {
        std::cerr << "Warning: Did not recognise variable name "
                  << selectedVarName << std::endl;
    }

    return value;
}

void TestDecayModel::compareHistos( const TestHistos& theHistos,
                                    const std::string& refFileName ) const
{
    // Compare histograms with the same name, calculating the chi-squared
    std::unique_ptr<TFile> refFile{ TFile::Open( refFileName.c_str(), "read" ) };

    if ( !refFile ) {
        std::cerr << "Could not open reference file " << refFileName << std::endl;
        return;
    }

    // TODO - should we plot the (signed) chisq histogram? and save it as pdf/png?

    // TODO - add comparison of mixedHist

    for ( auto& [_, hist] : theHistos.get1DHistos() ) {
        const std::string histName{ hist->GetName() };
        // Get equivalent reference histogram
        const TH1* refHist{
            dynamic_cast<TH1*>( refFile->Get( histName.c_str() ) ) };

        if ( refHist ) {
            double chiSq{ 0.0 };
            int nDof{ 0 };
            int iGood{ 0 };
            const double pValue{
                refHist->Chi2TestX( hist.get(), chiSq, nDof, iGood, "WW" ) };
            const double integral{ refHist->Integral() };
            std::cout << "Histogram " << histName << " chiSq/nDof = " << chiSq
                      << "/" << nDof << ", pValue = " << pValue
                      << ", integral = " << integral << std::endl;

        } else {
            std::cerr << "Could not find reference histogram " << histName
                      << std::endl;
        }
    }

    for ( auto& [_, hist] : theHistos.get2DHistos() ) {
        const std::string histName{ hist->GetName() };
        // Get equivalent reference histogram
        const TH2* refHist{
            dynamic_cast<TH2*>( refFile->Get( histName.c_str() ) ) };

        if ( refHist ) {
            double chiSq{ 0.0 };
            int nDof{ 0 };
            int iGood{ 0 };
            const double pValue{
                refHist->Chi2TestX( hist.get(), chiSq, nDof, iGood, "WW" ) };
            const double integral{ refHist->Integral() };
            std::cout << "Histogram " << histName << " chiSq/nDof = " << chiSq
                      << "/" << nDof << ", pValue = " << pValue
                      << ", integral = " << integral << std::endl;

        } else {
            std::cerr << "Could not find reference histogram " << histName
                      << std::endl;
        }
    }

    refFile->Close();
}

int TestDecayModel::getPartGroup( const int PDGId ) const
{
    int group( -1 );

    const int absPDGId = std::abs( PDGId );

    if ( absPDGId >= 11 && absPDGId <= 16 ) {
        group = 0;    // leptons
    } else if ( absPDGId == 22 ) {
        group = 1;    // photon
    } else if ( absPDGId == 211 ) {
        group = 2;    // pi+-
    } else if ( absPDGId == 111 ) {
        group = 3;    // pi0
    } else if ( absPDGId == 321 ) {
        group = 4;    // K+-
    } else if ( absPDGId == 311 || absPDGId == 130 || absPDGId == 310 ) {
        group = 5;    // K0
    } else if ( absPDGId == 411 ) {
        group = 6;    // D+-
    } else if ( absPDGId == 421 ) {
        group = 7;    // D0
    } else if ( absPDGId == 2212 || absPDGId == 2112 || absPDGId == 2224 ||
                absPDGId == 2214 || absPDGId == 2114 || absPDGId == 1114 ) {
        group = 8;    // light baryons
    } else if ( absPDGId >= 3112 && absPDGId <= 3334 ) {
        group = 9;    // strange baryons
    } else if ( absPDGId != 0 ) {
        group = 10;    // other particles
    }

    return group;
}

bool TestDecayModel::isPartType( const int group,
                                 const std::string& particleType ) const
{
    if ( ( !particleType.compare( "LeptonOrPhoton" ) &&
           ( group == 0 || group == 1 ) ) ||
         ( !particleType.compare( "Pion" ) && ( group == 2 || group == 3 ) ) ||
         ( !particleType.compare( "Kaon" ) && ( group == 4 || group == 5 ) ) ||
         ( !particleType.compare( "Dmeson" ) && ( group == 6 || group == 7 ) ) ||
         ( !particleType.compare( "Baryon" ) && ( group == 8 || group == 9 ) ) ||
         ( !particleType.compare( "Other" ) && group == 10 ) ) {
        return true;
    } else {
        return false;
    }
}

double TestDecayModel::getCosAcoplanarityAngle( const EvtParticle* selectedParent,
                                                const int sel_NDaugMax,
                                                const int d1, const int d2 ) const
{
    // Given a two-body decay, the acoplanarity angle is defined as the angle between the
    // two decay planes (normal vectors) in the reference frame of the mother.
    // Each normal vector is the cross product of the momentum of one daughter (in the frame of the mother)
    // and the momentum of one of the granddaughters (in the reference frame of the daughter).
    // In case of 3 daughters, we build an intermediate daughter (2+3) from the last 2 daughters
    // (which are then treated as grand daughters), and in case of 4 daughters, we build
    // 2 intermediate daughters (1+2) and (3+4)

    double cosAco{ 0.0 };

    if ( sel_NDaugMax < 2 || sel_NDaugMax > 4 ) {
        return cosAco;
    }

    const EvtParticle* daughter1{ selectedParent->getDaug( 0 ) };
    const EvtParticle* daughter2{ selectedParent->getDaug( 1 ) };
    const EvtParticle* daughter3{ selectedParent->getDaug( 2 ) };
    const EvtParticle* daughter4{ selectedParent->getDaug( 3 ) };

    const EvtParticle* grandDaughter1{ daughter1->getDaug( d1 - 1 ) };
    const EvtParticle* grandDaughter2{ daughter2->getDaug( d2 - 1 ) };

    const EvtVector4R parent4Vector{ selectedParent->getP4Lab() };

    const EvtVector4R daughter4Vector1{
        sel_NDaugMax <= 3 ? daughter1->getP4Lab()
                          : daughter1->getP4Lab() + daughter2->getP4Lab() };

    const EvtVector4R daughter4Vector2{
        sel_NDaugMax == 2   ? daughter2->getP4Lab()
        : sel_NDaugMax == 3 ? daughter2->getP4Lab() + daughter3->getP4Lab()
                            : daughter3->getP4Lab() + daughter4->getP4Lab() };

    const EvtVector4R grandDaughter4Vector1{
        sel_NDaugMax <= 3 ? grandDaughter1->getP4Lab() : daughter1->getP4Lab() };

    const EvtVector4R grandDaughter4Vector2{
        sel_NDaugMax == 2   ? grandDaughter2->getP4Lab()
        : sel_NDaugMax == 3 ? daughter2->getP4Lab()
                            : daughter3->getP4Lab() };

    const EvtVector4R parentBoost{ parent4Vector.get( 0 ),
                                   -parent4Vector.get( 1 ),
                                   -parent4Vector.get( 2 ),
                                   -parent4Vector.get( 3 ) };

    const EvtVector4R daughter1Boost{ daughter4Vector1.get( 0 ),
                                      -daughter4Vector1.get( 1 ),
                                      -daughter4Vector1.get( 2 ),
                                      -daughter4Vector1.get( 3 ) };

    const EvtVector4R daughter2Boost{ daughter4Vector2.get( 0 ),
                                      -daughter4Vector2.get( 1 ),
                                      -daughter4Vector2.get( 2 ),
                                      -daughter4Vector2.get( 3 ) };

    // Boosting daughters to reference frame of the mother
    const EvtVector4R daughter4Vector1_boosted{
        boostTo( daughter4Vector1, parentBoost ) };
    const EvtVector4R daughter4Vector2_boosted{
        boostTo( daughter4Vector2, parentBoost ) };

    // Boosting each granddaughter to reference frame of its mother
    const EvtVector4R grandDaughter4Vector1_boosted{
        boostTo( grandDaughter4Vector1, daughter1Boost ) };
    const EvtVector4R grandDaughter4Vector2_boosted{
        boostTo( grandDaughter4Vector2, daughter2Boost ) };

    // We calculate the normal vectors of the decay two planes
    const EvtVector4R normalVector1{
        daughter4Vector1_boosted.cross( grandDaughter4Vector1_boosted ) };
    const EvtVector4R normalVector2{
        daughter4Vector2_boosted.cross( grandDaughter4Vector2_boosted ) };

    const double normalVector1Mag{ normalVector1.d3mag() };
    const double normalVector2Mag{ normalVector2.d3mag() };

    if ( normalVector1Mag > 0.0 && normalVector2Mag > 0.0 ) {
        cosAco = normalVector1.dot( normalVector2 ) /
                 ( normalVector1Mag * normalVector2Mag );
    }

    return cosAco;
}

int main( int argc, char* argv[] )
{
    if ( argc < 2 || argc > 3 ) {
        std::cerr << "Expecting at least one argument: test configuration json file"
                  << "\nAdditional argument supported for: general configuration json file"
                  << std::endl;
        return 1;
    }

    // Tweak ROOT behaviour
    ROOT::EnableThreadSafety();
    TH1::AddDirectory( kFALSE );

    // Load input file in json format
    const std::string testConfigFileName{ argv[1] };
    const std::string generalConfigFileName{
        ( argc > 2 ) ? argv[2] : "jsonFiles/config/default.json" };

    json generalConfig;
    json testConfig;

    {
        std::ifstream inputStr{ generalConfigFileName };
        inputStr >> generalConfig;
    }

    {
        std::ifstream inputStr{ testConfigFileName };
        inputStr >> testConfig;
    }

    if ( testConfig.is_array() ) {
        for ( auto& cc : testConfig ) {
            cc.merge_patch( generalConfig );
            TestDecayModel test{ cc };
            test.run();
        }
    } else {
        testConfig.merge_patch( generalConfig );
        TestDecayModel test{ testConfig };
        test.run();
    }

    return 0;
}
