
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

#include <array>
#include <cmath>
#include <mutex>
#include <string>
#include <utility>
#include <vector>

class EvtGen;
class EvtParticle;

/**
   \enum FSRGenerator
   Enumeration of possible FSR generators and their settings
*/
enum class FSRGenerator
{
    PHOTOS,            ///< PHOTOS generator
    SherpaPhotons1,    ///< Sherpa Photons generator, mode 1
    SherpaPhotons20,    ///< Sherpa Photons generator, mode 2, exact ME calculations off
    SherpaPhotons21    ///< Sherpa Photons generator, mode 2, exact ME calculations on
};

/**
   \enum ThreadModel
   Enumeration of possible models for multi-threading
*/
enum class ThreadModel
{
    StdLib,    ///< C++ standard library model
    TBB        ///< Intel Thread Building Blocks model
};

//! \cond DOXYGEN_IGNORE

/*
   Here we define a custom version of the NLOHMANN_JSON_SERIALIZE_ENUM macro.
   The macro supplied by the nlohmann_json package will convert any
   unrecognised string to the first enum state, which here we throw an
   exception in such cases.
   In addition, we define the to_string utility function.
*/

#define CUSTOM_NLOHMANN_JSON_SERIALIZE_ENUM( ENUM_TYPE, ... )                     \
    template <typename BasicJsonType>                                             \
    inline void to_json( BasicJsonType& j, const ENUM_TYPE& e )                   \
    {                                                                             \
        static_assert( std::is_enum<ENUM_TYPE>::value,                            \
                       #ENUM_TYPE " must be an enum!" );                          \
        static const std::pair<ENUM_TYPE, BasicJsonType> m[] = __VA_ARGS__;       \
        auto it = std::find_if(                                                   \
            std::begin( m ), std::end( m ),                                       \
            [e]( const std::pair<ENUM_TYPE, BasicJsonType>& ej_pair ) -> bool {   \
                return ej_pair.first == e;                                        \
            } );                                                                  \
        if ( it == std::end( m ) ) {                                              \
            throw std::runtime_error{                                             \
                "ERROR in conversion to JSON, invalid state for " #ENUM_TYPE };   \
        }                                                                         \
        j = it->second;                                                           \
    }                                                                             \
    template <typename BasicJsonType>                                             \
    inline void from_json( const BasicJsonType& j, ENUM_TYPE& e )                 \
    {                                                                             \
        static_assert( std::is_enum<ENUM_TYPE>::value,                            \
                       #ENUM_TYPE " must be an enum!" );                          \
        static const std::pair<ENUM_TYPE, BasicJsonType> m[] = __VA_ARGS__;       \
        auto it = std::find_if(                                                   \
            std::begin( m ), std::end( m ),                                       \
            [&j]( const std::pair<ENUM_TYPE, BasicJsonType>& ej_pair ) -> bool {  \
                return ej_pair.second == j;                                       \
            } );                                                                  \
        if ( it == std::end( m ) ) {                                              \
            throw std::runtime_error{                                             \
                "ERROR in conversion from JSON, invalid state for " #ENUM_TYPE }; \
        }                                                                         \
        e = it->first;                                                            \
    }                                                                             \
    inline std::string to_string( const ENUM_TYPE& e )                            \
    {                                                                             \
        static_assert( std::is_enum<ENUM_TYPE>::value,                            \
                       #ENUM_TYPE " must be an enum!" );                          \
        nlohmann::json j;                                                         \
        to_json( j, e );                                                          \
        return j.get<std::string>();                                              \
    }

// Generate functions for serialisation of the FSRGenerator and ThreadModel enums

CUSTOM_NLOHMANN_JSON_SERIALIZE_ENUM(
    FSRGenerator, {
                      { FSRGenerator::PHOTOS, "PHOTOS" },
                      { FSRGenerator::SherpaPhotons1, "SherpaPhotons1" },
                      { FSRGenerator::SherpaPhotons20, "SherpaPhotons20" },
                      { FSRGenerator::SherpaPhotons21, "SherpaPhotons21" },
                  } )

CUSTOM_NLOHMANN_JSON_SERIALIZE_ENUM( ThreadModel,
                                     {
                                         { ThreadModel::StdLib, "StdLib" },
                                         { ThreadModel::TBB, "TBB" },
                                     } )

//! \endcond

/**
   \class HistInfo
   Information about a given histogram: variable name(s), daughter indices
*/
class HistInfo {
  public:
    /**
       \enum Axis
       Enumeration of histogram axes
    */
    enum class Axis
    {
        X,
        Y
    };

    /// Constructor for 1D histogram
    HistInfo( const std::string& name, const int d1, const int d2 ) :
        m_nameX{ name }, m_d1X{ d1 }, m_d2X{ d2 }
    {
    }

    /// Constructor for 2D histogram
    HistInfo( const std::string& nameX, const int d1X, const int d2X,
              const std::string& nameY, const int d1Y, const int d2Y ) :
        m_nameX{ nameX }, m_nameY{ nameY }, m_d1X{ d1X }, m_d2X{ d2X }, m_d1Y{ d1Y }, m_d2Y{ d2Y }
    {
    }

    /// Retrieve variable name for given axis
    const std::string getName( const Axis var = Axis::X ) const
    {
        return var == Axis::Y ? m_nameY : m_nameX;
    }

    /// Retrieve daughter 1 index for given axis
    int getd1( const Axis var = Axis::X ) const
    {
        return var == Axis::Y ? m_d1Y : m_d1X;
    }

    /// Retrieve daughter 2 index for given axis
    int getd2( const Axis var = Axis::X ) const
    {
        return var == Axis::Y ? m_d2Y : m_d2X;
    }

  private:
    /// x-axis variable name
    std::string m_nameX;
    /// y-axis variable name
    std::string m_nameY;
    /// x-axis daughter 1 index
    int m_d1X{ 0 };
    /// x-axis daughter 2 index
    int m_d2X{ 0 };
    /// y-axis daughter 1 index
    int m_d1Y{ 0 };
    /// y-axis daughter 2 index
    int m_d2Y{ 0 };
};

/**
   \class TestHistos
   \brief The set of histograms to be populated for a given test

   Defined such that each thread can have its own set of histograms, which can
   then be combined by the main thread after synchronisation and saved to the
   output file.

   To work properly with tbb::parallel_reduce, a default-constructed TestHistos
   instance must act like a zero in a summation operation.
*/
class TestHistos {
  public:
    /// Default constructor
    TestHistos() = default;
    /// Constructor - creates the histograms specified in the JSON configuration
    TestHistos( const std::string& parentName, const nlohmann::json& config );
    /// Copy constructor
    TestHistos( const TestHistos& rhs );
    /// Move constructor
    TestHistos( TestHistos&& rhs ) noexcept;
    /// Copy asignment operator
    TestHistos& operator=( const TestHistos& rhs );
    /// Move asignment operator
    TestHistos& operator=( TestHistos&& rhs ) noexcept;
    /// Destructor
    ~TestHistos() = default;

    /// Add the rhs histograms to our own
    void add( const TestHistos& rhs );

    /// Normalise all histograms
    void normalise();

    /// Write all histograms to the supplied ROOT file
    void save( TFile* outputFile );

    /// typedef for 1D histogram store
    using HistStore1D = std::vector<std::pair<HistInfo, std::unique_ptr<TH1>>>;
    /// typedef for 2D histogram store
    using HistStore2D = std::vector<std::pair<HistInfo, std::unique_ptr<TH2>>>;

    /// Retrieve a const reference to the 1D histogram store
    const HistStore1D& get1DHistos() const { return m_1DhistVect; }
    /// Retrieve a reference to the 1D histogram store
    HistStore1D& get1DHistos() { return m_1DhistVect; }
    /// Retrieve a const reference to the 2D histogram store
    const HistStore2D& get2DHistos() const { return m_2DhistVect; }
    /// Retrieve a reference to the 2D histogram store
    HistStore2D& get2DHistos() { return m_2DhistVect; }
    /// Retrieve a const pointer to the mixed histogram
    const TH1* getMixedHist() const { return m_mixedHist.get(); }
    /// Retrieve a pointer to the mixed histogram
    TH1* getMixedHist() { return m_mixedHist.get(); }

  private:
    /// Swap the histogram stores between two instances
    void swap( TestHistos& rhs ) noexcept;

    /// The store of 1D histograms
    HistStore1D m_1DhistVect;
    /// The store of 2D histograms
    HistStore2D m_2DhistVect;
    /// The mixed histogram
    std::unique_ptr<TH1> m_mixedHist;
};

/**
  \struct TestConfig
  The settings for the test, as read from the JSON file
*/
struct TestConfig {
    /// Name of the parent particle
    std::string parentName;
    /// Names of all first-generation decay products
    std::vector<std::string> daughterNames;
    /// Names of all second-generation decay products
    std::vector<std::vector<std::string>> grandDaughterNames;

    /// Names of the decay models
    std::vector<std::string> modelNames;
    /// Parameters of the decay models
    std::vector<std::vector<std::string>> modelParameters;

    /// Flag to control whether the charge-conjugate decay shoudl be generated
    std::vector<bool> doConjDecay;

    /// Extra global commands, e.g. noFSR
    std::vector<std::string> extraCommands;

    /// The FSR generator to use
    FSRGenerator fsrGenerator;

    /// The threading model to use
    ThreadModel threadModel;

    /// The name of the output ROOT file
    std::string outFileName;
    /// The name of the reference ROOT file
    std::string refFileName;
    /// The name of the decay file
    std::string decFileName;

    /// The number of events to generate
    std::size_t nEvents;
    /// The number of threads to use
    std::size_t nThreads;
    /// The base of the random number seed
    std::size_t rngSeed;

    /// Flag to control printing of debug output
    bool debugFlag;

    /// The set of histograms to be filled
    TestHistos testHistograms;
};

/**
   \class TestDecayModel
   \brief Setup and run a test of a particular decay chain

   The test configuration is loaded from a JSON object.
   This contains the configuration of the model to test and the histograms to
   fill, plus the configuration of how to run the test (FSR settings, thread
   model, number of threads).
*/
class TestDecayModel {
  public:
    /// Constructor
    explicit TestDecayModel( const nlohmann::json& config );

    /// Run the test
    void run();

  private:
    /// Check that all mandatory fields are present in the JSON record
    static bool checkMandatoryFields( const nlohmann::json& config );

    /// Populate the test configuration from the JSON record
    static TestConfig readConfig( const nlohmann::json& config );

    /// Run multi-threaded with standard library threading model
    TestHistos runStdThreads() const;

    /// Run multi-threaded with Intel Thread Building Blocks threading model
    TestHistos runTBBThreads() const;

    /**
       For the specified range of events, run the main body: initialise,
       create and read the decay file, generate events
    */
    TestHistos runDecayBody( const std::size_t firstEvent,
                             const std::size_t nEvents ) const;

    /// Create the decay file for the decay chain under test
    void createDecFile() const;

    /// Generate the specified range of events and fill the histograms
    void generateEvents( EvtGen& theGen, TestHistos& theHistos,
                         const std::size_t firstEvent,
                         const std::size_t nEvents ) const;

    /// Calculate the value of the given variable for the generated decay
    double getValue( const EvtParticle* rootPart, const std::string& varName,
                     const int d1, const int d2 ) const;

    /// Compare generated histograms with those from a reference file
    void compareHistos( const TestHistos& theHistos,
                        const std::string& refFileName ) const;

    /**
       \brief Find the index of the daughter with the largest energy

       Sign convention: take negative daughter if mother is neutral,
       otherwise the one with the same sign as the mother.
    */
    int findChargedDaughterWithMaxE( const EvtParticle* parent ) const;

    /// Find the particle group for the given particle ID
    int getPartGroup( const int PDGId ) const;

    /// Check whether the given particle group is of the requested type
    bool isPartType( const int group, const std::string& particleType ) const;

    /// Calculate the cosine of the acoplanarity angle
    double getCosAcoplanarityAngle( const EvtParticle* selectedParent,
                                    const int sel_NDaugMax, const int d1,
                                    const int d2 ) const;

    /// The test configuration
    const TestConfig m_config;

    /// Thread safety mutex
    static std::once_flag m_createDecFile_threadlock;
};

#endif
