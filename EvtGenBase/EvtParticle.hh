
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

#ifndef EVTPARTICLE_HH
#define EVTPARTICLE_HH

#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtSpinDensity.hh"
#include "EvtGenBase/EvtSpinType.hh"
#include "EvtGenBase/EvtVector4R.hh"

#include <assert.h>
#include <map>
#include <string>
#include <vector>

class EvtDiracSpinor;
class EvtVector4C;
class EvtTensor4C;
class EvtStdHep;
class EvtSecondary;
class EvtRaritaSchwinger;

const int MAX_DAUG = 100;
const int MAX_LEVEL = 10;
const int MAX_TRIES = 10000;

class EvtParticle {
  public:
    /**
  * Default constructor.
  */
    EvtParticle();

    /**
  * Destructor.
  */
    virtual ~EvtParticle();

    /**
  * Returns polarization vector in the parents restframe.
  */
    virtual EvtVector4C epsParent( int i ) const;

    /**
  * Returns polarization vector in the particles own restframe.
  */
    virtual EvtVector4C eps( int i ) const;

    /**
  * Returns polarization vector in the parents restframe for a photon.
  */
    virtual EvtVector4C epsParentPhoton( int i ) const;

    /**
  * Returns polarization vector in the particles own restframe for a photon.
  */
    virtual EvtVector4C epsPhoton( int i ) const;

    /**
  * Returns Dirac spinor in the parents restframe for a Dirac particle.
  */
    virtual EvtDiracSpinor spParent( int ) const;

    /**
  * Returns Dirac spinor in the particles own restframe for a Dirac particle.
  */
    virtual EvtDiracSpinor sp( int ) const;

    /**
  * Returns Dirac spinor in the parents restframe for a Neutrino particle.
  */
    virtual EvtDiracSpinor spParentNeutrino() const;

    /**
  * Returns Dirac spinor in the particles own restframe for a
  * Neutrino particle.
  */
    virtual EvtDiracSpinor spNeutrino() const;

    /**
  * Returns tensor in the parents restframe for a spin 2 particle.
  */
    virtual EvtTensor4C epsTensorParent( int i ) const;

    /**
  * Returns tensor in the particles own restframe for a spin 2 particle.
  */
    virtual EvtTensor4C epsTensor( int i ) const;

    /**
   * Returns Rarita-Schwinger spinor in the parents restframe for a 
   * Rarita-Schwinger particle.
   */
    virtual EvtRaritaSchwinger spRSParent( int ) const;

    /**
   * Returns Rarita-Schwinger spinor in the particles own restframe for a 
   * Rarita-Schwinger particle.
   */
    virtual EvtRaritaSchwinger spRS( int ) const;

    /**
  * Initialiaze particle with id and 4momentum.
  */
    virtual void init( EvtId part_n, const EvtVector4R& p4 ) = 0;

    /**
  * Add another daughter to the particle
  */
    void addDaug( EvtParticle* node );

    /**
  * Decay particle
  */
    void decay();

    /** 
  * Delete a decay chain
  */
    void deleteTree();
    void deleteDaughters( bool keepChannel = false );

    /**
  * Should only be used internally.
  */
    void setChannel( int i );

    /**
  * Creates the daughters in the list of ids and 
  * adds them to the parent. Note that momentum
  * is left uninitialized, this is _only_ creation.
  */
    void makeDaughters( size_t ndaug, EvtId* id );

    /**
  * Creates the daughters in the list of ids and 
  * adds them to the parent. Note that momentum
  * is left uninitialized, this is _only_ creation.
  */
    void makeDaughters( size_t ndaug, std::vector<EvtId> idVector );

    /**
  * Similar to the routine above except that here 
  * momentum is generated according to phase space 
  * daughters are filled with this momentum.
  */
    double initializePhaseSpace( size_t numdaughter, EvtId* daughters,
                                 bool forceResetMasses = false,
                                 double poleSize = -1., int whichTwo1 = 0,
                                 int whichTwo2 = 1 );

    /**
  * Get pointer the the i:th daugther.
  */
    EvtParticle* getDaug( const int i ) { return m_daug[i]; }

    /**
  * Get const pointer the the i:th daugther.
  */
    const EvtParticle* getDaug( const int i ) const { return m_daug[i]; }

    /**
  * Iterates over the particles in a decay chain.
  */
    EvtParticle* nextIter( EvtParticle* rootOfTree = nullptr );

    /**
  * Makes stdhep list
  */
    void makeStdHep( EvtStdHep& stdhep, EvtSecondary& secondary,
                     EvtId* stable_parent_ihep );
    void makeStdHep( EvtStdHep& stdhep );

    /**
  * Gets 4vector in the labframe, i.e., the frame in which the root
  * particles momentum is measured.
  */
    EvtVector4R getP4Lab() const;

    /**
  * Gets 4vector in the labframe for the 4-momentum before FSR was 
  * generated in the parents decay. The lab frame is where the root
  * particles momentum is measured.
  */
    EvtVector4R getP4LabBeforeFSR() const;

    /**
  * Gets 4vector in the particles restframe, i.e. this functiont will
  * return (m,0,0,0)
  */
    EvtVector4R getP4Restframe() const;

    /**
  * Returns the 4position of the particle in the lab frame.
  */
    EvtVector4R get4Pos() const;

    /**
  * Returns pointer to parent particle.
  */
    EvtParticle* getParent() const;

    /**
  * Makes partptr the idaug:th daugther.
  */
    void insertDaugPtr( int idaug, EvtParticle* partptr )
    {
        m_daug[idaug] = partptr;
        partptr->m_parent = this;
    }
    /**
  * Returns mass of particle.
  */
    double mass() const;

    /**
  * Used internally to decide if first time particle is decayed.
  */
    int firstornot() const;
    void setFirstOrNot();
    void resetFirstOrNot();

    /**
  * Returns Id of particle.
  */
    EvtId getId() const;

    /**
  * Returns the PDG id of the particle
  */

    int getPDGId() const;

    /** 
  * Returns particle type.
  */

    EvtSpinType::spintype getSpinType() const;

    /**
  * Returns number of spin states of the particle.
  */
    int getSpinStates() const;

    /**
  * Returns 4momentum in parents restframe.
  */
    const EvtVector4R& getP4() const;

    /**
  * Sets the 4momentum in the parents restframe.
  */
    void setP4( const EvtVector4R& p4 )
    {
        m_p = p4;
        m_pBeforeFSR = p4;
    }

    void setP4WithFSR( const EvtVector4R& p4 ) { m_p = p4; }

    void setFSRP4toZero() { m_pBeforeFSR.set( 0.0, 0.0, 0.0, 0.0 ); }

    /**
  * Retunrs the decay channel.
  */
    int getChannel() const;

    /**
  * Returns number of daugthers.
  */
    size_t getNDaug() const;
    void resetNDaug()
    {
        m_ndaug = 0;
        return;
    }

    /**
  * Prints out the particle "tree" of a given particle.  The
  * tree consists of all daughters (and their daughters, etc)
  * and their properties.
  */
    void printTree() const;

    void printTreeRec( unsigned int level ) const;

    std::string treeStr() const;
    std::string treeStrRec( unsigned int level ) const;

    /**
  * Prints information for the particle.
  */
    void printParticle() const;

    /**
  * Set lifetime of the particle in parents restframe.
  */
    void setLifetime( double tau );

    /**
  * Generate lifetime according to pure exponential.
  */
    void setLifetime();

    /**
  * Returns the lifetime.
  */
    double getLifetime() const;

    /** 
  * Set diagonal spindensity matrix.
  */
    void setDiagonalSpinDensity();

    /** 
  * Set spindensity matrix for e+e- -> V
  */
    void setVectorSpinDensity();

    /**
  * Set forward spin density matrix.
  */
    void setSpinDensityForward( const EvtSpinDensity& rho )
    {
        m_rhoForward = rho;
    }

    /**
  * Set forward spin density matrix according to the density matrix
  * rho in the helicity amplitude basis.
  */
    void setSpinDensityForwardHelicityBasis( const EvtSpinDensity& rho );
    void setSpinDensityForwardHelicityBasis( const EvtSpinDensity& rho,
                                             double alpha, double beta,
                                             double gamma );

    /**
  * Returns a rotation matrix need to rotate the basis state
  * to the helicity basis. The EvtSpinDensity matrix is just use
  * as a matrix here. This function is to be implemented in each
  * derived class.
  */
    virtual EvtSpinDensity rotateToHelicityBasis() const = 0;
    virtual EvtSpinDensity rotateToHelicityBasis( double alpha, double beta,
                                                  double gamma ) const = 0;

    /**
  * Get forward spin density matrix.
  */
    EvtSpinDensity getSpinDensityForward() { return m_rhoForward; }

    /**
  * Set backward spin density matrix.
  */
    void setSpinDensityBackward( const EvtSpinDensity& rho )
    {
        m_rhoBackward = rho;
    }

    /**
  * Get backward spin density matrix.
  */
    EvtSpinDensity getSpinDensityBackward() { return m_rhoBackward; }

    //Hacks will be removed when better solutions are thought of!
    //This is used to suppress use of random numbers when doing initialization
    //of some models.
    void noLifeTime() { m_genlifetime = 0; }

    //lange - April 29, 2002
    void setId( EvtId id ) { m_id = id; }
    void initDecay( bool useMinMass = false );
    bool generateMassTree();

    double compMassProb();

    //setMass will blow away any existing 4vector
    void setMass( double m ) { m_p = EvtVector4R( m, 0.0, 0.0, 0.0 ); }

    //void setMixed() {m_mix=true;}
    //void setUnMixed() {m_mix=false;}
    //bool getMixed() {return m_mix;}

    //void takeCConj() {EvtGenReport(EVTGEN_INFO,"EvtGen") << "should take conj\n";}

    //this means that the particle has gone through initDecay
    // and thus has a mass
    bool isInitialized() { return m_isInit; }
    bool hasValidP4() { return m_validP4; }
    bool isDecayed() { return m_isDecayed; }

    // decay prob - only relevent if already decayed
    // and is a scalar particle
    // returned is a double* that should be prob/probMax
    // FIXME - this should probably be changed to std::optional
    const double* decayProb() const { return m_decayProb; }
    void setDecayProb( double p );

    // Return the name of the particle (from the EvtId number)
    std::string getName() const;

    // Specify whether the particle has a special named attribute with
    // a set value. By default, nothing is set, but derived classes
    // can set this to mean something specific, e.g. if a photon is FSR
    void setAttribute( std::string attName, int attValue )
    {
        m_intAttributes[attName] = attValue;
    }

    // Retrieve the integer value for the given attribute name
    int getAttribute( std::string attName ) const;

    // Specify if the particle has a double attribute value, e.g. amplitude weight.
    // By default, nothing is set, but derived classes can set this to mean something specific
    void setAttributeDouble( std::string attName, double attValue )
    {
        m_dblAttributes[attName] = attValue;
    }

    // Retrieve the double value for the given attribute name
    double getAttributeDouble( std::string attName ) const;

  protected:
    void setp( double e, double px, double py, double pz )
    {
        m_p.set( e, px, py, pz );
        m_pBeforeFSR = m_p;
    }

    void setp( const EvtVector4R& p4 )
    {
        m_p = p4;
        m_pBeforeFSR = m_p;
    }

    void setpart_num( EvtId particle_number )
    {
        assert( m_channel == -10 || m_id.getId() == particle_number.getId() ||
                m_id.getId() == -1 );
        m_id = particle_number;
    }
    bool m_validP4;

    // A typedef to define the attribute (name, integer) map
    typedef std::map<std::string, int> EvtAttIntMap;
    EvtAttIntMap m_intAttributes;

    // A typedef to define the attribute (name, double) map
    typedef std::map<std::string, double> EvtAttDblMap;
    EvtAttDblMap m_dblAttributes;

  private:
    EvtParticle* m_daug[MAX_DAUG];
    size_t m_ndaug;
    EvtParticle* m_parent;
    int m_channel;
    int m_first;
    EvtId m_id;
    EvtVector4R m_p;
    EvtVector4R m_pBeforeFSR;
    double m_t;
    bool m_isInit;
    bool m_isDecayed;

    //bool m_mix;

    EvtSpinDensity m_rhoForward;
    EvtSpinDensity m_rhoBackward;

    void makeStdHepRec( int firstparent, int lastparent, EvtStdHep& stdhep,
                        EvtSecondary& secondary, EvtId* stable_parent_ihep );
    void makeStdHepRec( int firstparent, int lastparent, EvtStdHep& stdhep );

    //This is a hack until things gets straightened out. (Ryd)
    int m_genlifetime;

    //should never be used, therefor is private.
    //these does _not_ have an implementation
    EvtParticle& operator=( const EvtParticle& p );
    EvtParticle( const EvtParticle& p );

    double* m_decayProb;
};

#endif
