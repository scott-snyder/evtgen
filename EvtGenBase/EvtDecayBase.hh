
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

#ifndef EVTDECAYBASE_HH
#define EVTDECAYBASE_HH

#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtSpinType.hh"

#include <stdlib.h>
#include <string>
#include <vector>
class EvtParticle;
class EvtSpinType;

class EvtDecayBase {
  public:
    //These pure virtual methods has to be implemented
    //by any derived class
    virtual std::string getName() = 0;
    virtual void decay( EvtParticle* p ) = 0;
    virtual void makeDecay( EvtParticle* p, bool recursive = true ) = 0;
    virtual EvtDecayBase* clone() = 0;

    //These virtual methods can be implemented by the
    //derived class to implement nontrivial functionality.
    virtual void init();
    virtual void initProbMax();
    virtual std::string commandName();
    virtual void command( std::string cmd );

    virtual std::string getParamName( int i );
    virtual std::string getParamDefault( int i );

    double getProbMax( double prob );
    double resetProbMax( double prob );

    EvtDecayBase() = default;
    virtual ~EvtDecayBase() = default;

    virtual bool matchingDecay( const EvtDecayBase& other ) const;

    EvtId getParentId() const { return m_parent; }
    double getBranchingFraction() const { return m_brfr; }
    void disableCheckQ() { m_chkCharge = false; };
    void checkQ();
    int getNDaug() const { return m_ndaug; }
    EvtId* getDaugs() { return m_daug.data(); }
    EvtId getDaug( int i ) const { return m_daug[i]; }
    int getNArg() const { return m_narg; }
    bool getFSR() const { return m_fsr; }
    void setFSR() { m_fsr = true; }
    void setVerbose() { m_verbose = true; }
    void setSummary() { m_summary = true; }
    double* getArgs();
    std::string* getArgsStr() { return m_args.data(); }
    double getArg( unsigned int j );
    double getStoredArg( int j ) const { return m_storedArgs.at( j ); }
    double getNStoredArg() const { return m_storedArgs.size(); }
    std::string getArgStr( int j ) const { return m_args[j]; }
    std::string getModelName() const { return m_modelname; }
    int getDSum() const { return m_dsum; }
    bool summary() const { return m_summary; }
    bool verbose() const { return m_verbose; }

    void saveDecayInfo( EvtId ipar, int ndaug, EvtId* daug, int narg,
                        std::vector<std::string>& args, std::string name,
                        double brfr );
    void printSummary() const;
    void printInfo() const;

    //Does not really belong here but I don't have a better place.
    static void findMasses( EvtParticle* p, int ndaugs, EvtId daugs[10],
                            double masses[10] );
    static void findMass( EvtParticle* p );
    static double findMaxMass( EvtParticle* p );

    //Methods to set the maximum probability.
    void setProbMax( double prbmx );
    void noProbMax();

    void checkNArg( int a1, int a2 = -1, int a3 = -1, int a4 = -1 );
    void checkNDaug( int d1, int d2 = -1 );

    void checkSpinParent( EvtSpinType::spintype sp );
    void checkSpinDaughter( int d1, EvtSpinType::spintype sp );

    // lange - some models can take more daughters
    // than they really have to fool aliases (VSSBMIX for example)
    virtual int nRealDaughters() { return m_ndaug; }

  protected:
    bool m_daugsDecayedByParentModel;
    bool daugsDecayedByParentModel() { return m_daugsDecayedByParentModel; }

  private:
    std::vector<double> m_storedArgs;
    std::vector<EvtId> m_daug;
    std::vector<double> m_argsD;
    std::vector<std::string> m_args;

    std::string m_modelname = "**********";

    EvtId m_parent = EvtId( -1, -1 );
    int m_ndaug = 0;
    int m_narg = 0;
    double m_brfr = 0;
    int m_dsum = 0;

    bool m_fsr = false;
    bool m_summary = false;
    bool m_verbose = false;

    // The default is that the user module does _not_ set any probmax.
    bool m_defaultprobmax = true;
    int m_ntimes_prob = 0;
    double m_probmax = 0.0;

    //Should charge conservation be checked when model is created?
    //Default is to check that charge is conserved
    bool m_chkCharge = true;

    //These are used for gathering statistics.
    double m_sum_prob = 0.0;
    double m_max_prob = 0.0;
};

#endif
