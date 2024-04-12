
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

    EvtId getParentId() const { return _parent; }
    double getBranchingFraction() const { return _brfr; }
    void disableCheckQ() { _chkCharge = false; };
    void checkQ();
    int getNDaug() const { return _ndaug; }
    EvtId* getDaugs() { return _daug.data(); }
    EvtId getDaug( int i ) const { return _daug[i]; }
    int getNArg() const { return _narg; }
    bool getFSR() const { return _fsr; }
    void setFSR() { _fsr = true; }
    void setVerbose() { _verbose = true; }
    void setSummary() { _summary = true; }
    double* getArgs();
    std::string* getArgsStr() { return _args.data(); }
    double getArg( unsigned int j );
    double getStoredArg( int j ) const { return _storedArgs.at( j ); }
    double getNStoredArg() const { return _storedArgs.size(); }
    std::string getArgStr( int j ) const { return _args[j]; }
    std::string getModelName() const { return _modelname; }
    int getDSum() const { return _dsum; }
    bool summary() const { return _summary; }
    bool verbose() const { return _verbose; }

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
    virtual int nRealDaughters() { return _ndaug; }

  protected:
    bool _daugsDecayedByParentModel;
    bool daugsDecayedByParentModel() { return _daugsDecayedByParentModel; }

  private:
    std::vector<double> _storedArgs;
    std::vector<EvtId> _daug;
    std::vector<double> _argsD;
    std::vector<std::string> _args;

    std::string _modelname = "**********";

    EvtId _parent = EvtId( -1, -1 );
    int _ndaug = 0;
    int _narg = 0;
    double _brfr = 0;
    int _dsum = 0;

    bool _fsr = false;
    bool _summary = false;
    bool _verbose = false;

    // The default is that the user module does _not_ set any probmax.
    bool defaultprobmax = true;
    int ntimes_prob = 0;
    double probmax = 0.0;

    //Default is to check that charge is conserved
    bool _chkCharge = true;

    //These are used for gathering statistics.
    double sum_prob = 0.0;
    double max_prob = 0.0;
};

#endif
