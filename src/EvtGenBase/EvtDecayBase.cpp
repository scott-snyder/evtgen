
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

#include "EvtGenBase/EvtDecayBase.hh"

#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtSpinType.hh"
#include "EvtGenBase/EvtStatus.hh"

#include <ctype.h>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <vector>
using std::endl;
using std::fstream;
void EvtDecayBase::checkQ()
{
    int i;
    int q = 0;
    int qpar;

    //If there are no daughters (jetset etc) then we do not
    //want to do this test.  Why?  Because sometimes the parent
    //will have a nonzero charge.

    if ( m_ndaug != 0 ) {
        for ( i = 0; i < m_ndaug; i++ ) {
            q += EvtPDL::chg3( m_daug[i] );
        }
        qpar = EvtPDL::chg3( m_parent );

        if ( q != qpar ) {
            EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                << m_modelname.c_str() << " generator expected "
                << " charge to be conserved, found:" << endl;
            EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                << "Parent charge of " << ( qpar / 3 ) << endl;
            EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                << "Sum of daughter charge of " << ( q / 3 ) << endl;
            EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                << "The parent is " << EvtPDL::name( m_parent ).c_str() << endl;
            for ( i = 0; i < m_ndaug; i++ ) {
                EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                    << "Daughter " << EvtPDL::name( m_daug[i] ).c_str() << endl;
            }
            EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                << "Will terminate execution!" << endl;

            ::abort();
        }
    }
}

double EvtDecayBase::getProbMax( double prob )
{
    int i;

    //diagnostics
    m_sum_prob += prob;
    if ( prob > m_max_prob )
        m_max_prob = prob;

    if ( m_defaultprobmax && m_ntimes_prob <= 500 ) {
        //We are building up m_probmax with this iteration
        m_ntimes_prob += 1;
        if ( prob > m_probmax ) {
            m_probmax = prob;
        }
        if ( m_ntimes_prob == 500 ) {
            m_probmax *= 1.2;
        }
        return 1000000.0 * prob;
    }

    if ( prob > m_probmax * 1.0001 ) {
        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << "prob > probmax:(" << prob << ">" << m_probmax << ")";
        EvtGenReport( EVTGEN_INFO, "" ) << "(" << m_modelname.c_str() << ") ";
        EvtGenReport( EVTGEN_INFO, "" )
            << EvtPDL::name( m_parent ).c_str() << " -> ";
        for ( i = 0; i < m_ndaug; i++ ) {
            EvtGenReport( EVTGEN_INFO, "" )
                << EvtPDL::name( m_daug[i] ).c_str() << " ";
        }
        EvtGenReport( EVTGEN_INFO, "" ) << endl;

        if ( m_defaultprobmax )
            m_probmax = prob;
    }

    m_ntimes_prob += 1;

    return m_probmax;

}    //getProbMax

double EvtDecayBase::resetProbMax( double prob )
{
    EvtGenReport( EVTGEN_INFO, "EvtGen" ) << "Reseting prob max\n";
    EvtGenReport( EVTGEN_INFO, "EvtGen" )
        << "prob > probmax:(" << prob << ">" << m_probmax << ")";
    EvtGenReport( EVTGEN_INFO, "" ) << "(" << m_modelname.c_str() << ")";
    EvtGenReport( EVTGEN_INFO, "" ) << EvtPDL::getStdHep( m_parent ) << "->";

    for ( int i = 0; i < m_ndaug; i++ ) {
        EvtGenReport( EVTGEN_INFO, "" ) << EvtPDL::getStdHep( m_daug[i] ) << " ";
    }
    EvtGenReport( EVTGEN_INFO, "" ) << endl;

    m_probmax = 0.0;
    m_defaultprobmax = false;
    m_ntimes_prob = 0;

    return prob;
}

std::string EvtDecayBase::commandName()
{
    return std::string( "" );
}

void EvtDecayBase::command( std::string )
{
    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "Should never call EvtDecayBase::command" << endl;
    ::abort();
}

std::string EvtDecayBase::getParamName( int i )
{
    switch ( i ) {
        case 0:
            return "param00";
        case 1:
            return "param01";
        case 2:
            return "param02";
        case 3:
            return "param03";
        case 4:
            return "param04";
        case 5:
            return "param05";
        case 6:
            return "param06";
        case 7:
            return "param07";
        case 8:
            return "param08";
        case 9:
            return "param09";
        default:
            return "";
    }
}

std::string EvtDecayBase::getParamDefault( int /*i*/ )
{
    return "";
}

void EvtDecayBase::init()
{
    //This default version of init does nothing;
    //A specialized version of this function can be
    //supplied for each decay model to do initialization.

    return;
}

void EvtDecayBase::initProbMax()
{
    //This function is called if the decay does not have a
    //specialized initialization.
    //The default is to set the maximum
    //probability to 0 and the number of times called to 0
    //and m_defaultprobmax to 'true' such that the decay will be
    //generated many many times
    //in order to generate a reasonable maximum probability
    //for the decay.

    m_defaultprobmax = true;
    m_ntimes_prob = 0;
    m_probmax = 0.0;

}    //initProbMax

void EvtDecayBase::saveDecayInfo( EvtId ipar, int ndaug, EvtId* daug, int narg,
                                  std::vector<std::string>& args,
                                  std::string name, double brfr )
{
    int i;

    m_brfr = brfr;
    m_ndaug = ndaug;
    m_narg = narg;
    m_parent = ipar;

    m_dsum = 0;

    if ( m_ndaug > 0 ) {
        m_daug.resize( m_ndaug );
        for ( i = 0; i < m_ndaug; i++ ) {
            m_daug[i] = daug[i];
            m_dsum += daug[i].getAlias();
        }
    } else {
        m_daug.clear();
    }

    if ( m_narg > 0 ) {
        m_args.resize( m_narg + 1 );
        for ( i = 0; i < m_narg; i++ ) {
            m_args[i] = args[i];
        }
    } else {
        m_args.clear();
    }

    m_modelname = name;

    this->init();
    this->initProbMax();

    if ( m_chkCharge ) {
        this->checkQ();
    }

    if ( m_defaultprobmax ) {
        EvtGenReport( EVTGEN_INFO, "EvtGen" ) << "No default probmax for ";
        EvtGenReport( EVTGEN_INFO, "" ) << "(" << m_modelname.c_str() << ") ";
        EvtGenReport( EVTGEN_INFO, "" )
            << EvtPDL::name( m_parent ).c_str() << " -> ";
        for ( i = 0; i < m_ndaug; i++ ) {
            EvtGenReport( EVTGEN_INFO, "" )
                << EvtPDL::name( m_daug[i] ).c_str() << " ";
        }
        EvtGenReport( EVTGEN_INFO, "" ) << endl;
        EvtGenReport( EVTGEN_INFO, "" )
            << "This is fine for development, but must be provided for production."
            << endl;
        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << "Never fear though - the decay will use the \n";
        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << "500 iterations to build up a good probmax \n";
        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << "before accepting a decay. " << endl;
    }
}

void EvtDecayBase::printSummary() const
{
    if ( m_ntimes_prob > 0 ) {
        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << "Calls = " << m_ntimes_prob
            << " eff: " << m_sum_prob / ( m_probmax * m_ntimes_prob )
            << " frac. max:" << m_max_prob / m_probmax;
        EvtGenReport( EVTGEN_INFO, "" )
            << " probmax:" << m_probmax << " max:" << m_max_prob << " : ";
    }

    printInfo();
}

void EvtDecayBase::printInfo() const
{
    EvtGenReport( EVTGEN_INFO, "" ) << EvtPDL::name( m_parent ).c_str() << " -> ";
    for ( int i = 0; i < m_ndaug; i++ ) {
        EvtGenReport( EVTGEN_INFO, "" )
            << EvtPDL::name( m_daug[i] ).c_str() << " ";
    }
    EvtGenReport( EVTGEN_INFO, "" ) << " (" << m_modelname.c_str() << ")" << endl;
}

void EvtDecayBase::setProbMax( double prbmx )
{
    m_defaultprobmax = false;
    m_probmax = prbmx;
}

void EvtDecayBase::noProbMax()
{
    m_defaultprobmax = false;
}

double EvtDecayBase::findMaxMass( EvtParticle* p )
{
    double maxOkMass = EvtPDL::getMaxMass( p->getId() );

    //protect against vphotons
    if ( maxOkMass < 0.0000000001 )
        return 10000000.;
    //and against already determined masses
    if ( p->hasValidP4() )
        maxOkMass = p->mass();

    EvtParticle* par = p->getParent();
    if ( par ) {
        double maxParMass = findMaxMass( par );
        size_t i;
        double minDaugMass = 0.;
        for ( i = 0; i < par->getNDaug(); i++ ) {
            EvtParticle* dau = par->getDaug( i );
            if ( dau != p ) {
                // it might already have a mass
                if ( dau->isInitialized() || dau->hasValidP4() )
                    minDaugMass += dau->mass();
                else
                    //give it a bit of phase space
                    minDaugMass += 1.000001 * EvtPDL::getMinMass( dau->getId() );
            }
        }
        if ( maxOkMass > ( maxParMass - minDaugMass ) )
            maxOkMass = maxParMass - minDaugMass;
    }
    return maxOkMass;
}

// given list of daughters ( by number ) returns a
// list of viable masses.

void EvtDecayBase::findMass( EvtParticle* p )
{
    //Need to also check that this mass does not screw
    //up the parent
    //This code assumes that for the ith daughter, 0..i-1
    //already have a mass
    double maxOkMass = findMaxMass( p );

    int count = 0;
    double mass;
    bool massOk = false;
    size_t i;
    while ( !massOk ) {
        count++;
        if ( count > 10000 ) {
            EvtGenReport( EVTGEN_INFO, "EvtGen" )
                << "Can not find a valid mass for: "
                << EvtPDL::name( p->getId() ).c_str() << endl;
            EvtGenReport( EVTGEN_INFO, "EvtGen" )
                << "Now printing parent and/or grandparent tree\n";
            if ( p->getParent() ) {
                if ( p->getParent()->getParent() ) {
                    p->getParent()->getParent()->printTree();
                    EvtGenReport( EVTGEN_INFO, "EvtGen" )
                        << p->getParent()->getParent()->mass() << endl;
                    EvtGenReport( EVTGEN_INFO, "EvtGen" )
                        << p->getParent()->mass() << endl;
                } else {
                    p->getParent()->printTree();
                    EvtGenReport( EVTGEN_INFO, "EvtGen" )
                        << p->getParent()->mass() << endl;
                }
            } else
                p->printTree();
            EvtGenReport( EVTGEN_INFO, "EvtGen" )
                << "maxokmass=" << maxOkMass << " "
                << EvtPDL::getMinMass( p->getId() ) << " "
                << EvtPDL::getMaxMass( p->getId() ) << endl;
            if ( p->getNDaug() ) {
                for ( i = 0; i < p->getNDaug(); i++ ) {
                    EvtGenReport( EVTGEN_INFO, "EvtGen" )
                        << p->getDaug( i )->mass() << " ";
                }
                EvtGenReport( EVTGEN_INFO, "EvtGen" ) << endl;
            }
            if ( maxOkMass >= EvtPDL::getMinMass( p->getId() ) ) {
                EvtGenReport( EVTGEN_INFO, "EvtGen" )
                    << "taking a default value\n";
                p->setMass( maxOkMass );
                return;
            }
            assert( 0 );
        }
        mass = EvtPDL::getMass( p->getId() );
        //Just need to check that this mass is > than
        //the mass of all daughters
        double massSum = 0.;
        if ( p->getNDaug() ) {
            for ( i = 0; i < p->getNDaug(); i++ ) {
                massSum += p->getDaug( i )->mass();
            }
        }
        //some special cases are handled with 0 (stable) or 1 (k0->ks/kl) daughters
        if ( p->getNDaug() < 2 )
            massOk = true;
        if ( p->getParent() ) {
            if ( p->getParent()->getNDaug() == 1 )
                massOk = true;
        }
        if ( !massOk ) {
            if ( massSum < mass )
                massOk = true;
            if ( mass > maxOkMass )
                massOk = false;
        }
    }

    p->setMass( mass );
}

void EvtDecayBase::findMasses( EvtParticle* p, int ndaugs, EvtId daugs[10],
                               double masses[10] )
{
    int i;
    double mass_sum;

    int count = 0;

    if ( !( p->firstornot() ) ) {
        for ( i = 0; i < ndaugs; i++ ) {
            masses[i] = p->getDaug( i )->mass();
        }    //for
    }        //if
    else {
        p->setFirstOrNot();
        // if only one daughter do it

        if ( ndaugs == 1 ) {
            masses[0] = p->mass();
            return;
        }

        //until we get a combo whose masses are less than the parent mass.
        do {
            mass_sum = 0.0;

            for ( i = 0; i < ndaugs; i++ ) {
                masses[i] = EvtPDL::getMass( daugs[i] );
                mass_sum = mass_sum + masses[i];
            }

            count++;

            if ( count == 10000 ) {
                EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                    << "Decaying particle:" << EvtPDL::name( p->getId() ).c_str()
                    << " (m=" << p->mass() << ")" << endl;
                EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                    << "To the following daugthers" << endl;
                for ( i = 0; i < ndaugs; i++ ) {
                    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                        << EvtPDL::name( daugs[i] ).c_str() << endl;
                }
                EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                    << "Has been rejected " << count
                    << " times, will now take minimal masses "
                    << " of daugthers" << endl;

                mass_sum = 0.;
                for ( i = 0; i < ndaugs; i++ ) {
                    masses[i] = EvtPDL::getMinMass( daugs[i] );
                    mass_sum = mass_sum + masses[i];
                }
                if ( mass_sum > p->mass() ) {
                    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                        << "Parent mass=" << p->mass()
                        << "to light for daugthers." << endl
                        << "Will throw the event away." << endl;
                    //dont terminate - start over on the event.
                    EvtStatus::setRejectFlag();
                    mass_sum = 0.;
                    //	  ::abort();
                }
            }
        } while ( mass_sum > p->mass() );
    }    //else

    return;
}

void EvtDecayBase::checkNArg( int a1, int a2, int a3, int a4 )
{
    if ( m_narg != a1 && m_narg != a2 && m_narg != a3 && m_narg != a4 ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << m_modelname.c_str() << " generator expected " << endl;
        EvtGenReport( EVTGEN_ERROR, "EvtGen" ) << a1 << endl;
        ;
        if ( a2 > -1 ) {
            EvtGenReport( EVTGEN_ERROR, "EvtGen" ) << " or " << a2 << endl;
        }
        if ( a3 > -1 ) {
            EvtGenReport( EVTGEN_ERROR, "EvtGen" ) << " or " << a3 << endl;
        }
        if ( a4 > -1 ) {
            EvtGenReport( EVTGEN_ERROR, "EvtGen" ) << " or " << a4 << endl;
        }
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << " arguments but found:" << m_narg << endl;
        printSummary();
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Will terminate execution!" << endl;
        ::abort();
    }
}
void EvtDecayBase::checkNDaug( int d1, int d2 )
{
    if ( m_ndaug != d1 && m_ndaug != d2 ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << m_modelname.c_str() << " generator expected ";
        EvtGenReport( EVTGEN_ERROR, "EvtGen" ) << d1;
        if ( d2 > -1 ) {
            EvtGenReport( EVTGEN_ERROR, "EvtGen" ) << " or " << d2;
        }
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << " daughters but found:" << m_ndaug << endl;
        printSummary();
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Will terminate execution!" << endl;
        ::abort();
    }
}

void EvtDecayBase::checkSpinParent( EvtSpinType::spintype sp )
{
    EvtSpinType::spintype parenttype = EvtPDL::getSpinType( getParentId() );
    if ( parenttype != sp ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << m_modelname.c_str() << " did not get the correct parent spin\n";
        printSummary();
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Will terminate execution!" << endl;
        ::abort();
    }
}

void EvtDecayBase::checkSpinDaughter( int d1, EvtSpinType::spintype sp )
{
    EvtSpinType::spintype parenttype = EvtPDL::getSpinType( getDaug( d1 ) );
    if ( parenttype != sp ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << m_modelname.c_str()
            << " did not get the correct daughter spin d=" << d1 << endl;
        printSummary();
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Will terminate execution!" << endl;
        ::abort();
    }
}

double* EvtDecayBase::getArgs()
{
    if ( !m_argsD.empty() )
        return m_argsD.data();
    //The user has asked for a list of doubles - the arguments
    //better all be doubles...
    if ( m_narg == 0 )
        return m_argsD.data();

    m_argsD.resize( m_narg );
    for ( int i = 0; i < m_narg; i++ ) {
        char* tc;
        m_argsD[i] = strtod( m_args[i].c_str(), &tc );
    }
    return m_argsD.data();
}

double EvtDecayBase::getArg( unsigned int j )
{
    // Verify string
    const char* str = m_args[j].c_str();
    int i = 0;
    while ( str[i] != 0 ) {
        if ( isalpha( str[i] ) && str[i] != 'e' ) {
            EvtGenReport( EVTGEN_INFO, "EvtGen" )
                << "String " << str << " is not a number" << endl;
            assert( 0 );
        }
        i++;
    }

    char** tc = nullptr;
    double result = strtod( m_args[j].c_str(), tc );

    if ( m_storedArgs.size() < j + 1 ) {    // then store the argument's value
        m_storedArgs.push_back( result );
    }

    return result;
}

bool EvtDecayBase::matchingDecay( const EvtDecayBase& other ) const
{
    if ( m_ndaug != other.m_ndaug )
        return false;
    if ( m_parent != other.m_parent )
        return false;

    std::vector<int> useDs;
    for ( int i = 0; i < m_ndaug; i++ )
        useDs.push_back( 0 );

    for ( int i = 0; i < m_ndaug; i++ ) {
        bool foundIt = false;
        for ( int j = 0; j < m_ndaug; j++ ) {
            if ( useDs[j] == 1 )
                continue;
            if ( m_daug[i] == other.m_daug[j] &&
                 m_daug[i].getAlias() == other.m_daug[j].getAlias() ) {
                foundIt = true;
                useDs[j] = 1;
                break;
            }
        }
        if ( foundIt == false )
            return false;
    }
    for ( int i = 0; i < m_ndaug; i++ )
        if ( useDs[i] == 0 )
            return false;

    return true;
}
