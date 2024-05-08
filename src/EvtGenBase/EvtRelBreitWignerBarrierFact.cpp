
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

#include "EvtGenBase/EvtRelBreitWignerBarrierFact.hh"

#include "EvtGenBase/EvtAmpPdf.hh"
#include "EvtGenBase/EvtIntervalFlatPdf.hh"
#include "EvtGenBase/EvtMassAmp.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtPredGen.hh"
#include "EvtGenBase/EvtPropBreitWignerRel.hh"
#include "EvtGenBase/EvtSpinType.hh"
#include "EvtGenBase/EvtTwoBodyVertex.hh"

#include <algorithm>

EvtRelBreitWignerBarrierFact::EvtRelBreitWignerBarrierFact(
    double mass, double width, double maxRange, EvtSpinType::spintype sp ) :
    EvtAbsLineShape( mass, width, maxRange, sp )
{    // double mDaug1, double mDaug2, int l) {

    m_includeDecayFact = true;
    m_includeBirthFact = true;
    m_mass = mass;
    m_width = width;
    m_spin = sp;
    m_blattDecay = 3.0;
    m_blattBirth = 1.0;
    m_maxRange = maxRange;
    m_errorCond = false;

    double maxdelta = 15.0 * width;

    if ( maxRange > 0.00001 ) {
        m_massMax = mass + maxdelta;
        m_massMin = mass - maxRange;
    } else {
        m_massMax = mass + maxdelta;
        m_massMin = mass - 15.0 * width;
    }

    m_massMax = mass + maxdelta;
    if ( m_massMin < 0. ) {
        if ( m_width > 0.0001 ) {
            m_massMin = 0.00011;
        } else {
            m_massMin = 0.;
        }
    }
}

EvtRelBreitWignerBarrierFact::EvtRelBreitWignerBarrierFact(
    const EvtRelBreitWignerBarrierFact& x ) :
    EvtAbsLineShape( x )
{
    m_massMax = x.m_massMax;
    m_massMin = x.m_massMin;
    m_blattDecay = x.m_blattDecay;
    m_blattBirth = x.m_blattBirth;
    m_maxRange = x.m_maxRange;
    m_includeDecayFact = x.m_includeDecayFact;
    m_includeBirthFact = x.m_includeBirthFact;
    m_errorCond = x.m_errorCond;
}

EvtRelBreitWignerBarrierFact& EvtRelBreitWignerBarrierFact::operator=(
    const EvtRelBreitWignerBarrierFact& x )
{
    m_mass = x.m_mass;
    m_width = x.m_width;
    m_spin = x.m_spin;
    m_massMax = x.m_massMax;
    m_massMin = x.m_massMin;
    m_blattDecay = x.m_blattDecay;
    m_blattBirth = x.m_blattBirth;
    m_maxRange = x.m_maxRange;
    m_includeDecayFact = x.m_includeDecayFact;
    m_includeBirthFact = x.m_includeBirthFact;
    m_errorCond = x.m_errorCond;

    return *this;
}

EvtAbsLineShape* EvtRelBreitWignerBarrierFact::clone()
{
    return new EvtRelBreitWignerBarrierFact( *this );
}

double EvtRelBreitWignerBarrierFact::getMassProb( double mass, double massPar,
                                                  int nDaug, double* massDau )
{
    m_errorCond = false;
    //return EvtAbsLineShape::getMassProb(mass,massPar,nDaug,massDau);
    if ( nDaug != 2 )
        return EvtAbsLineShape::getMassProb( mass, massPar, nDaug, massDau );

    double dTotMass = 0.;

    int i;
    for ( i = 0; i < nDaug; i++ ) {
        dTotMass += massDau[i];
    }
    //EvtGenReport(EVTGEN_INFO,"EvtGen") << mass << " " << massPar << " " << dTotMass << " "<< endl;
    //    if ( (mass-dTotMass)<0.0001 ) return 0.;
    //EvtGenReport(EVTGEN_INFO,"EvtGen") << mass << " " << dTotMass << endl;
    if ( ( mass < dTotMass ) )
        return 0.;

    if ( m_width < 0.0001 )
        return 1.;

    if ( massPar > 0.0000000001 ) {
        if ( mass > massPar )
            return 0.;
    }

    if ( m_errorCond )
        return 0.;

    // we did all the work in getRandMass
    return 1.;
}

double EvtRelBreitWignerBarrierFact::getRandMass( EvtId* parId, int nDaug,
                                                  EvtId* dauId, EvtId* othDaugId,
                                                  double maxMass,
                                                  double* dauMasses )
{
    if ( nDaug != 2 )
        return EvtAbsLineShape::getRandMass( parId, nDaug, dauId, othDaugId,
                                             maxMass, dauMasses );

    if ( m_width < 0.00001 )
        return m_mass;

    //first figure out L - take the lowest allowed.

    EvtSpinType::spintype spinD1 = EvtPDL::getSpinType( dauId[0] );
    EvtSpinType::spintype spinD2 = EvtPDL::getSpinType( dauId[1] );

    int t1 = EvtSpinType::getSpin2( spinD1 );
    int t2 = EvtSpinType::getSpin2( spinD2 );
    int t3 = EvtSpinType::getSpin2( m_spin );

    int Lmin = -10;

    // the user has overridden the partial wave to use.
    for ( unsigned int vC = 0; vC < m_userSetPW.size(); vC++ ) {
        if ( dauId[0] == m_userSetPWD1[vC] && dauId[1] == m_userSetPWD2[vC] ) {
            Lmin = 2 * m_userSetPW[vC];
        }
        if ( dauId[0] == m_userSetPWD2[vC] && dauId[1] == m_userSetPWD1[vC] ) {
            Lmin = 2 * m_userSetPW[vC];
        }
    }

    // allow for special cases.
    if ( Lmin < -1 ) {
        //There are some things I don't know how to deal with
        if ( t3 > 4 )
            return EvtAbsLineShape::getRandMass( parId, nDaug, dauId, othDaugId,
                                                 maxMass, dauMasses );
        if ( t1 > 4 )
            return EvtAbsLineShape::getRandMass( parId, nDaug, dauId, othDaugId,
                                                 maxMass, dauMasses );
        if ( t2 > 4 )
            return EvtAbsLineShape::getRandMass( parId, nDaug, dauId, othDaugId,
                                                 maxMass, dauMasses );

        //figure the min and max allowwed "spins" for the daughters state
        Lmin = std::max( t3 - t2 - t1, std::max( t2 - t3 - t1, t1 - t3 - t2 ) );
        if ( Lmin < 0 )
            Lmin = 0;
        assert( Lmin == 0 || Lmin == 2 || Lmin == 4 );
    }

    //double massD1=EvtPDL::getMeanMass(dauId[0]);
    //double massD2=EvtPDL::getMeanMass(dauId[1]);
    double massD1 = dauMasses[0];
    double massD2 = dauMasses[1];

    // I'm not sure how to define the vertex factor here - so retreat to nonRel code.
    if ( ( massD1 + massD2 ) > m_mass )
        return EvtAbsLineShape::getRandMass( parId, nDaug, dauId, othDaugId,
                                             maxMass, dauMasses );

    //parent vertex factor not yet implemented
    double massOthD = -10.;
    double massParent = -10.;
    int birthl = -10;
    if ( othDaugId ) {
        EvtSpinType::spintype spinOth = EvtPDL::getSpinType( *othDaugId );
        EvtSpinType::spintype spinPar = EvtPDL::getSpinType( *parId );

        int tt1 = EvtSpinType::getSpin2( spinOth );
        int tt2 = EvtSpinType::getSpin2( spinPar );
        int tt3 = EvtSpinType::getSpin2( m_spin );

        //figure the min and max allowwed "spins" for the daughters state
        if ( ( tt1 <= 4 ) && ( tt2 <= 4 ) ) {
            birthl = std::max( tt3 - tt2 - tt1,
                               std::max( tt2 - tt3 - tt1, tt1 - tt3 - tt2 ) );
            if ( birthl < 0 )
                birthl = 0;

            massOthD = EvtPDL::getMeanMass( *othDaugId );
            massParent = EvtPDL::getMeanMass( *parId );
        }

        // allow user to override
        for ( size_t vC = 0; vC < m_userSetBirthPW.size(); vC++ ) {
            if ( *othDaugId == m_userSetBirthOthD[vC] &&
                 *parId == m_userSetBirthPar[vC] ) {
                birthl = 2 * m_userSetBirthPW[vC];
            }
        }
    }
    double massM = m_massMax;
    if ( ( maxMass > -0.5 ) && ( maxMass < massM ) )
        massM = maxMass;

    //special case... if the parent mass is _fixed_ we can do a little better
    //and only for a two body decay as that seems to be where we have problems

    // Define relativistic propagator amplitude

    EvtTwoBodyVertex vd( massD1, massD2, m_mass, Lmin / 2 );
    vd.set_f( m_blattDecay );
    EvtPropBreitWignerRel bw( m_mass, m_width );
    EvtMassAmp amp( bw, vd );

    if ( m_includeDecayFact ) {
        amp.addDeathFact();
        amp.addDeathFactFF();
    }
    if ( massParent > -1. ) {
        if ( m_includeBirthFact ) {
            EvtTwoBodyVertex vb( m_mass, massOthD, massParent, birthl / 2 );
            vb.set_f( m_blattBirth );
            amp.setBirthVtx( vb );
            amp.addBirthFact();
            amp.addBirthFactFF();
        }
    }

    EvtAmpPdf<EvtPoint1D> pdf( amp );

    // Estimate maximum and create predicate for accept reject

    double tempMaxLoc = m_mass;
    if ( maxMass > -0.5 && maxMass < m_mass )
        tempMaxLoc = maxMass;
    double tempMax = m_massMax;
    if ( maxMass > -0.5 && maxMass < m_massMax )
        tempMax = maxMass;
    double tempMinMass = m_massMin;
    if ( massD1 + massD2 > m_massMin )
        tempMinMass = massD1 + massD2;

    //redo sanity check - is there a solution to our problem.
    //if not return an error condition that is caught by the
    //mass prob calculation above.
    if ( tempMinMass > tempMax ) {
        m_errorCond = true;
        return tempMinMass;
    }

    if ( tempMaxLoc < tempMinMass )
        tempMaxLoc = tempMinMass;

    double safetyFactor = 1.2;

    EvtPdfMax<EvtPoint1D> max(
        safetyFactor *
        pdf.evaluate( EvtPoint1D( tempMinMass, tempMax, tempMaxLoc ) ) );

    EvtPdfPred<EvtPoint1D> pred( pdf );
    pred.setMax( max );

    EvtIntervalFlatPdf flat( tempMinMass, tempMax );
    EvtPdfGen<EvtPoint1D> gen( flat );
    EvtPredGen<EvtPdfGen<EvtPoint1D>, EvtPdfPred<EvtPoint1D>> predgen( gen, pred );

    EvtPoint1D point = predgen();
    return point.value();
}
