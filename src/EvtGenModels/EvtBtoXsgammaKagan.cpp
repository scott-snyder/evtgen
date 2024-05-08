
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

#include "EvtGenModels/EvtBtoXsgammaKagan.hh"

#include "EvtGenBase/EvtConst.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtReport.hh"

#include "EvtGenModels/EvtBtoXsgamma.hh"
#include "EvtGenModels/EvtBtoXsgammaFermiUtil.hh"
#include "EvtGenModels/EvtItgAbsIntegrator.hh"
#include "EvtGenModels/EvtItgFourCoeffFcn.hh"
#include "EvtGenModels/EvtItgFunction.hh"
#include "EvtGenModels/EvtItgPtrFunction.hh"
#include "EvtGenModels/EvtItgSimpsonIntegrator.hh"
#include "EvtGenModels/EvtItgThreeCoeffFcn.hh"
#include "EvtGenModels/EvtItgTwoCoeffFcn.hh"

#include <fstream>
#include <stdlib.h>
#include <string>
using std::endl;
using std::fstream;

bool EvtBtoXsgammaKagan::m_bbprod = false;
double EvtBtoXsgammaKagan::m_intervalMH = 0;

void EvtBtoXsgammaKagan::init( int nArg, double* args )
{
    if ( ( nArg ) > 12 || ( nArg > 1 && nArg < 10 ) || nArg == 11 ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "EvtBtoXsgamma generator model "
            << "EvtBtoXsgammaKagan expected "
            << "either 1(default config) or "
            << "10 (default mass range) or "
            << "12 (user range) arguments but found: " << nArg << endl;
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Will terminate execution!" << endl;
        ::abort();
    }

    if ( nArg == 1 ) {
        m_bbprod = true;
        getDefaultHadronicMass();
    } else {
        m_bbprod = false;
        computeHadronicMass( nArg, args );
    }

    double mHminLimit = 0.6373;
    double mHmaxLimit = 4.5;

    if ( nArg > 10 ) {
        m_mHmin = args[10];
        m_mHmax = args[11];
        if ( m_mHmin > m_mHmax ) {
            EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                << "Minimum hadronic mass exceeds maximum " << endl;
            EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                << "Will terminate execution!" << endl;
            ::abort();
        }
        if ( m_mHmin < mHminLimit ) {
            EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                << "Minimum hadronic mass below K pi threshold" << endl;
            EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                << "Resetting to K pi threshold" << endl;
            m_mHmin = mHminLimit;
        }
        if ( m_mHmax > mHmaxLimit ) {
            EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                << "Maximum hadronic mass above 4.5 GeV/c^2" << endl;
            EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                << "Resetting to 4.5 GeV/c^2" << endl;
            m_mHmax = mHmaxLimit;
        }
    } else {
        m_mHmin = mHminLimit;    //  usually just above K pi threshold for Xsd/u
        m_mHmax = mHmaxLimit;
    }
}

void EvtBtoXsgammaKagan::getDefaultHadronicMass()
{
    m_massHad = { 0,        0.0625995, 0.125199, 0.187798, 0.250398, 0.312997,
                  0.375597, 0.438196,  0.500796, 0.563395, 0.625995, 0.688594,
                  0.751194, 0.813793,  0.876392, 0.938992, 1.00159,  1.06419,
                  1.12679,  1.18939,   1.25199,  1.31459,  1.37719,  1.43979,
                  1.50239,  1.56499,   1.62759,  1.69019,  1.75278,  1.81538,
                  1.87798,  1.94058,   2.00318,  2.06578,  2.12838,  2.19098,
                  2.25358,  2.31618,   2.37878,  2.44138,  2.50398,  2.56658,
                  2.62918,  2.69178,   2.75438,  2.81698,  2.87958,  2.94217,
                  3.00477,  3.06737,   3.12997,  3.19257,  3.25517,  3.31777,
                  3.38037,  3.44297,   3.50557,  3.56817,  3.63077,  3.69337,
                  3.75597,  3.81857,   3.88117,  3.94377,  4.00637,  4.06896,
                  4.13156,  4.19416,   4.25676,  4.31936,  4.38196,  4.44456,
                  4.50716,  4.56976,   4.63236,  4.69496,  4.75756,  4.82016,
                  4.88276,  4.94536,   5.00796 };
    m_brHad = { 0,           1.03244e-09, 3.0239e-08,  1.99815e-07, 7.29392e-07,
                1.93129e-06, 4.17806e-06, 7.86021e-06, 1.33421e-05, 2.09196e-05,
                3.07815e-05, 4.29854e-05, 5.74406e-05, 7.3906e-05,  9.2003e-05,
                0.000111223, 0.000130977, 0.000150618, 0.000169483, 0.000186934,
                0.000202392, 0.000215366, 0.000225491, 0.000232496, 0.000236274,
                0.000236835, 0.000234313, 0.000228942, 0.000221042, 0.000210994,
                0.000199215, 0.000186137, 0.000172194, 0.000157775, 0.000143255,
                0.000128952, 0.000115133, 0.000102012, 8.97451e-05, 7.84384e-05,
                6.81519e-05, 5.89048e-05, 5.06851e-05, 4.34515e-05, 3.71506e-05,
                3.1702e-05,  2.70124e-05, 2.30588e-05, 1.96951e-05, 1.68596e-05,
                1.44909e-05, 1.25102e-05, 1.08596e-05, 9.48476e-06, 8.34013e-06,
                7.38477e-06, 6.58627e-06, 5.91541e-06, 5.35022e-06, 4.87047e-06,
                4.46249e-06, 4.11032e-06, 3.80543e-06, 3.54051e-06, 3.30967e-06,
                3.10848e-06, 2.93254e-06, 2.78369e-06, 2.65823e-06, 2.55747e-06,
                2.51068e-06, 2.57179e-06, 2.74684e-06, 3.02719e-06, 3.41182e-06,
                3.91387e-06, 4.56248e-06, 5.40862e-06, 6.53915e-06, 8.10867e-06,
                1.04167e-05 };
    m_massHad.resize( 81 );
    m_brHad.resize( 81 );

    m_intervalMH = 80;
}

void EvtBtoXsgammaKagan::computeHadronicMass( int /*nArg*/, double* args )
{
    //Input parameters
    int fermiFunction = (int)args[1];
    m_mB = args[2];
    m_mb = args[3];
    m_mu = args[4];
    m_lam1 = args[5];
    m_delta = args[6];
    m_z = args[7];
    m_nIntervalS = args[8];
    m_nIntervalmH = args[9];
    std::vector<double> mHVect( int( m_nIntervalmH + 1.0 ) );
    m_massHad.clear();
    m_massHad.resize( int( m_nIntervalmH + 1.0 ) );
    m_brHad.clear();
    m_brHad.resize( int( m_nIntervalmH + 1.0 ) );
    m_intervalMH = m_nIntervalmH;

    //Going to have to add a new entry into the data file - takes ages...
    EvtGenReport( EVTGEN_WARNING, "EvtGen" )
        << "EvtBtoXsgammaKagan: calculating new hadronic mass spectra. This takes a while..."
        << endl;

    //Now need to compute the mHVect vector for
    //the current parameters

    //A few more parameters
    double _mubar = m_mu;
    m_mW = 80.33;
    m_mt = 175.0;
    m_alpha = 1. / 137.036;
    m_lambdabar = m_mB - m_mb;
    m_kappabar = 3.382 - 4.14 * ( sqrt( m_z ) - 0.29 );
    m_fz = Fz( m_z );
    m_rer8 = ( 44. / 9. ) - ( 8. / 27. ) * pow( EvtConst::pi, 2. );
    m_r7 = ( -10. / 3. ) - ( 8. / 9. ) * pow( EvtConst::pi, 2. );
    m_rer2 = -4.092 + 12.78 * ( sqrt( m_z ) - .29 );
    m_gam77 = 32. / 3.;
    m_gam27 = 416. / 81.;
    m_gam87 = -32. / 9.;
    m_lam2 = .12;
    m_beta0 = 23. / 3.;
    m_beta1 = 116. / 3.;
    m_alphasmZ = .118;
    m_mZ = 91.187;
    m_ms = m_mb / 50.;

    double eGammaMin = 0.5 * m_mB * ( 1. - m_delta );
    double eGammaMax = 0.5 * m_mB;
    double yMin = 2. * eGammaMin / m_mB;
    double yMax = 2. * eGammaMax / m_mB;
    double _CKMrat = 0.976;
    double Nsl = 1.0;

    //Calculate alpha the various scales
    m_alphasmW = CalcAlphaS( m_mW );
    m_alphasmt = CalcAlphaS( m_mt );
    m_alphasmu = CalcAlphaS( m_mu );
    m_alphasmubar = CalcAlphaS( _mubar );

    //Calculate the Wilson Coefficients and Delta
    m_etamu = m_alphasmW / m_alphasmu;
    m_kSLemmu = ( 12. / 23. ) * ( ( 1. / m_etamu ) - 1. );
    CalcWilsonCoeffs();
    CalcDelta();

    //Build s22 and s27 vector - saves time because double
    //integration is required otherwise
    std::vector<double> s22Coeffs( int( m_nIntervalS + 1.0 ) );
    std::vector<double> s27Coeffs( int( m_nIntervalS + 1.0 ) );
    std::vector<double> s28Coeffs( int( m_nIntervalS + 1.0 ) );

    double dy = ( yMax - yMin ) / m_nIntervalS;
    double yp = yMin;

    std::vector<double> sCoeffs( 1 );
    sCoeffs[0] = m_z;

    //Define s22 and s27 functions
    auto mys22Func = EvtItgPtrFunction{ &s22Func, 0., yMax + 0.1, sCoeffs };
    auto mys27Func = EvtItgPtrFunction{ &s27Func, 0., yMax + 0.1, sCoeffs };

    //Use a simpson integrator
    auto mys22Simp = EvtItgSimpsonIntegrator{ mys22Func, 1.0e-4, 20 };
    auto mys27Simp = EvtItgSimpsonIntegrator{ mys27Func, 1.0e-4, 50 };

    int i;

    for ( i = 0; i < int( m_nIntervalS + 1.0 ); i++ ) {
        s22Coeffs[i] = ( 16. / 27. ) * mys22Simp.evaluate( 1.0e-20, yp );
        s27Coeffs[i] = ( -8. / 9. ) * m_z * mys27Simp.evaluate( 1.0e-20, yp );
        s28Coeffs[i] = -s27Coeffs[i] / 3.;
        yp = yp + dy;
    }

    //Define functions and vectors used to calculate mHVect. Each function takes a set
    //of vectors which are used as the function coefficients
    std::vector<double> FermiCoeffs( 6 );
    std::vector<double> varCoeffs( 3 );
    std::vector<double> DeltaCoeffs( 1 );
    std::vector<double> s88Coeffs( 2 );
    std::vector<double> sInitCoeffs( 3 );

    varCoeffs[0] = m_mB;
    varCoeffs[1] = m_mb;
    varCoeffs[2] = 0.;

    DeltaCoeffs[0] = m_alphasmu;

    s88Coeffs[0] = m_mb;
    s88Coeffs[1] = m_ms;

    sInitCoeffs[0] = m_nIntervalS;
    sInitCoeffs[1] = yMin;
    sInitCoeffs[2] = yMax;

    FermiCoeffs[0] = fermiFunction;
    FermiCoeffs[1] = 0.0;
    FermiCoeffs[2] = 0.0;
    FermiCoeffs[3] = 0.0;
    FermiCoeffs[4] = 0.0;
    FermiCoeffs[5] = 0.0;

    //Coefficients for gamma function
    std::vector<double> gammaCoeffs( 6 );
    gammaCoeffs[0] = 76.18009172947146;
    gammaCoeffs[1] = -86.50532032941677;
    gammaCoeffs[2] = 24.01409824083091;
    gammaCoeffs[3] = -1.231739572450155;
    gammaCoeffs[4] = 0.1208650973866179e-2;
    gammaCoeffs[5] = -0.5395239384953e-5;

    //Calculate quantities for the fermi function to be used
    //Distinguish among the different shape functions
    if ( fermiFunction == 1 ) {
        FermiCoeffs[1] = m_lambdabar;
        FermiCoeffs[2] = ( -3. * pow( m_lambdabar, 2. ) / m_lam1 ) - 1.;
        FermiCoeffs[3] = m_lam1;
        FermiCoeffs[4] = 1.0;

        auto myNormFunc = std::make_unique<EvtItgPtrFunction>(
            &EvtBtoXsgammaFermiUtil::FermiExpFunc, -m_mb, m_mB - m_mb,
            FermiCoeffs );
        auto myNormSimp =
            std::make_unique<EvtItgSimpsonIntegrator>( *myNormFunc, 1.0e-4, 40 );
        FermiCoeffs[4] = myNormSimp->normalisation();

    } else if ( fermiFunction == 2 ) {
        double a = EvtBtoXsgammaFermiUtil::FermiGaussFuncRoot( m_lambdabar,
                                                               m_lam1, m_mb,
                                                               gammaCoeffs );
        FermiCoeffs[1] = m_lambdabar;
        FermiCoeffs[2] = a;
        FermiCoeffs[3] =
            EvtBtoXsgammaFermiUtil::Gamma( ( 2.0 + a ) / 2., gammaCoeffs ) /
            EvtBtoXsgammaFermiUtil::Gamma( ( 1.0 + a ) / 2., gammaCoeffs );
        FermiCoeffs[4] = 1.0;

        auto myNormFunc = std::make_unique<EvtItgPtrFunction>(
            &EvtBtoXsgammaFermiUtil::FermiGaussFunc, -m_mb, m_mB - m_mb,
            FermiCoeffs );
        auto myNormSimp =
            std::make_unique<EvtItgSimpsonIntegrator>( *myNormFunc, 1.0e-4, 40 );
        FermiCoeffs[4] = myNormSimp->normalisation();

    } else if ( fermiFunction == 3 ) {
        double rho = EvtBtoXsgammaFermiUtil::FermiRomanFuncRoot( m_lambdabar,
                                                                 m_lam1 );
        FermiCoeffs[1] = m_mB;
        FermiCoeffs[2] = m_mb;
        FermiCoeffs[3] = rho;
        FermiCoeffs[4] = m_lambdabar;
        FermiCoeffs[5] = 1.0;

        auto myNormFunc = std::make_unique<EvtItgPtrFunction>(
            &EvtBtoXsgammaFermiUtil::FermiRomanFunc, -m_mb, m_mB - m_mb,
            FermiCoeffs );
        auto myNormSimp =
            std::make_unique<EvtItgSimpsonIntegrator>( *myNormFunc, 1.0e-4, 40 );
        FermiCoeffs[5] = myNormSimp->normalisation();
    }

    //Define functions
    auto myDeltaFermiFunc = EvtItgThreeCoeffFcn{ &DeltaFermiFunc, -m_mb,
                                                 m_mB - m_mb,     FermiCoeffs,
                                                 varCoeffs,       DeltaCoeffs };
    auto mys88FermiFunc = EvtItgThreeCoeffFcn{ &s88FermiFunc, -m_mb,
                                               m_mB - m_mb,   FermiCoeffs,
                                               varCoeffs,     s88Coeffs };
    auto mys77FermiFunc = EvtItgTwoCoeffFcn{ &s77FermiFunc, -m_mb, m_mB - m_mb,
                                             FermiCoeffs, varCoeffs };
    auto mys78FermiFunc = EvtItgTwoCoeffFcn{ &s78FermiFunc, -m_mb, m_mB - m_mb,
                                             FermiCoeffs, varCoeffs };
    auto mys22FermiFunc = EvtItgFourCoeffFcn{ &sFermiFunc, -m_mb,
                                              m_mB - m_mb, FermiCoeffs,
                                              varCoeffs,   sInitCoeffs,
                                              s22Coeffs };
    auto mys27FermiFunc = EvtItgFourCoeffFcn{ &sFermiFunc, -m_mb,
                                              m_mB - m_mb, FermiCoeffs,
                                              varCoeffs,   sInitCoeffs,
                                              s27Coeffs };
    auto mys28FermiFunc = EvtItgFourCoeffFcn{ &sFermiFunc, -m_mb,
                                              m_mB - m_mb, FermiCoeffs,
                                              varCoeffs,   sInitCoeffs,
                                              s28Coeffs };

    //Define integrators
    auto myDeltaFermiSimp = EvtItgSimpsonIntegrator{ myDeltaFermiFunc, 1.0e-4,
                                                     40 };
    auto mys77FermiSimp = EvtItgSimpsonIntegrator{ mys77FermiFunc, 1.0e-4, 40 };
    auto mys88FermiSimp = EvtItgSimpsonIntegrator{ mys88FermiFunc, 1.0e-4, 40 };
    auto mys78FermiSimp = EvtItgSimpsonIntegrator{ mys78FermiFunc, 1.0e-4, 40 };
    auto mys22FermiSimp = EvtItgSimpsonIntegrator{ mys22FermiFunc, 1.0e-4, 40 };
    auto mys27FermiSimp = EvtItgSimpsonIntegrator{ mys27FermiFunc, 1.0e-4, 40 };
    auto mys28FermiSimp = EvtItgSimpsonIntegrator{ mys28FermiFunc, 1.0e-4, 40 };

    //Finally calculate mHVect for the range of hadronic masses
    double mHmin = sqrt( m_mB * m_mB - 2. * m_mB * eGammaMax );
    double mHmax = sqrt( m_mB * m_mB - 2. * m_mB * eGammaMin );
    double dmH = ( mHmax - mHmin ) / m_nIntervalmH;

    double mH = mHmin;

    //Calculating the Branching Fractions
    for ( i = 0; i < int( m_nIntervalmH + 1.0 ); i++ ) {
        double ymH = 1. - ( ( mH * mH ) / ( m_mB * m_mB ) );

        //Need to set ymH as one of the input parameters
        myDeltaFermiFunc.setCoeff( 2, 2, ymH );
        mys77FermiFunc.setCoeff( 2, 2, ymH );
        mys88FermiFunc.setCoeff( 2, 2, ymH );
        mys78FermiFunc.setCoeff( 2, 2, ymH );
        mys22FermiFunc.setCoeff( 2, 2, ymH );
        mys27FermiFunc.setCoeff( 2, 2, ymH );
        mys28FermiFunc.setCoeff( 2, 2, ymH );

        //Integrate

        double deltaResult = myDeltaFermiSimp.evaluate( ( m_mB * ymH - m_mb ),
                                                        m_mB - m_mb );
        double s77Result = mys77FermiSimp.evaluate( ( m_mB * ymH - m_mb ),
                                                    m_mB - m_mb );
        double s88Result = mys88FermiSimp.evaluate( ( m_mB * ymH - m_mb ),
                                                    m_mB - m_mb );
        double s78Result = mys78FermiSimp.evaluate( ( m_mB * ymH - m_mb ),
                                                    m_mB - m_mb );
        double s22Result = mys22FermiSimp.evaluate( ( m_mB * ymH - m_mb ),
                                                    m_mB - m_mb );
        double s27Result = mys27FermiSimp.evaluate( ( m_mB * ymH - m_mb ),
                                                    m_mB - m_mb );
        mys28FermiSimp.evaluate( ( m_mB * ymH - m_mb ), m_mB - m_mb );

        double py =
            ( pow( _CKMrat, 2. ) * ( 6. / m_fz ) * ( m_alpha / EvtConst::pi ) *
              ( deltaResult * m_cDeltatot +
                ( m_alphasmu / EvtConst::pi ) *
                    ( s77Result * pow( m_c70mu, 2. ) +
                      s27Result * m_c2mu * ( m_c70mu - m_c80mu / 3. ) +
                      s78Result * m_c70mu * m_c80mu + s22Result * m_c2mu * m_c2mu +
                      s88Result * m_c80mu * m_c80mu ) ) );

        mHVect[i] = 2. * ( mH / ( m_mB * m_mB ) ) * 0.105 * Nsl * py;

        m_massHad[i] = mH;
        m_brHad[i] = 2. * ( mH / ( m_mB * m_mB ) ) * 0.105 * Nsl * py;

        mH = mH + dmH;
    }
}

double EvtBtoXsgammaKagan::GetMass( int /*Xscode*/ )
{
    //  Get hadronic mass for the event according to the hadronic mass spectra computed in computeHadronicMass
    double mass = 0.0;
    double min = m_mHmin;
    if ( m_bbprod )
        min = 1.1;
    //  double max=4.5;
    double max = m_mHmax;
    double xbox( 0 ), ybox( 0 );
    double boxheight( 0 );
    double trueHeight( 0 );
    double boxwidth = max - min;
    double wgt( 0. );

    for ( int i = 0; i < int( m_intervalMH + 1.0 ); i++ ) {
        if ( m_brHad[i] > boxheight )
            boxheight = m_brHad[i];
    }
    while ( ( mass > max ) || ( mass < min ) ) {
        xbox = EvtRandom::Flat( boxwidth ) + min;
        ybox = EvtRandom::Flat( boxheight );
        trueHeight = 0.0;
        // Correction by Peter Richardson
        for ( int i = 1; i < int( m_intervalMH + 1.0 ); ++i ) {
            if ( ( m_massHad[i] >= xbox ) && ( 0.0 == trueHeight ) ) {
                wgt = ( xbox - m_massHad[i - 1] ) /
                      ( m_massHad[i] - m_massHad[i - 1] );
                trueHeight = m_brHad[i - 1] +
                             wgt * ( m_brHad[i] - m_brHad[i - 1] );
            }
        }

        if ( ybox > trueHeight ) {
            mass = 0.0;
        } else {
            mass = xbox;
        }
    }

    return mass;
}

double EvtBtoXsgammaKagan::CalcAlphaS( double scale )
{
    double v = 1. - m_beta0 * ( m_alphasmZ / ( 2. * EvtConst::pi ) ) *
                        ( log( m_mZ / scale ) );
    return ( m_alphasmZ / v ) *
           ( 1. - ( ( m_beta1 / m_beta0 ) *
                    ( m_alphasmZ / ( 4. * EvtConst::pi ) ) * ( log( v ) / v ) ) );
}

void EvtBtoXsgammaKagan::CalcWilsonCoeffs()
{
    double mtatmw = m_mt * pow( ( m_alphasmW / m_alphasmt ), ( 12. / 23. ) ) *
                    ( 1 +
                      ( 12. / 23. ) * ( ( 253. / 18. ) - ( 116. / 23. ) ) *
                          ( ( m_alphasmW - m_alphasmt ) / ( 4.0 * EvtConst::pi ) ) -
                      ( 4. / 3. ) * ( m_alphasmt / EvtConst::pi ) );
    double xt = pow( mtatmw, 2. ) / pow( m_mW, 2. );

    /////LO
    m_c2mu = .5 * pow( m_etamu, ( -12. / 23. ) ) +
             .5 * pow( m_etamu, ( 6. / 23. ) );

    double c7mWsm = ( ( 3. * pow( xt, 3. ) - 2. * pow( xt, 2. ) ) /
                      ( 4. * pow( ( xt - 1. ), 4. ) ) ) *
                        log( xt ) +
                    ( ( -8. * pow( xt, 3. ) - 5. * pow( xt, 2. ) + 7. * xt ) /
                      ( 24. * pow( ( xt - 1. ), 3. ) ) );

    double c8mWsm = ( ( -3. * pow( xt, 2. ) ) / ( 4. * pow( ( xt - 1. ), 4. ) ) ) *
                        log( xt ) +
                    ( ( -pow( xt, 3. ) + 5. * pow( xt, 2. ) + 2. * xt ) /
                      ( 8. * pow( ( xt - 1. ), 3. ) ) );

    double c7constmu = ( 626126. / 272277. ) * pow( m_etamu, ( 14. / 23. ) ) -
                       ( 56281. / 51730. ) * pow( m_etamu, ( 16. / 23. ) ) -
                       ( 3. / 7. ) * pow( m_etamu, ( 6. / 23. ) ) -
                       ( 1. / 14. ) * pow( m_etamu, ( -12. / 23. ) ) -
                       .6494 * pow( m_etamu, .4086 ) -
                       .038 * pow( m_etamu, -.423 ) -
                       .0186 * pow( m_etamu, -.8994 ) -
                       .0057 * pow( m_etamu, .1456 );

    m_c70mu = c7mWsm * pow( m_etamu, ( 16. / 23. ) ) +
              ( 8. / 3. ) *
                  ( pow( m_etamu, ( 14. / 23. ) ) -
                    pow( m_etamu, ( 16. / 23. ) ) ) *
                  c8mWsm +
              c7constmu;

    double c8constmu = ( 313063. / 363036. ) * pow( m_etamu, ( 14. / 23. ) ) -
                       .9135 * pow( m_etamu, .4086 ) +
                       .0873 * pow( m_etamu, -.423 ) -
                       .0571 * pow( m_etamu, -.8994 ) +
                       .0209 * pow( m_etamu, .1456 );

    m_c80mu = c8mWsm * pow( m_etamu, ( 14. / 23. ) ) + c8constmu;

    //Compute the dilogarithm (PolyLog(2,x)) with the Simpson integrator
    //The dilogarithm is defined as: Li_2(x)=Int_0^x(-log(1.-z)/z)
    //however, Mathematica implements it as  Sum[z^k/k^2,{k,1,Infinity}], so, althought the two
    //results are similar and both implemented in the program, we prefer to use the
    //one closer to the Mathematica implementation as identical to what used by the theorists.

    // EvtItgFunction *myDiLogFunc = new EvtItgFunction(&diLogFunc, 0., 1.-1./xt);
    //EvtItgAbsIntegrator *myDiLogSimp = new EvtItgSimpsonIntegrator(*myDiLogFunc, 1.0e-4, 50);
    //double li2 = myDiLogSimp->evaluate(1.0e-20,1.-1./xt);

    double li2 = diLogMathematica( 1. - 1. / xt );

    double c7mWsm1 =
        ( ( -16. * pow( xt, 4. ) - 122. * pow( xt, 3. ) + 80. * pow( xt, 2. ) -
            8. * xt ) /
              ( 9. * pow( ( xt - 1. ), 4. ) ) * li2 +
          ( 6. * pow( xt, 4. ) + 46. * pow( xt, 3. ) - 28. * pow( xt, 2. ) ) /
              ( 3. * pow( ( xt - 1. ), 5. ) ) * pow( log( xt ), 2. ) +
          ( -102. * pow( xt, 5. ) - 588. * pow( xt, 4. ) -
            2262. * pow( xt, 3. ) + 3244. * pow( xt, 2. ) - 1364. * xt + 208. ) /
              ( 81. * pow( ( xt - 1 ), 5. ) ) * log( xt ) +
          ( 1646. * pow( xt, 4. ) + 12205. * pow( xt, 3. ) -
            10740. * pow( xt, 2. ) + 2509. * xt - 436. ) /
              ( 486. * pow( ( xt - 1 ), 4. ) ) );

    double c8mWsm1 =
        ( ( -4. * pow( xt, 4. ) + 40. * pow( xt, 3. ) + 41. * pow( xt, 2. ) + xt ) /
              ( 6. * pow( ( xt - 1. ), 4. ) ) * li2 +
          ( -17. * pow( xt, 3. ) - 31. * pow( xt, 2. ) ) /
              ( 2. * pow( ( xt - 1. ), 5. ) ) * pow( log( xt ), 2. ) +
          ( -210. * pow( xt, 5. ) + 1086. * pow( xt, 4. ) +
            4893. * pow( xt, 3. ) + 2857. * pow( xt, 2. ) - 1994. * xt + 280. ) /
              ( 216. * pow( ( xt - 1 ), 5. ) ) * log( xt ) +
          ( 737. * pow( xt, 4. ) - 14102. * pow( xt, 3. ) -
            28209. * pow( xt, 2. ) + 610. * xt - 508. ) /
              ( 1296. * pow( ( xt - 1 ), 4. ) ) );

    double E1 = ( xt * ( 18. - 11. * xt - pow( xt, 2. ) ) /
                      ( 12. * pow( ( 1. - xt ), 3. ) ) +
                  pow( xt, 2. ) * ( 15. - 16. * xt + 4. * pow( xt, 2. ) ) /
                      ( 6. * pow( ( 1. - xt ), 4. ) ) * log( xt ) -
                  2. / 3. * log( xt ) );

    double e1 = 4661194. / 816831.;
    double e2 = -8516. / 2217.;
    double e3 = 0.;
    double e4 = 0.;
    double e5 = -1.9043;
    double e6 = -.1008;
    double e7 = .1216;
    double e8 = .0183;

    double f1 = -17.3023;
    double f2 = 8.5027;
    double f3 = 4.5508;
    double f4 = .7519;
    double f5 = 2.004;
    double f6 = .7476;
    double f7 = -.5385;
    double f8 = .0914;

    double g1 = 14.8088;
    double g2 = -10.809;
    double g3 = -.874;
    double g4 = .4218;
    double g5 = -2.9347;
    double g6 = .3971;
    double g7 = .1600;
    double g8 = .0225;

    double c71constmu =
        ( ( e1 * m_etamu * E1 + f1 + g1 * m_etamu ) *
              pow( m_etamu, ( 14. / 23. ) ) +
          ( e2 * m_etamu * E1 + f2 + g2 * m_etamu ) *
              pow( m_etamu, ( 16. / 23. ) ) +
          ( e3 * m_etamu * E1 + f3 + g3 * m_etamu ) * pow( m_etamu, ( 6. / 23. ) ) +
          ( e4 * m_etamu * E1 + f4 + g4 * m_etamu ) *
              pow( m_etamu, ( -12. / 23. ) ) +
          ( e5 * m_etamu * E1 + f5 + g5 * m_etamu ) * pow( m_etamu, .4086 ) +
          ( e6 * m_etamu * E1 + f6 + g6 * m_etamu ) * pow( m_etamu, ( -.423 ) ) +
          ( e7 * m_etamu * E1 + f7 + g7 * m_etamu ) * pow( m_etamu, ( -.8994 ) ) +
          ( e8 * m_etamu * E1 + f8 + g8 * m_etamu ) * pow( m_etamu, .1456 ) );

    double c71pmu = ( ( ( 297664. / 14283. * pow( m_etamu, ( 16. / 23. ) ) -
                          7164416. / 357075. * pow( m_etamu, ( 14. / 23. ) ) +
                          256868. / 14283. * pow( m_etamu, ( 37. / 23. ) ) -
                          6698884. / 357075. * pow( m_etamu, ( 39. / 23. ) ) ) *
                        ( c8mWsm ) ) +
                      37208. / 4761. *
                          ( pow( m_etamu, ( 39. / 23. ) ) -
                            pow( m_etamu, ( 16. / 23. ) ) ) *
                          ( c7mWsm ) +
                      c71constmu );

    m_c71mu = ( m_alphasmW / m_alphasmu *
                    ( pow( m_etamu, ( 16. / 23. ) ) * c7mWsm1 +
                      8. / 3. *
                          ( pow( m_etamu, ( 14. / 23. ) ) -
                            pow( m_etamu, ( 16. / 23. ) ) ) *
                          c8mWsm1 ) +
                c71pmu );

    m_c7emmu = ( ( 32. / 75. * pow( m_etamu, ( -9. / 23. ) ) -
                   40. / 69. * pow( m_etamu, ( -7. / 23. ) ) +
                   88. / 575. * pow( m_etamu, ( 16. / 23. ) ) ) *
                     c7mWsm +
                 ( -32. / 575. * pow( m_etamu, ( -9. / 23. ) ) +
                   32. / 1449. * pow( m_etamu, ( -7. / 23. ) ) +
                   640. / 1449. * pow( m_etamu, ( 14. / 23. ) ) -
                   704. / 1725. * pow( m_etamu, ( 16. / 23. ) ) ) *
                     c8mWsm -
                 190. / 8073. * pow( m_etamu, ( -35. / 23. ) ) -
                 359. / 3105. * pow( m_etamu, ( -17. / 23. ) ) +
                 4276. / 121095. * pow( m_etamu, ( -12. / 23. ) ) +
                 350531. / 1009125. * pow( m_etamu, ( -9. / 23. ) ) +
                 2. / 4347. * pow( m_etamu, ( -7. / 23. ) ) -
                 5956. / 15525. * pow( m_etamu, ( 6. / 23. ) ) +
                 38380. / 169533. * pow( m_etamu, ( 14. / 23. ) ) -
                 748. / 8625. * pow( m_etamu, ( 16. / 23. ) ) );

    // Wilson coefficients values as according to Kagan's program
    // m_c2mu=1.10566;
    //m_c70mu=-0.314292;
    // m_c80mu=-0.148954;
    // m_c71mu=0.480964;
    // m_c7emmu=0.0323219;
}

void EvtBtoXsgammaKagan::CalcDelta()
{
    double cDelta77 = ( 1. +
                        ( m_alphasmu / ( 2. * EvtConst::pi ) ) *
                            ( m_r7 - ( 16. / 3. ) + m_gam77 * log( m_mb / m_mu ) ) +
                        ( ( pow( ( 1. - m_z ), 4. ) / m_fz ) - 1. ) *
                            ( 6. * m_lam2 / pow( m_mb, 2. ) ) +
                        ( m_alphasmubar / ( 2. * EvtConst::pi ) ) * m_kappabar ) *
                      pow( m_c70mu, 2. );

    double cDelta27 = ( ( m_alphasmu / ( 2. * EvtConst::pi ) ) *
                            ( m_rer2 + m_gam27 * log( m_mb / m_mu ) ) -
                        ( m_lam2 / ( 9. * m_z * pow( m_mb, 2. ) ) ) ) *
                      m_c2mu * m_c70mu;

    double cDelta78 = ( m_alphasmu / ( 2. * EvtConst::pi ) ) *
                      ( m_rer8 + m_gam87 * log( m_mb / m_mu ) ) * m_c70mu *
                      m_c80mu;

    m_cDeltatot = cDelta77 + cDelta27 + cDelta78 +
                  ( m_alphasmu / ( 2. * EvtConst::pi ) ) * m_c71mu * m_c70mu +
                  ( m_alpha / m_alphasmu ) * ( 2. * m_c7emmu * m_c70mu -
                                               m_kSLemmu * pow( m_c70mu, 2. ) );
}

double EvtBtoXsgammaKagan::Delta( double y, double alphasMu )
{
    //Fix for singularity at endpoint
    if ( y >= 1.0 )
        y = 0.9999999999;

    return ( -4. * ( alphasMu / ( 3. * EvtConst::pi * ( 1. - y ) ) ) *
             ( log( 1. - y ) + 7. / 4. ) *
             exp( -2. * ( alphasMu / ( 3. * EvtConst::pi ) ) *
                  ( pow( log( 1. - y ), 2 ) + ( 7. / 2. ) * log( 1. - y ) ) ) );
}

double EvtBtoXsgammaKagan::s77( double y )
{
    //Fix for singularity at endpoint
    if ( y >= 1.0 )
        y = 0.9999999999;

    return ( ( 1. / 3. ) *
             ( 7. + y - 2. * pow( y, 2 ) - 2. * ( 1. + y ) * log( 1. - y ) ) );
}

double EvtBtoXsgammaKagan::s88( double y, double mb, double ms )
{
    //Fix for singularity at endpoint
    if ( y >= 1.0 )
        y = 0.9999999999;

    return ( ( 1. / 27. ) * ( ( 2. * ( 2. - 2. * y + pow( y, 2 ) ) / y ) *
                                  ( log( 1. - y ) + 2. * log( mb / ms ) ) -
                              2. * pow( y, 2 ) - y - 8. * ( ( 1. - y ) / y ) ) );
}

double EvtBtoXsgammaKagan::s78( double y )
{
    //Fix for singularity at endpoint
    if ( y >= 1.0 )
        y = 0.9999999999;

    return ( ( 8. / 9. ) * ( ( ( 1. - y ) / y ) * log( 1. - y ) + 1. +
                             ( pow( y, 2 ) / 4. ) ) );
}

double EvtBtoXsgammaKagan::ReG( double y )
{
    if ( y < 4. )
        return -2. * pow( atan( sqrt( y / ( 4. - y ) ) ), 2. );
    else {
        return 2. * ( pow( log( ( sqrt( y ) + sqrt( y - 4. ) ) / 2. ), 2. ) ) -
               ( 1. / 2. ) * pow( EvtConst::pi, 2. );
    }
}

double EvtBtoXsgammaKagan::ImG( double y )
{
    if ( y < 4. )
        return 0.0;
    else {
        return ( -2. * EvtConst::pi * log( ( sqrt( y ) + sqrt( y - 4. ) ) / 2. ) );
    }
}

double EvtBtoXsgammaKagan::s22Func( double y, const std::vector<double>& coeffs )
{
    //coeffs[0]=z
    return ( 1. - y ) *
           ( ( pow( coeffs[0], 2. ) / pow( y, 2. ) ) *
                 ( pow( ReG( y / coeffs[0] ), 2. ) +
                   pow( ImG( y / coeffs[0] ), 2. ) ) +
             ( coeffs[0] / y ) * ReG( y / coeffs[0] ) + ( 1. / 4. ) );
}

double EvtBtoXsgammaKagan::s27Func( double y, const std::vector<double>& coeffs )
{
    //coeffs[0] = z
    return ( ReG( y / coeffs[0] ) + y / ( 2. * coeffs[0] ) );
}

double EvtBtoXsgammaKagan::DeltaFermiFunc( double y,
                                           const std::vector<double>& coeffs1,
                                           const std::vector<double>& coeffs2,
                                           const std::vector<double>& coeffs3 )
{
    //coeffs1=fermi function coeffs, coeffs2[0]=mB, coeffs2[1]=mb,
    //coeffs2[2]=ymH, coeffs3[0]=DeltaCoeff (alphasmu)

    return FermiFunc( y, coeffs1 ) * ( coeffs2[0] / ( coeffs2[1] + y ) ) *
           Delta( ( coeffs2[0] * coeffs2[2] ) / ( coeffs2[1] + y ), coeffs3[0] );
}

double EvtBtoXsgammaKagan::s77FermiFunc( double y,
                                         const std::vector<double>& coeffs1,
                                         const std::vector<double>& coeffs2 )
{
    //coeffs1=fermi function coeffs, coeffs2[0]=mB, coeffs2[1]=mb,
    //coeffs2[2]=ymH
    return FermiFunc( y, coeffs1 ) * ( coeffs2[0] / ( coeffs2[1] + y ) ) *
           s77( ( coeffs2[0] * coeffs2[2] ) / ( coeffs2[1] + y ) );
}

double EvtBtoXsgammaKagan::s88FermiFunc( double y,
                                         const std::vector<double>& coeffs1,
                                         const std::vector<double>& coeffs2,
                                         const std::vector<double>& coeffs3 )
{
    //coeffs1=fermi function coeffs, coeffs2[0]=mB, coeffs2[1]=mb,
    //coeffs2[2]=ymH, coeffs3=s88 coeffs
    return FermiFunc( y, coeffs1 ) * ( coeffs2[0] / ( coeffs2[1] + y ) ) *
           s88( ( coeffs2[0] * coeffs2[2] ) / ( coeffs2[1] + y ), coeffs3[0],
                coeffs3[1] );
}

double EvtBtoXsgammaKagan::s78FermiFunc( double y,
                                         const std::vector<double>& coeffs1,
                                         const std::vector<double>& coeffs2 )
{
    //coeffs1=fermi function coeffs, coeffs2[0]=mB, coeffs2[1]=mb,
    //coeffs2[2]=ymH
    return FermiFunc( y, coeffs1 ) * ( coeffs2[0] / ( coeffs2[1] + y ) ) *
           s78( ( coeffs2[0] * coeffs2[2] ) / ( coeffs2[1] + y ) );
}

double EvtBtoXsgammaKagan::sFermiFunc( double y,
                                       const std::vector<double>& coeffs1,
                                       const std::vector<double>& coeffs2,
                                       const std::vector<double>& coeffs3,
                                       const std::vector<double>& coeffs4 )
{
    //coeffs1=fermi function coeffs, coeffs2[0]=mB, coeffs2[1]=mb,
    //coeffs2[2]=ymH, coeffs3[0]=nIntervals in s22 or s27 array, coeffs3[1]=yMin,
    //coeffs3[2]=yMax, coeffs4=s22 or s27 array
    return FermiFunc( y, coeffs1 ) * ( coeffs2[0] / ( coeffs2[1] + y ) ) *
           GetArrayVal( coeffs2[0] * coeffs2[2] / ( coeffs2[1] + y ),
                        coeffs3[0], coeffs3[1], coeffs3[2], coeffs4 );
}

double EvtBtoXsgammaKagan::Fz( double z )
{
    return ( 1. - 8. * z + 8. * pow( z, 3. ) - pow( z, 4. ) -
             12. * pow( z, 2. ) * log( z ) );
}

double EvtBtoXsgammaKagan::GetArrayVal( double xp, double nInterval, double xMin,
                                        double xMax, std::vector<double> array )
{
    double dx = ( xMax - xMin ) / nInterval;
    int bin1 = int( ( ( xp - xMin ) / ( xMax - xMin ) ) * nInterval );

    double x1 = double( bin1 ) * dx + xMin;

    if ( xp == x1 )
        return array[bin1];

    int bin2( 0 );
    if ( xp > x1 ) {
        bin2 = bin1 + 1;
    } else if ( xp < x1 ) {
        bin2 = bin1 - 1;
    }

    if ( bin1 <= 0 ) {
        bin1 = 0;
        bin2 = 1;
    }

    //If xp is in the last bin, always interpolate between the last two bins
    if ( bin1 == (int)nInterval ) {
        bin2 = (int)nInterval;
        bin1 = (int)nInterval - 1;
        x1 = double( bin1 ) * dx + xMin;
    }

    double x2 = double( bin2 ) * dx + xMin;
    double y1 = array[bin1];

    double y2 = array[bin2];
    double m = ( y2 - y1 ) / ( x2 - x1 );
    double c = y1 - m * x1;
    double result = m * xp + c;

    return result;
}

double EvtBtoXsgammaKagan::FermiFunc( double y, const std::vector<double>& coeffs )
{
    //Fermi shape functions :1=exponential, 2=gaussian, 3=roman
    if ( int( coeffs[0] ) == 1 )
        return EvtBtoXsgammaFermiUtil::FermiExpFunc( y, coeffs );
    if ( int( coeffs[0] ) == 2 )
        return EvtBtoXsgammaFermiUtil::FermiGaussFunc( y, coeffs );
    if ( int( coeffs[0] ) == 3 )
        return EvtBtoXsgammaFermiUtil::FermiRomanFunc( y, coeffs );
    return 1.;
}

double EvtBtoXsgammaKagan::diLogFunc( double y )
{
    return -log( fabs( 1. - y ) ) / y;
}

double EvtBtoXsgammaKagan::diLogMathematica( double y )
{
    double li2( 0 );
    for ( int i = 1; i < 1000;
          i++ ) {    //the value 1000 should actually be Infinite...
        li2 += pow( y, i ) / ( i * i );
    }
    return li2;
}
