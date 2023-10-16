#!/usr/bin/env python

import json
import itertools

class JsonFileChecker :

    particleNames = {
        'B0' : [ 'Bd', 'B0' ],
        'anti-B0' : [ 'Bd', 'Bdbar', 'B0', 'B0bar' ],
        'B_s0' : [ 'Bs', 'Bs0', 'B0s' ],
        'anti-B_s0' : [ 'Bs', 'Bsbar', 'Bs0', 'Bs0bar', 'B0s', 'B0sbar' ],
        'B+' : [ 'Bu', 'Bp', 'B+' ],
        'B-' : [ 'Bu', 'Bm', 'B-' ],
        'B_c+' : [ 'Bc', 'Bcp', 'Bc+' ],
        'B_c-' : [ 'Bc', 'Bcm', 'Bc-' ],
        'B*0' : [ 'Bstd', 'Bst0' ],
        'anti-B*0' : [ 'Bstd', 'Bst0' ],
        'B_s*0' : [ 'Bsst', 'Bs0st', 'Bs0star' ],
        'anti-B_s*0' : [ 'Bsst', 'Bs0st', 'Bs0star' ],
        'D+' : [ 'D+', 'Dp', 'D' ],
        'D-' : [ 'D-', 'Dm', 'D' ],
        'D0' : [ 'D0', 'Dz' ],
        'anti-D0' : [ 'D0', 'D0bar', 'anti-D0', 'Dz', 'Dzbar', 'anti-Dz' ],
        'D_s+' : [ 'Ds+', 'Dsp', 'Ds' ],
        'D_s-' : [ 'Ds-', 'Dsm', 'Ds' ],
        'D*+' : [ 'Dst+', 'Dstp', 'Dst', 'Dstar+', 'Dstarp', 'Dstar' ],
        'D*-' : [ 'Dst-', 'Dstm', 'Dst', 'Dstar-', 'Dstarm', 'Dstar' ],
        'D*0' : [ 'Dst0', 'Dstz', 'Dstar0', 'Dstarz' ],
        'anti-D*0' : [ 'Dst0', 'Dstz', 'Dstar0', 'Dstarz' ],
        'D_2*0' : [ 'D2st0', 'D2stz', 'D2star0', 'D2starz' ],
        'anti-D_2*0' : [ 'D2st0', 'D2stz', 'D2star0', 'D2starz' ],
        'D*(2S)+' : [ 'Dst2S+', 'Dst2Sp', 'Dst2S', 'Dstar2S+', 'Dstar2Sp', 'Dstar2S' ],
        'D*(2S)-' : [ 'Dst2S-', 'Dst2Sm', 'Dst2S', 'Dstar2S-', 'Dstar2Sm', 'Dstar2S' ],
        'K+' : [ 'K+', 'Kp', 'K' ],
        'K-' : [ 'K-', 'Km', 'K' ],
        'K0' : [ 'K0', 'Kz' ],
        'anti-K0' : [ 'K0', 'Kz', 'anti-K0' ],
        'K_S0' : [ 'KS', 'KS0', 'K0S', 'Ks', 'Ks0', 'K0s' ],
        'K_L0' : [ 'KL', 'KL0', 'K0L', 'Kl', 'Kl0', 'K0l' ],
        'K*+' : [ 'Kst+', 'Kstp', 'Kst', 'Kstar' ],
        'K*-' : [ 'Kst-', 'Kstm', 'Kst', 'Kstar' ],
        'K*0' : [ 'Kst0', 'Kstz', 'Kst', 'Kstar0' ],
        'anti-K*0' : [ 'Kst0', 'Kstz', 'Kst', 'Kstar0' ],
        'K_1+' : [ 'K1+', 'K1p', 'K1' ],
        'K_1-' : [ 'K1-', 'K1m', 'K1' ],
        'K\'_1+' : [ 'Kprime1+', 'Kprime1p', 'Kprime1' ],
        'K\'_1-' : [ 'Kprime1-', 'Kprime1m', 'Kprime1' ],
        'pi+' : [ 'pi+', 'pip', 'pi' ],
        'pi-' : [ 'pi-', 'pim', 'pi' ],
        'pi0' : [ 'pi0', 'piz' ],
        'rho+' : [ 'rho+', 'rhop', 'rho' ],
        'rho-' : [ 'rho-', 'rhom', 'rho' ],
        'rho0' : [ 'rho0', 'rhoz' ],
        'a_1+' : [ 'a1+', 'a1p', 'a1' ],
        'a_1-' : [ 'a1-', 'a1m', 'a1' ],
        'a_2+' : [ 'a2+', 'a2p', 'a2' ],
        'a_2-' : [ 'a2-', 'a2m', 'a2' ],
        'f_2' : [ 'f2' ],
        'p+' : [ 'p+', 'p' ],
        'anti-p-' : [ 'p-', 'p', 'pbar' ],
        'Lambda_b0' : [ 'Lambdab', 'Lambdab0', 'Lb', 'Lbz', 'Lb0' ],
        'Lambda_c+' : [ 'Lambdac', 'Lambdac+', 'Lambdacp', 'Lc', 'Lc+', 'Lcp' ],
        'Lambda_c(2593)+' : [ 'Lambdac2593', 'Lambdac2593+', 'Lambdac2593p', 'Lc2593', 'Lc2593+', 'Lc2593p' ],
        'Lambda_c(2625)+' : [ 'Lambdac2625', 'Lambdac2625+', 'Lambdac2625p', 'Lc2625', 'Lc2625+', 'Lc2625p' ],
        'Lambda0' : [ 'Lambda', 'Lambda0', 'Lz', 'L0' ],
        'Delta+' : [ 'Delta+', 'Deltap' ],
        'J/psi' : [ 'Jpsi' ],
        'psi(2S)' : [ 'psi2S' ],
        'psi(3770)' : [ 'psi3770' ],
        'psi(4160)' : [ 'psi4160' ],
        'psi(4415)' : [ 'psi4415' ],
        'X_2(3872)' : [ 'X3872' ],
        'chi_c0' : [ 'chic0' ],
        'chi_c1' : [ 'chic1' ],
        'chi_c2' : [ 'chic2' ],
        'Upsilon' : [ 'Upsilon1S', 'Y1S', 'Upsilon' ],
        'Upsilon(2S)' : [ 'Upsilon2S', 'Y2S' ],
        'Upsilon(3S)' : [ 'Upsilon3S', 'Y3S' ],
        'Upsilon(4S)' : [ 'Upsilon4S', 'Y4S' ],
        'Upsilon(5S)' : [ 'Upsilon5S', 'Y5S' ],
        'omega' : [ 'omega' ],
        'phi' : [ 'phi' ],
        'eta' : [ 'eta' ],
        'eta\'' : [ 'etap', 'etapr', 'etaprime' ],
        'e+' : [ 'e', 'ep', 'e+' ],
        'e-' : [ 'e', 'em', 'e-' ],
        'mu+' : [ 'mu', 'mup', 'mu+' ],
        'mu-' : [ 'mu', 'mum', 'mu-' ],
        'tau+' : [ 'tau', 'taup', 'tau+' ],
        'tau-' : [ 'tau', 'taum', 'tau-' ],
        'nu_e' : [ 'nue', 'nu-e', 'nu' ],
        'anti-nu_e' : [ 'anti-nue', 'nue', 'nu-e', 'nu' ],
        'nu_mu' : [ 'numu', 'nu-mu', 'nu' ],
        'anti-nu_mu' : [ 'anti-numu', 'numu', 'nu-mu', 'nu' ],
        'nu_tau' : [ 'nutau', 'nu-tau', 'nu' ],
        'anti-nu_tau' : [ 'anti-nutau', 'nutau', 'nu-tau', 'nu' ],
        'gamma' : [ 'gamma' ],
        'Xsd' : [ 'Xsd' ],
        'Xsu' : [ 'Xsu' ],
        'Xu-' : [ 'Xu' ],
        }

    allowedSuffixes = [
            'cuts',
            'CLEO',
            'BaBar',
            'Longitudinal',
            'Transverse',
            'tree',
            'nlo',
            'pWave',
            'sWave',
            'yesPhotos'
            ]

    def __init__( self, jsonFileName ) :
        self.j = None
        with open( 'jsonFiles/'+jsonFileName ) as jsonFile :
            self.j = json.load( jsonFile )

        self.jsonFileName = jsonFileName
        self.modelString = jsonFileName[ : jsonFileName.index('__') ]
        self.particleString = jsonFileName[ jsonFileName.index('__') + 2 :  jsonFileName.index('.json') ]
        self.jsonFileNameBase = jsonFileName[ : jsonFileName.index('.json') ]

    def checkModels( self ) :

        modelsFromString = self.modelString.split('=')

        jsonModelsList = [ i for i in itertools.filterfalse( lambda x : x == "", self.j['models'] ) ]

        if len(modelsFromString) != len(jsonModelsList) :
            print(f'ERROR : different numbers of models: {len(modelsFromString)} in filename, {len(jsonModelsList)} in JSON')
            return False

        allOK = True

        for model1, model2 in zip( modelsFromString, jsonModelsList ) :
            if  model1 != model2 :
                print(f'ERROR : different model name : {model1} in filename, {model2} in JSON')
                allOK = False

        return allOK


    def checkParticles( self ) :

        allOK = True

        models = self.j['models']

        generations = self.particleString.split('_')
        if generations[-1] in self.allowedSuffixes :
            generations = generations[:-1]

        parent = self.j['parent']

        daughters = self.j['daughters']
        if 'VSS_BMIX' == models[0] and len(daughters) == 4 and daughters[0] == daughters[2] and daughters[1] == daughters[3] :
            daughters = daughters[:2]

        grandDaughters = None
        if 'grand_daughters' in self.j :
            grandDaughters = self.j['grand_daughters']
        if ( grandDaughters and len(generations) != 3 ) or ( not grandDaughters and len(generations) != 2 ) :
            print(f'ERROR : incorrect number of generations in particle string: {generations}')
            return False

        if generations[0] not in self.particleNames[ parent ] :
            print(f'ERROR : parent name mis-match : {generations[0]} in filename, {parent} in JSON')
            allOK = False

        tauolaIndex = -1
        tauolaDecay = -1
        if 'TAUOLA' in models :
            tauolaIndex = models.index('TAUOLA')
            tauolaDecay = self.j['parameters'][tauolaIndex][0]

        if tauolaIndex == 0 :
            print(f'WARNING : cannot check accuracy of daughter string for TAUOLA decay {tauolaDecay}')
        else :
            daugList = [ self.particleNames[d] for d in daughters ]

            combinations = itertools.product( *daugList )

            string_combinations = []
            for comb in combinations :
                string = ''
                for name in comb :
                    string += name
                string_combinations.append( string )

            if generations[1] not in  string_combinations :
                print(f'ERROR : daughter name mis-match : {generations[1]} in filename, {daughters} in JSON')
                print(string_combinations)
                allOK = False

        if grandDaughters :
            gdStrings = generations[2].split(',')

            gdLists = [ i for i in itertools.filterfalse( lambda x : len(x) == 0, grandDaughters ) ]

            if len(gdStrings) != len(gdLists) :
                print(f'ERROR : mis-match in number of grand-daughter sets: {gdStrings} in filename, {gdLists} in JSON')
                return False

            for index, gdList in enumerate(gdLists) :
                if 'VSS_BMIX' == models[index+1] and len(gdList) == 4 and gdList[0] == gdList[2] and gdList[1] == gdList[3] :
                    gdList = gdList[:2]

            if tauolaIndex > 0 :
                gdStrings.pop(tauolaIndex-1)
                gdLists.pop(tauolaIndex-1)
                print(f'WARNING : cannot check accuracy of grand-daughter string for TAUOLA decay {tauolaDecay}')

            for gdString, gdList in zip( gdStrings, gdLists ) :

                grandDaugList = [ self.particleNames[gd] for gd in gdList ]

                combinations = itertools.product( *grandDaugList )

                string_combinations = []
                for comb in combinations :
                    string = ''
                    for name in comb :
                        string += name
                    string_combinations.append( string )

                if gdString not in  string_combinations :
                    print(f'ERROR : grand-daughter name mis-match : {gdString} in filename, {gdList} in JSON')
                    print(string_combinations)
                    allOK = False

        return allOK


    def checkFileNames( self ) :

        rootFileName = self.j['outfile']
        if rootFileName != self.jsonFileNameBase+'.root' :
            print(f'ERROR : ROOT file name mis-match : {self.jsonFileNameBase}.root from JSON filename, {rootFileName} in JSON field')
            return False

        refFileName = self.j['reference']
        if refFileName != 'Ref'+self.jsonFileNameBase+'.root' :
            print(f'ERROR : reference file name mis-match : Ref{self.jsonFileNameBase}.root from JSON filename, {refFileName} in JSON field')
            return False

        return True


    def runChecks( self ) :

        modelsOK = checker.checkModels()

        particlesOK = checker.checkParticles()

        namesOK = checker.checkFileNames()

        return modelsOK and particlesOK and namesOK



if __name__ == '__main__' :

    import os, sys

    jsonFiles = os.listdir('jsonFiles')

    allOK = True

    for jsonFileName in jsonFiles :

        print(f'Checking file: {jsonFileName}')

        checker = JsonFileChecker( jsonFileName )

        if not checker.runChecks() :
            allOK = False

    if not allOK :
        sys.exit(1)

