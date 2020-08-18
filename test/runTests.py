#!/usr/bin/env python

########################################################################
# Copyright 1998-2022 CERN for the benefit of the EvtGen authors       #
#                                                                      #
# This file is part of EvtGen.                                         #
#                                                                      #
# EvtGen is free software: you can redistribute it and/or modify       #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation, either version 3 of the License, or    #
# (at your option) any later version.                                  #
#                                                                      #
# EvtGen is distributed in the hope that it will be useful,            #
# but WITHOUT ANY WARRANTY; without even the implied warranty of       #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        #
# GNU General Public License for more details.                         #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with EvtGen.  If not, see <https://www.gnu.org/licenses/>.     #
########################################################################

import os
import sys
import json
import subprocess

# Usage: runTests.py <model json file> <dependency json file> <changed source file 1> [[changed source file 2] ...]
if len(sys.argv) < 4 :
    print('ERROR: must provide model and dependency json files and one or more modified files')
    sys.exit(1)

modelJson = sys.argv[1]
depJson = sys.argv[2]
changedFiles = sys.argv[3:]

# load the json files that contain:
# - the associations between each model and the header/source pair that define the model class
# - the associations between each source file and the other header/source files on which it depends
with open(modelJson) as f1 :
    models = json.load(f1)

with open(depJson) as f2 :
    srcdeps = json.load(f2)

if not models or not srcdeps or len(models) == 0 or len(srcdeps) == 0 :
    print('ERROR: problem loading model dependency information')
    sys.exit(1)

# turn the lists into sets, automatically removes duplicates
for srcfile in srcdeps :
    depset = set(srcdeps[srcfile])
    srcdeps[srcfile] = depset

def getDeps( srcfile, srcdeps, founddeps ) :
    """
    Find all the dependencies for a model
    Arguments:
    srcfile - the model's source file name (or, in recursion, a source file on which it depends)
    srcdeps - the dictionary of all source files' dependencies
    founddeps - a set of dependencies already found for the model
    """
    founddeps.add(srcfile) # make sure our source file is already marked as found
    deps = set(srcdeps[srcfile]) # NB need to make a copy here

    # remove already found dependencies from deps,
    # then add new ones to founddeps
    deps.difference_update(founddeps)
    founddeps.update(deps)

    # make list of source files in the dependencies
    cpps = []
    for f in deps :
        if f.endswith('.cpp') :
            cpps.append(f)

    for f in cpps :
        # add the dependencies of this source to our set
        deps.update( getDeps( f, srcdeps, founddeps ) )

    return deps

# create the inverse look-up dictionary
fileDeps = {}
for model in models :
    hdrfile = models[model]['header']
    srcfile = models[model]['source']

    founddeps = set([])
    modeldeps = getDeps(srcfile,srcdeps,founddeps)
    modeldeps.add(hdrfile)
    modeldeps.add(srcfile)

    for modeldep in modeldeps :
        if modeldep in fileDeps :
            fileDeps[modeldep].append( model )
        else :
            fileDeps[modeldep] = [ model ]

testAll = False
modelsToTest = []
skippedFiles = []
jsonFilesChanged = []

for changedFile in changedFiles :
    # if one of the test json files has changed we need to re-run that test
    if changedFile.startswith('test/jsonFiles/') and changedFile.endswith('.json') :
        jsonFilesChanged.append( changedFile )

    # otherwise only consider source code files (TODO - any others that we should consider? evt.pdl, DECAY.DEC and similar?)
    if not ( changedFile.endswith('.hh') or changedFile.endswith('.cpp') or changedFile.endswith('.cc') ) :
        skippedFiles.append( changedFile )
        continue

    # the EvtGen, EvtModelReg, or EvtNoRadCorr objects; or the test code itself should trigger a rerun of all tests
    # TODO - any others that should be in here?
    if changedFile == 'EvtGen/EvtGen.hh' or changedFile == 'src/EvtGen.cpp' or changedFile == 'EvtGenModels/EvtModelReg.hh' or changedFile == 'src/EvtGenModels/EvtModelReg.cpp' or changedFile == 'EvtGenModels/EvtNoRadCorr.hh' or changedFile.startswith('test/') :
        testAll = True
        modelsToTest = []
        break

    # TODO temporarily skip EvtGenExternal stuff
    if 'EvtGenExternal/' in changedFile or changedFile == 'EvtGenModels/EvtAbsExternalGen.hh' :
        skippedFiles.append( changedFile )
        continue

    # otherwise check in the model dependency dict
    if changedFile not in fileDeps :
        print(f'WARNING: no dependency information for modified file: {changedFile}')
        continue

    models = fileDeps[changedFile]
    if len(models) == 0 :
        print(f'WARNING: no models listed for modified file: {changedFile}')
        continue

    modelsToTest.extend( models )

if testAll :
    print('Need to test all models')
    modelsToTest = list(models.keys())

if len(modelsToTest) == 0 and len(jsonFilesChanged) == 0 :
    print('No models to test, exiting...')
    if len(skippedFiles) == len(changedFiles) :
        sys.exit(0)
    else :
        sys.exit(1)

if len(jsonFilesChanged) != 0 :
    print(f'The following test configs have changed and so need to be run: {jsonFilesChanged}')

modelsToTest = set(modelsToTest)
if len(modelsToTest) != 0 :
    print(f'Need to test the following models: {modelsToTest}')

testFilesToRun = []
allTestFiles = os.listdir('jsonFiles')
for model in modelsToTest :
    testFilesForThisModel = []
    for testFile in allTestFiles :
        if not testFile.endswith('.json') :
            continue
        testModelNames = testFile[:testFile.index('__')].split('=')
        if model in testModelNames :
            testFilesForThisModel.append( testFile )
    if len(testFilesForThisModel) == 0 :
        print(f'WARNING: no tests available for model: {model}')
    else :
        testFilesToRun.extend( testFilesForThisModel )

for jsonFile in jsonFilesChanged :
    testFilesToRun.append( jsonFile.split('/')[-1] )

if len(testFilesToRun) == 0 :
    print(f'No tests to run, exiting...')
    sys.exit(1)

testFilesToRun = set(testFilesToRun)
print(f'Will run the following tests: {testFilesToRun}')

joinOutput = True
allOK = True
for testFile in testFilesToRun :
    if joinOutput :
        result = subprocess.run(['./testDecayModel', 'jsonFiles/'+testFile], text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    else :
        result = subprocess.run(['./testDecayModel', 'jsonFiles/'+testFile], text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    outFileName = testFile.replace('.json','.out')
    with open(outFileName,'w') as outFile :
        outFile.write(result.stdout)

    if not joinOutput :
        errFileName = testFile.replace('.json','.err')
        with open(errFileName,'w') as errFile :
            errFile.write(result.stderr)

    if result.returncode == 0 :
        print(f'Successfully ran test: {testFile}')
    else :
        print(f'ERROR: problem running test: {testFile}')
        if joinOutput :
            print(f'     : check the log file: {outFileName}')
        else :
            print(f'     : check the log files: {outFileName} and {errFileName}')
        allOK = False

if not allOK :
    sys.exit(1)
