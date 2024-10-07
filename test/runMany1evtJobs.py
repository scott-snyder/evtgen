#!/usr/bin/env python
#SBATCH --partition=epp
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3997
#SBATCH --time=12:00:00

import os
import json
import subprocess
import tempfile

nevents = 10000
testname = 'HELAMP=TSS__Bd_Kst2gamma_Kpi'

j = None
with open( 'jsonFiles/config/default.json' ) as jsonFile :
    j = json.load( jsonFile )

with tempfile.TemporaryDirectory() as tmpdirname :
    for i in range(nevents) :

        if i % 100 == 0 :
            print(f'Events remaining: {nevents - i}', flush=True)

        j['rng_seed'] = 1000000*(i+1)
        j['outfile'] = f'{tmpdirname}/{testname}_{i+1}.root'
        j['events'] = 1

        with open( f'{tmpdirname}/tmp.json', 'w' ) as jsonFile :
            json.dump( j, jsonFile )

        result = subprocess.run(['./testDecayModel', 'jsonFiles/{testname}.json', f'{tmpdirname}/tmp.json'], text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

        #logFileName = j['outfile'].replace('.root', '.out')
        #with open( logFileName, 'w' ) as logFile :
            #logFile.write( result.stdout )

    cmd_args = ['hadd', '-ff', '{testname}.root']

    tmpdirlisting = os.listdir(f'{tmpdirname}')
    inputfiles = [ f'{tmpdirname}/{file}' for file in tmpdirlisting if file.endswith('.root') ]
    cmd_args.extend( inputfiles )

    result = subprocess.run(cmd_args, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    with open( 'hadd.log', 'w' ) as logFile :
        logFile.write( result.stdout )
