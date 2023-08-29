#!/bin/bash

#For condor
python3 -m pepper.runproc tgamma_processor.py config/test_config.json --condor 30 --condorinit environment.sh --dataset TGJets_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8 -R --out condor_output --eventdir condor_output
