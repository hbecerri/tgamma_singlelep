#!/bin/bash
source environment.sh
python3 -m pepper.runproc example_processor.py example_config.json --out output --dataset TGJets_TuneCP5_13TeV-amcatnlo-madspin-pythia8
