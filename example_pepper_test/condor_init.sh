#!/bin/bash
source ~/.bashrc
#conda activate pepenv

# Load CMS grid environment
source /cvmfs/grid.desy.de/etc/profile.d/grid-ui-env.sh > /dev/null
# Set custom VOMS proxy path. This needs to be accessible from Condor
# Please do not forget to run voms-proxy-init --voms cms --out $X509_USER_PROXY
export X509_USER_PROXY=~/.globus/x509up

# Load LCG
#source /cvmfs/sft.cern.ch/lcg/views/setupViews.sh LCG_102 x86_64-centos7-gcc11-opt

# Make sure python libs installed in the user directory are prefered over system-wide ones
#export PYTHONPATH=`python3 -c 'import site; print(site.getusersitepackages())'`:$PYTHONPATH

# Parsl installs some of its commands into ~/.local/bin if installed as user
export PATH=~/.local/bin:$PATH

export PYTHONPATH=/afs/desy.de/user/y/yian/cms/topgamma/pepper/test/:$PYTHONPATH
