#!/bin/bash
LCG_RELEASE=LCG_106 # includes ROOT 6.32, like CMSSW_14_1_0_pre4
# LCG_RELEASE=dev3/latest # includes nightly build of ROOT master, useful for development
LCG_PATH=/cvmfs/sft.cern.ch/lcg/views/$LCG_RELEASE/x86_64-el9-gcc13-opt

source $LCG_PATH/setup.sh
source $LCG_PATH/bin/thisroot.sh