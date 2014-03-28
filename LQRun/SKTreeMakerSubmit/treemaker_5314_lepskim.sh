#!/bin/sh
### sets all configurable variables to defaul values
source functions.sh
###
###########################################################
## CONFIGURE JOB ####
###########################################################
#
#  THIS SECTION IS FOR USER:
#
############################################################

cycle="SKTreeMaker"

stream="egamma"
njobs=30
data_lumi="AtoD"
loglevel="INFO"
logstep=1000
declare -a input_samples=("A" "B" "C" "D" "DY10to50" "DY50plus" "ttbar" "WZ" "ZZ" "WW" "Wjets" "QCD_mumu" "QCD_30-40_EM2" "QCD_40_EM2" "Wgamma" )

### submit this configured job (uses bin/submit.sh)
source submit.sh