#!/bin/sh

rundata=true
runmc=true
runsig=true
runfakes=true
runflips=true

if [[ $1  == "NP" ]];
    then
    runfakes=true
    rundata=false
    runmc=false
    runsig=false
    runflips=false
fi

if [[ $1  == "CF" ]];
    then
    runfakes=false
    rundata=false
    runmc=false
    runsig=false
    runflips=true
fi

if [[ $1  == "DATA" ]];
    then
    runfakes=false
    rundata=true
    runmc=false
    runsig=false
    runflips=false
fi

if [[ $1  == "MC" ]];
    then
    runfakes=false
    rundata=false
    runmc=true
    runsig=false
    runflips=false
fi




if [[ $runmc  == "true" ]];

then

    source functions.sh
    
    cycle="HNDiElectron"
    skinput="True"
    
    njobs=30
    data_lumi="AtoD"
    loglevel="INFO"

    loglevel="INFO"
    logstep=1000
    
    declare -a input_samples=("WZtollqq_mg" "WZtoqqln_mg" "WZtollln_mg" "ZZtollnn_mg" "ZZtollqq_mg" "ZZtollll_mg" "SSWmWm" "SSWpWp" "WW_dp" "ttW" "ttZ" "Wjets" "ttbar" "stbar_sch" "stbar_tch" "stbar_tW" "st_sch" "st_tch" "st_tW" "Wgamma" "HtoZZ" "Zgamma" "DY10to50" "DY50plus" "Zbb" "WW_mg" "ggHtoZZ" "HtoTauTau" "HtoWW" "WgammaE" "WgammaTau" "WWW" "TTWW" "TTG" "ZZZ" "WZZ" "WWZ" "WWG" "WZ_py" "ZZ_py" "WW_py" "QCD_20_30_EM" "QCD_20_30_BCtoE" "QCD_30_80_EM" "QCD_30_80_BCtoE" "QCD_80_170_EM" "QCD_80_170_BCtoE" "QCD_170_250_EM" "QCD_170_250_BCtoE" "QCD_250_350_EM" "QCD_250_350_BCtoE" "QCD_30-40_EM2" "QCD_40_EM2")
    
    outputdir=$LQANALYZER_DIR"/data/output/SSElectronSingleEl/"
    ### submit this configured job (uses bin/submit.sh)
    source submit.sh $1
fi


if [[ $rundata  == "true" ]];
then

    source functions.sh

    cycle="HNDiElectron"
    skinput="True"

    njobs=30
    data_lumi="AtoD"
    loglevel="INFO"

    loglevel="INFO"
    logstep=1000

    declare -a input_samples=( "A" "B" "C" "D")
    stream="egamma"
    
    outputdir=$LQANALYZER_DIR"/data/output/SSElectronSingleEl/"
    ### submit this configured job (uses bin/submit.sh)
    source submit.sh $1
    source hadd.sh ${LQANALYZER_DIR}/data/output/SSElectronSingleEl/  HNDiElectron_data_5_3_14.root HNDiElectron_period*
fi


if [[ $runfakes  == "true" ]];
    then
    source functions.sh
    
    cycle="HNDiElectron"
    skinput="True"
    loglevel="INFO"

    njobs=30
    data_lumi="AtoD"
    
    loglevel="INFO"
    logstep=1000
    
    runnp="True"
    declare -a input_samples=("A" "B" "C" "D") 
    
    stream="egamma"
    outputdir=$LQANALYZER_DIR"/data/output/SSElectronSingleEl/"


    ### submit this configured job (uses bin/submit.sh)
    source submit.sh $1
    source hadd.sh ${LQANALYZER_DIR}/data/output/SSElectronSingleEl/  HNDiElectron_SKnonprompt_5_3_14.root HNDiElectron_nonprompt_period*
fi

if [[ $runflips  == "true" ]];
    then
    source functions.sh
    
    cycle="HNDiElectron"
    skinput="True"
    loglevel="INFO"
    
    njobs=30
    data_lumi="AtoD"

    loglevel="INFO"
    logstep=1000

    runcf="True"
    declare -a input_samples=("A" "B" "C" "D" )
    
    
    stream="egamma"
    outputdir=$LQANALYZER_DIR"/data/output/SSElectronSingleEl/"
    
    
    ### submit this configured job (uses bin/submit.sh)
    source submit.sh $1
    source hadd.sh ${LQANALYZER_DIR}/data/output/SSElectronSingleEl/  HNDiElectron_SKchargeflip_5_3_14.root HNDiElectron_chargeflip_*
fi



echo ""
echo "End of example_submit.sh script."