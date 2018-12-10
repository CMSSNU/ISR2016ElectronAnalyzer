#!/bin/sh

########################
### SAMPLE LIST ########## 
#######################


source ${LQANALYZER_DIR}/LQRun/txt/list_user_mc.sh


declare -a input_samples=('WZ')

declare -a all_mc=('') 

declare -a mc_noqcd=('') 

declare -a ch_wa=('') 


declare -a ch_wz=('') 

declare -a hn_ll_ee=('') 

declare -a hn_ll_mm=('') 

declare -a hn_ll_em=('') 

declare -a hn_mmm=('') 

declare -a hn_mmm_official=('') 


declare -a hn_moriond_ll=('') 

declare -a hn_ll_tchann_ee=('') 

declare -a hn_ll_tchann_mm=('') 


declare -a diboson_pythia=('WZ' 'ZZ' 'WW')
declare -a diboson_powheg=('WZTo3LNu_powheg' 'ZZTo4L_powheg' 'WWTo2L2Nu' 'WWToLNuQQ')
declare -a diboson_powheg=('WZTo3LNu_amcatnlo')

declare -a dy_mcatnlo=('DYJets_10to50' 'DYJets') 

declare -a dilepton_list=('DYJets_10to50' 'DYJets' 'WJets' 'WpWpEWK' 'WpWpQCD' 'TT_powheg'  'SingleTop_s' 'SingleTbar_t' 'SingleTop_t'  'SingleTbar_tW' 'SingleTop_tW' 'WWW' 'ttW' 'ttZ' 'ttH_nonbb' 'ttH_bb' 'ZZZ' 'WZZ' 'WWZ' 'VBF_HToMuMu' 'WGtoLNuG'  'ZGto2LG' 'WZTo3LNu_powheg' 'ZZTo4L_powheg' 'WWTo2L2Nu' 'WWToLNuQQ' 'TTTT' 'TG' 'TTG' 'ggHtoWW' 'ggHtoZZ' 'vbfHtoWW' 'vbfHtoZZ' 'tZq' 'ww_ds')
declare -a trilepton_list=('DYJets_10to50' 'DYJets' 'WJets' 'WZ' 'ZZ' 'WW'  'TT_powheg' 'SingleTop_s' 'SingleTbar_t' 'SingleTop_t' 'SingleTbar_tW_noHadron'  'SingleTbar_tW' 'SingleTop_tW_noHadron' 'SingleTop_tW' 'WZZ' 'ZZZ' 'WWZ' 'WWW'  'WZTo3LNu_powheg' 'ZZTo4L_powheg'  'WWTo2L2Nu' 'WWToLNuQQ')

declare -a qcd=('') 

declare -a qcd_mu=('') 

declare -a qcd_eg=('') 


declare -a hn_mm=( 'WJets' 'WpWpEWK' 'WpWpQCD' 'TT_powheg'  'SingleTop_s' 'SingleTbar_t' 'SingleTop_t'  'SingleTbar_tW' 'SingleTop_tW' 'WWW' 'ttW' 'ttZ' 'ttH_nonbb' 'ttH_bb' 'ttbb' 'ZZZ' 'WZZ'  'VBF_HToMuMu' 'WGtoLNuG'  'ZGto2LG' 'WZTo3LNu_powheg' 'ZZTo4L_powheg'  'WWTo2L2Nu' 'WWToLNuQQ')

declare -a hn_ee=('DYJets_10to50' 'DYJets''WJets' 'WpWpEWK' 'WpWpQCD' 'TT_powheg'  'SingleTop_s' 'SingleTbar_t' 'SingleTop_t'  'SingleTbar_tW' 'SingleTop_tW' 'WWW' 'ttW' 'ttZ' 'ttH_nonbb' 'ttH_bb' 'ttbb' 'ZZZ' 'WZZ'  'VBF_HToMuMu' 'WGtoLNuG'  'ZGto2LG' 'WZTo3LNu_powheg' 'ZZTo4L_powheg'  'WWTo2L2Nu' 'WWToLNuQQ')

declare -a hn_fakeee=('DYJets_10to50' 'DYJets' 'WJets'  'TT_powheg')


declare -a singletop=('') 

