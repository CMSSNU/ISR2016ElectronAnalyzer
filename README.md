# Steps to install and use this analyser

## Install
> git clone git@github.com:CMSSNU/ISR2016ElectronAnalyzer.git

## Basic setup 

> cd ISR2016ElectronAnalyzer/

> source setup.sh 

After the above line, follow the instructions to register an email address.

## Run analyser

> cd LQRun/ExampleSubmit

> source dielectron_submit_new.sh

Below instruction is from John.

//###########################
//   Author John Almond
//###########################
Steps to run analysis:
1) Download LQanalyzer:
git clone git@github.com:jalmond/LQanalyzer.git
cd LQanalyzer/
git checkout -b BRANCH_NAME origin/BRANCH_NAME
 
2) in LQanalyzer main directory type:
source setup.sh

## NOTE


NEW for CATAnalyzer::

_______________________
Triggers:
_______________________
The following triggers are included in new SKtrees:

example: 
std::vector<TString> triggerslist;
triggerslist.push_back("HLT_IsoMu24_eta2p1_v");
if(!PassTrigger(triggerslist, prescale)) return;


_______________________
TriggerMatching:
_______________________
example:
KMuon lep; || KElectron lep;
lep.TriggerMatched("HLT_IsoMu24_eta2p1") 

can use any of these:
 ///HLT_IsoMu24_eta2p1_v
 ///HLT_Mu17_Mu8_DZ_v
 ///HLT_Mu17_TkMu8_DZ_v
 ///HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v
 ///HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v
 ///HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v
 ///HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v
 ///HLT_Ele12_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v
 ///HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v
 ///HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v
 ///HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v
 ///HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v
 ///HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v
 ///HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v
 ///HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v
 ///HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_
 ///HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v
 ///HLT_Ele27_eta2p1_WPLoose_Gsf_TriCentralPFJet30_v




/// BTAGGING : https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation76X
pfCombinedInclusiveSecondaryVertexV2BJetTags	
CombinedSecondaryVertex v2	CSVv2L	see KJet.h for WP
CombinedSecondaryVertex v2	CSVv2M	
CombinedSecondaryVertex v2	CSVv2T	



/// Pileup ID jets: use 13 TeV ID
https://twiki.cern.ch/twiki/bin/view/CMS/JPTPileupJetID
https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2014
