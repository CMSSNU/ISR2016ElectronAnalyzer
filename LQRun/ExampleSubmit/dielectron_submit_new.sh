sktree -a DYee -S DoubleEG  -list isr_bkg_list -s SKTree_DiLepSkim -n 50 
sktree -a DYee -list isr_dy_list -s FLATCAT -n 50 

# above line is same as the follwing 6 lines together
#sktree -a ExampleAnalyzerDiElectron -list dy_mcatnlo -s SKTree_DiLepSkim -n 15 
#sktree -a ExampleAnalyzerDiElectron -list diboson_pythia -s SKTree_DiLepSkim -n 15 
#sktree -a ExampleAnalyzerDiElectron -list singletop -s SKTree_DiLepSkim -n 15 
#sktree -a ExampleAnalyzerDiElectron -i TT_MG5  -s SKTree_DiLepSkim -n 15 
#sktree -a ExampleAnalyzerDiElectron -i WJets_MCatNLO -s SKTree_DiLepSkim -n 15
#sktree -a ExampleAnalyzerDiElectron  -S DoubleEG  -s SKTree_DiLepSkim -n 15
