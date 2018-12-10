#sktree -a ISRee -S DoubleEG  -list isr_bkg_list -s SKTree_DiLepSkim -n 50 
#sktree -a ISRee -S DoubleEG -list isr_bkg_list -s SKTree_DiLepSkim -n 50 
sktree -a ISRee -list isr_dy_list -s FLATCAT -n 1 

#sktree -a ISRmumu -list isr_bkg_list -s SKTree_DiLepSkim -n 50 
#sktree -a ISRmumu -S DoubleMuon  -list isr_bkg_list -s SKTree_DiLepSkim -n 50 
#sktree -a ISRmumu -list isr_dy_list -s FLATCAT -n 50 

#sktree -a ISRee_unfolding -S DoubleEG  -list isr_bkg_list -s SKTree_DiLepSkim -n 50 
#sktree -a ISRee_unfolding -list isr_dy_list -s FLATCAT -n 50 

#sktree -a ISRee_FSRunfolding -list isr_dy_list -s FLATCAT -n 50 

# above line is same as the follwing 6 lines together
#sktree -a ExampleAnalyzerDiElectron -list dy_mcatnlo -s SKTree_DiLepSkim -n 15 
#sktree -a ExampleAnalyzerDiElectron -list diboson_pythia -s SKTree_DiLepSkim -n 15 
#sktree -a ExampleAnalyzerDiElectron -list singletop -s SKTree_DiLepSkim -n 15 
#sktree -a ExampleAnalyzerDiElectron -i TT_MG5  -s SKTree_DiLepSkim -n 15 
#sktree -a ExampleAnalyzerDiElectron -i WJets_MCatNLO -s SKTree_DiLepSkim -n 15
#sktree -a ExampleAnalyzerDiElectron  -S DoubleEG  -s SKTree_DiLepSkim -n 15
