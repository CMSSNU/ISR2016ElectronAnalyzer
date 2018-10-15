// $Id: ISRee.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQISRee Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "ISRee.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (ISRee);


 /**
  *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
  *
  */
ISRee::ISRee() :  AnalyzerCore(), out_electrons(0)  {
  
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("ISRee");
  
  Message("In ISRee constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();
  //  MakeCleverHistograms(sighist_mm,"DiElectron");


}


void ISRee::InitialiseAnalysis() throw( LQError ) {
  
  /// Initialise histograms
  MakeHistograms();  
  //
  // You can out put messages simply with Message function. Message( "comment", output_level)   output_level can be VERBOSE/INFO/DEBUG/WARNING 
  // You can also use m_logger << level << "comment" << int/double  << LQLogger::endmsg;
  //
  
  Message("Making clever hists for Z ->ll test code", INFO);
  
  /// only available in v7-6-X branch and newer
  //// default lumimask is silver ////
  //// In v7-6-2-(current) the default is changed to gold (since METNoHF bug)
  ///When METNoHF isfixed the default will be back to silver
  /// set to gold if you want to use gold json in analysis
  /// To set uncomment the line below:


  return;
}


void ISRee::ExecuteEvents()throw( LQError ){

  // double electron trigger
  TString dielectron_trig="HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v";//

  // Now you should do an OR of 4 triggers 
  vector<TString> trignames;
  trignames.push_back(dielectron_trig);

  bool event_m80to100 = false;
  double weigtbytrig = WeightByTrigger(dielectron_trig,TargetLumi);

  // Mass boundaries  
  double massBinBoundaries[] = {20., 30., 40., 50., 60., 70., 80., 100., 200., 350};
  const int nBoundaries = sizeof(massBinBoundaries)/sizeof(double); 

  vector<TString> massbins;
  massbins.push_back("m20to30");
  massbins.push_back("m30to40");
  massbins.push_back("m40to50");
  massbins.push_back("m50to60");
  massbins.push_back("m60to70");
  massbins.push_back("m70to80");
  massbins.push_back("m80to100");
  massbins.push_back("m100to200");
  massbins.push_back("m200to350");

  //double genmassBinBoundaries[] = {10.,20.,30.,40.,50.,60.,70.,80.,90.,100.};
  //const int gennBoundaries = sizeof(genmassBinBoundaries)/sizeof(double);

  int genmassbin = 90;
  double genmassBinBoundaries[genmassbin+1];

  for(int i = 0; i < genmassbin+1; i++){
     genmassBinBoundaries[i] = (double) (i + 10);
  }
  const int gennBoundaries = sizeof(genmassBinBoundaries)/sizeof(double);

  vector<TString> genmassbins;

  for(int i = 0; i < genmassbin+1; i++){
     TString lowMass, highMass;
     lowMass.Form("%d", (int) genmassBinBoundaries[i]);
     highMass.Form("%d", (int) genmassBinBoundaries[i+1]);
     TString name = "m"+lowMass+"to"+highMass;
     genmassbins.push_back(name);
  }

  //genmassbins.push_back("m10to20");
  //genmassbins.push_back("m20to30");
  //genmassbins.push_back("m30to40");
  //genmassbins.push_back("m40to50");
  //genmassbins.push_back("m50to60");
  //genmassbins.push_back("m60to70");
  //genmassbins.push_back("m70to80");
  //genmassbins.push_back("m80to90");
  //genmassbins.push_back("m90to100");

  /// Apply the gen weight 
  if(!isData) weight*=MCweight;

  // Get generator information for gen/reco ratio correction (only for Drell-Yan MC).
  // For this, we needs full-phase space generator level information.
  // So, You should run this with FLATCAT.
  if(k_sample_name.Contains("DY")){  //only for Drell-Yan MC
    std::vector<snu::KTruth> truthcol =  eventbase->GetTruth();   //get truth particles
    TLorentzVector gendy,genl1,genl2;   //gendy: Z/gamma* fourvector, genl1: electron fourvector, genl2: anti-electron fourvector
    TLorentzVector gendy_postFSR,genl1_postFSR,genl2_postFSR;   
    TLorentzVector gendy_preFSR,genl1_preFSR,genl2_preFSR;   
    int truthsize=truthcol.size();

    //loop for collect dimuon Drell-Yan product
    for(int i=0;i<truthsize;i++){
      snu::KTruth truth=truthcol.at(i);
      if(truth.GenStatus()!=1) continue;  //stable-particle-requirement
      if(truth.PdgId()==11&&truth.StatusFlag(snu::KTruth::fromhardprocess)){
	genl1+=truth;
	gendy+=truth;
        genl1_postFSR+=truth;
	gendy_postFSR+=truth;

        genl1_preFSR+=truth;
      }else if(truth.PdgId()==-11&&truth.StatusFlag(snu::KTruth::fromhardprocess)){
	genl2+=truth;
	gendy+=truth;
	genl2_postFSR+=truth;
	gendy_postFSR+=truth;

        genl2_preFSR+=truth;
      }else if(truth.PdgId()==22){   //collect photons from DY muons by FSR
	int imother=truth.IndexMother();
	if(imother>=truthsize) continue;
	snu::KTruth mother=truthcol.at(truth.IndexMother());
	if(abs(mother.PdgId())==11&&mother.StatusFlag(snu::KTruth::fromhardprocess)){
	  gendy+=truth;
	}

        if(mother.PdgId() == 11 && mother.StatusFlag(snu::KTruth::fromhardprocess)){
          genl1_preFSR+=truth;
        }
        if(mother.PdgId() == -11 && mother.StatusFlag(snu::KTruth::fromhardprocess)){
          genl2_preFSR+=truth;
        }
      }

    }
   
    //Fill generator level hists, if we find both electron  and anti-electron
    if(genl1!=TLorentzVector(0,0,0,0)&&genl2!=TLorentzVector(0,0,0,0)){  
      double genptlep1 = genl1.Pt();
      double genptlep2 = genl2.Pt();
      double genetalep1 = genl1.Eta();
      double genetalep2 = genl2.Eta();
      double gendipt = gendy.Pt();
      double gendimass = gendy.M();  

      if(gendimass > 80. && gendimass < 100.){ 
         FillCutFlow("NoCut", 1.); 
        if((genl1_preFSR.Pt()> 25 && genl2_preFSR.Pt()> 15)||(genl2_preFSR.Pt()> 25 && genl1_preFSR.Pt()> 15) && fabs(genl1_preFSR.Eta()) < 2.4 && fabs(genl2_preFSR.Eta() && gendipt < 100.) < 2.4){
          event_m80to100 = true;
          FillCutFlow("GenKinCut", 1.); 
        }
      }
      
      FillHist("gendipt",gendipt,weight*weigtbytrig,0.,100.,100);
      FillHist("gendipt_check",(genl1_preFSR+genl2_preFSR).Pt(),weight*weigtbytrig,0.,100.,100);
      FillHist("gendimass",gendimass,weight*weigtbytrig,0.,500.,1000);

      TString postfix_preFSR="";
      for(int i = 0; i < gennBoundaries-1; i++){
         if(gendy.M() > genmassBinBoundaries[i] && gendy.M() < genmassBinBoundaries[i+1]){
           postfix_preFSR=genmassbins[i];
         }
      }// loop for mass bins

      FillHist("gendipt_preFSR",gendy.Pt(),weight*weigtbytrig,0.,100.,100);
      if(postfix_preFSR.Contains("m")) FillHist("gendipt_preFSR_"+postfix_preFSR,gendy.Pt(),weight*weigtbytrig,0.,100.,100);
      if(gendy.Pt() < 100 && postfix_preFSR.Contains("m")) FillHist("gendimass_preFSR_"+postfix_preFSR,gendy.M(),weight*weigtbytrig,0.,500.,1000);

      if( fabs(genl1.Eta()) < 2.4 && fabs(genl2.Eta()) < 2.4 ){
         if( (genl1.Pt() > 25 && genl2.Pt() > 15) || (genl1.Pt() > 15 && genl2.Pt() > 25)){
            FillHist("gendipt_etaCut_25_15_preFSR",gendy.Pt(),weight*weigtbytrig,0.,100.,100);
            if(postfix_preFSR.Contains("m")) FillHist("gendipt_etaCut_25_15_preFSR_"+postfix_preFSR,gendy.Pt(),weight*weigtbytrig,0.,100.,100);
            if(gendy.Pt() < 100 && postfix_preFSR.Contains("m")) FillHist("gendimass_etaCut_25_15_preFSR_"+postfix_preFSR,gendy.M(),weight*weigtbytrig,0.,500.,1000);
         }
         if( (genl1.Pt() > 20 && genl2.Pt() > 10) || (genl1.Pt() > 10 && genl2.Pt() > 20)){
            FillHist("gendipt_etaCut_25_15_preFSR",gendy.Pt(),weight*weigtbytrig,0.,100.,100);
            if(postfix_preFSR.Contains("m")) FillHist("gendipt_etaCut_20_10_preFSR_"+postfix_preFSR,gendy.Pt(),weight*weigtbytrig,0.,100.,100);
            if(gendy.Pt() < 100 && postfix_preFSR.Contains("m")) FillHist("gendimass_etaCut_20_10_preFSR_"+postfix_preFSR,gendy.M(),weight*weigtbytrig,0.,500.,1000);
         }
      }

      TString postfix="";
      for(int i = 0; i < nBoundaries-1; i++){
         if(gendy_postFSR.M() > massBinBoundaries[i] && gendy_postFSR.M() < massBinBoundaries[i+1]){
           postfix=massbins[i];
         }
      }// loop for mass bins
      
      FillHist("gendipt_postFSR",gendy_postFSR.Pt(),weight*weigtbytrig,0.,100.,100);
      if(postfix.Contains("m")) FillHist("gendipt_postFSR_"+postfix,gendy_postFSR.Pt(),weight*weigtbytrig,0.,100.,100);
      if(gendy_postFSR.Pt() < 100) FillHist("gendimass_postFSR",gendy_postFSR.M(),weight*weigtbytrig,0.,500.,1000);
 
      if( fabs(genl1_postFSR.Eta()) < 2.4 && fabs(genl2_postFSR.Eta()) < 2.4 ){

        FillHist("gendipt_etaCut_postFSR",gendy_postFSR.Pt(),weight*weigtbytrig,0.,100.,100);
        if(postfix.Contains("m")) FillHist("gendipt_etaCut_postFSR_"+postfix,gendy_postFSR.Pt(),weight*weigtbytrig,0.,100.,100);
        if(gendy_postFSR.Pt() < 100) FillHist("gendimass_etaCut_postFSR",gendy_postFSR.M(),weight*weigtbytrig,0.,500.,1000);

         if( (genl1_postFSR.Pt() > 23 && genl2_postFSR.Pt() > 12) || (genl1_postFSR.Pt() > 12 && genl2_postFSR.Pt() > 23)){
            FillHist("gendipt_etaCut_23_12_postFSR",gendy_postFSR.Pt(),weight*weigtbytrig,0.,100.,100);
            if(postfix.Contains("m")) FillHist("gendipt_etaCut_23_12_postFSR_"+postfix,gendy_postFSR.Pt(),weight*weigtbytrig,0.,100.,100);
            if(gendy_postFSR.Pt() < 100) FillHist("gendimass_etaCut_23_12_postFSR",gendy_postFSR.M(),weight*weigtbytrig,0.,500.,1000);
         }  

         if( (genl1_postFSR.Pt() > 25 && genl2_postFSR.Pt() > 15) || (genl1_postFSR.Pt() > 15 && genl2_postFSR.Pt() > 25)){
            FillHist("gendipt_etaCut_25_15_postFSR",gendy_postFSR.Pt(),weight*weigtbytrig,0.,100.,100);
            if(postfix.Contains("m")) FillHist("gendipt_etaCut_25_15_postFSR_"+postfix,gendy_postFSR.Pt(),weight*weigtbytrig,0.,100.,100);
            if(gendy_postFSR.Pt() < 100) FillHist("gendimass_etaCut_25_15_postFSR",gendy_postFSR.M(),weight*weigtbytrig,0.,500.,1000);
         }

         if( (genl1_postFSR.Pt() > 25 && genl2_postFSR.Pt() > 25) ){
            FillHist("gendipt_etaCut_25_25_postFSR",gendy_postFSR.Pt(),weight*weigtbytrig,0.,100.,100);
            if(postfix.Contains("m")) FillHist("gendipt_etaCut_25_25_postFSR_"+postfix,gendy_postFSR.Pt(),weight*weigtbytrig,0.,100.,100);
            if(gendy_postFSR.Pt() < 100) FillHist("gendimass_etaCut_25_25_postFSR",gendy_postFSR.M(),weight*weigtbytrig,0.,500.,1000);
         }  

         if( (genl1_postFSR.Pt() > 28 && genl2_postFSR.Pt() > 17) || (genl1_postFSR.Pt() > 17 && genl2_postFSR.Pt() > 28)){
            FillHist("gendipt_etaCut_28_17_postFSR",gendy_postFSR.Pt(),weight*weigtbytrig,0.,100.,100);
            if(postfix.Contains("m")) FillHist("gendipt_etaCut_28_17_postFSR_"+postfix,gendy_postFSR.Pt(),weight*weigtbytrig,0.,100.,100);
            if(gendy_postFSR.Pt() < 100) FillHist("gendimass_etaCut_28_17_postFSR",gendy_postFSR.M(),weight*weigtbytrig,0.,500.,1000);
         }  
      }

      if(genptlep1>genptlep2){
	FillHist("genleadingpt",genptlep1,weight*weigtbytrig,0.,500.,500);
	FillHist("gensubleadingpt",genptlep2,weight*weigtbytrig,0.,500.,500);
	FillHist("genleadingeta",genetalep1,weight*weigtbytrig,-4.,4.,80);
	FillHist("gensubleadingeta",genetalep2,weight*weigtbytrig,-4.,4.,80);
      }else{
	FillHist("genleadingpt",genptlep2,weight*weigtbytrig,0.,500.,500);
	FillHist("gensubleadingpt",genptlep1,weight*weigtbytrig,0.,500.,50);
	FillHist("genleadingeta",genetalep2,weight*weigtbytrig,-4.,4.,80);
	FillHist("gensubleadingeta",genetalep1,weight*weigtbytrig,-4.,4.,80);
      }    

     for(int i = 0; i < nBoundaries-1; i++){
        if(gendimass > massBinBoundaries[i] && gendimass < massBinBoundaries[i+1]){
          FillHist("gendielectronpt_"+massbins[i],gendipt,weight*weigtbytrig,0.,100.,100);
          FillHist("gendielectronmass_"+massbins[i],gendimass,weight*weigtbytrig,0.,500.,1000);
          if(genptlep1 > genptlep2){
            if(genptlep1 > 28 && genptlep2 > 17 && fabs(genetalep1) < 2.4 && fabs(genetalep2) < 2.4){
              FillHist("gendielectronpt_kincut_"+massbins[i],gendipt,weight*weigtbytrig,0.,100.,100);
              FillHist("gendielectronmass_kincut_"+massbins[i],gendimass,weight*weigtbytrig,0.,500.,1000);
            }
          }
          if(genptlep2 > genptlep1){
            if(genptlep2 > 28 && genptlep1 > 17 && fabs(genetalep1) < 2.4 && fabs(genetalep2) < 2.4){
              FillHist("gendielectronpt_kincut_"+massbins[i],gendipt,weight*weigtbytrig,0.,100.,100);
              FillHist("gendielectronmass_kincut_"+massbins[i],gendimass,weight*weigtbytrig,0.,500.,1000);
            }
          }
        }
     }// loop for mass bins

    }// di-electron found at gen level
  } // gen information for DY MC
    
  m_logger << DEBUG << "RunNumber/Event Number = "  << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;
  m_logger << DEBUG << "isData = " << isData << LQLogger::endmsg;
   
  FillHist("Nocut_event", 1., weight, 0.,2.,2);
  FillHist("isLumiMaskGold_nocut",eventbase->GetEvent().LumiMask(), weight,0.,2.,2);

  float pileup_reweight=(1.0);
  float pileup_reweight_up=(1.0);
  float pileup_reweight_down=(1.0);
  if (!k_isdata) {
    pileup_reweight=mcdata_correction->CatPileupWeight(eventbase->GetEvent(),0);
    pileup_reweight_up=mcdata_correction->CatPileupWeight(eventbase->GetEvent(),1);
    pileup_reweight_down=mcdata_correction->CatPileupWeight(eventbase->GetEvent(),-1);
  }
  
  if(isData) FillHist("Nvtx_nocut_data",  eventbase->GetEvent().nVertices() ,weight, 0. , 50., 50);
  else{
    FillHist("Nvtx_nocut_mc",  eventbase->GetEvent().nVertices() ,weight, 0. , 50., 50);
    FillHist("Nvtx_nocut_mc_pureweight",  eventbase->GetEvent().nVertices() ,weight*pileup_reweight, 0. , 50., 50);
  }  
  if(!PassMETFilter()) return;     /// Initial event cuts : 
  if(event_m80to100) FillCutFlow("METfilter", 1.);
  FillHist("Basic_METFilter_Cut", 1, weight, 0.,2.,2);
  /// #### CAT::: triggers stored are all HLT_Ele/HLT_DoubleEle/HLT_Mu/HLT_TkMu/HLT_Photon/HLT_DoublePhoton
  
  if(!eventbase->GetEvent().HasGoodPrimaryVertex()) return; //// Make cut on event wrt vertex                                                                               
  if(event_m80to100) FillCutFlow("GoodPV", 1.);
  FillHist("NonZero_Nvtx", 1, weight, 0.,2.,2);
  
  std::vector<snu::KElectron> electrons =  GetElectrons(true,true, "ELECTRON_POG_MEDIUM"); // Cut Based POG Medium WP 
  std::vector<snu::KElectron> electrons_ip =  GetElectrons(true,true, "ELECTRON_POG_MEDIUM_IP"); // + IP cuts

   std::vector<snu::KMuon> muons =GetMuons("MUON_POG_TIGHT",true); 

   bool trig_pass= PassTriggerOR(trignames);

   double idsf = mcdata_correction->ElectronScaleFactor("ELECTRON_POG_MEDIUM", electrons, 0);
   double idsf_up = mcdata_correction->ElectronScaleFactor("ELECTRON_POG_MEDIUM", electrons, 1);
   double idsf_down = mcdata_correction->ElectronScaleFactor("ELECTRON_POG_MEDIUM", electrons, -1);

   double recosf = mcdata_correction->ElectronRecoScaleFactor(electrons, 0);
   double recosf_up = mcdata_correction->ElectronRecoScaleFactor(electrons, 1);
   double recosf_down = mcdata_correction->ElectronRecoScaleFactor(electrons, -1);

   bool is_doubleelectron = (electrons.size() == 2);

   if(trig_pass) {
     if(event_m80to100) FillCutFlow("Trigger",1.);
   }
 
   // medium working point pog id without ip cuts
   if(is_doubleelectron){
   if(event_m80to100 && trig_pass) FillCutFlow("TwoPOGEle",1.); 
   
   bool is_os = (electrons.at(0).Charge() == (-electrons.at(1).Charge()));
   FillHist("os_cut", 1, weight, 0.,2.,2);

   double ptlep1 = electrons[0].Pt();
   double ptlep2 = electrons[1].Pt();
   double etalep1 = electrons[0].Eta();
   double etalep2 = electrons[1].Eta();
   double dipt = (electrons[0]+electrons[1]).Pt();
   double dimass = (electrons[0]+electrons[1]).M();

   bool TriggerMatch1 = (electrons[0].TriggerMatched(dielectron_trig) && electrons[1].TriggerMatched(dielectron_trig));

   bool mcfromtau = (electrons[0].MCFromTau()||electrons[1].MCFromTau()); // is there case only one of the lepton decaying from tau in DY sample?
   TString prefix="";
   if(mcfromtau&&k_sample_name.Contains("DY")) prefix="tau_";

   // set histogram range
   double dimass_min = 0., dimass_max = 500.;
   int ndimass = 1000;
   double dipt_min = 0., dipt_max = 100.;
   int ndipt = 100;
   double leppt_min = 0., leppt_max = 500.;
   int npt = 500;
   double lepeta_min = -2.5, lepeta_max = 2.5;
   int neta = 50;

     if( trig_pass && TriggerMatch1){
       if(event_m80to100) FillCutFlow("trigger_match",1.);
     }

     if( is_os && trig_pass && TriggerMatch1){
       if(event_m80to100) FillCutFlow("OS",1.);
     }
   // temporary kinematic cut

   TString postfix1="_asymptcut";

   // Kinematic cuts from 2016 DY x-section
   if(trig_pass && is_os && TriggerMatch1 && ptlep1 > 25. && ptlep2 > 15. && fabs(etalep1) < 2.4 && fabs(etalep2) < 2.4 && dipt < 100.){
       if(event_m80to100) FillCutFlow("recoKincuts",1.);

     for(int i = 0; i < nBoundaries-1; i++){
        if(dimass > massBinBoundaries[i] && dimass < massBinBoundaries[i+1]){

          // dilepton pt for each mass bin
          FillHist(prefix+"dielectronpt_"+massbins[i]+postfix1,dipt,weight*idsf*pileup_reweight*recosf*weigtbytrig, dipt_min, dipt_max, ndipt);
          FillHist(prefix+"dielectronmass_"+massbins[i]+postfix1,dimass,weight*idsf*pileup_reweight*recosf*weigtbytrig, dimass_min, dimass_max, ndimass);
          // leading and subleading pt for each mass bin
          FillHist(prefix+"leadingpt_"+massbins[i]+postfix1,ptlep1,weight*idsf*pileup_reweight*recosf*weigtbytrig, leppt_min, leppt_max, npt);
          FillHist(prefix+"subleadingpt_"+massbins[i]+postfix1,ptlep2,weight*idsf*pileup_reweight*recosf*weigtbytrig, leppt_min, leppt_max, npt);

          FillHist(prefix+"leadingeta_"+massbins[i]+postfix1,etalep1,weight*idsf*pileup_reweight*recosf*weigtbytrig, lepeta_min, lepeta_max, neta);
          FillHist(prefix+"subleadingeta_"+massbins[i]+postfix1,etalep2,weight*idsf*pileup_reweight*recosf*weigtbytrig, lepeta_min, lepeta_max, neta);
        }
     }

     FillHist(prefix+"Mass"+postfix1,dimass,weight*idsf*pileup_reweight*recosf*weigtbytrig, dimass_min, dimass_max, ndimass);
     FillHist(prefix+"Pt"+postfix1,dipt,weight*idsf*pileup_reweight*recosf*weigtbytrig, dipt_min, dipt_max, ndipt);
     FillHist(prefix+"Met"+postfix1, eventbase->GetEvent().MET(),weight*idsf*pileup_reweight*recosf*weigtbytrig, dimass_min, dimass_max, ndimass);

     FillHist(prefix+"DZ"+postfix1, electrons[0].dz(), weight*idsf*pileup_reweight*recosf*weigtbytrig, -.3, .3, 300);
     FillHist(prefix+"DZ"+postfix1, electrons[1].dz(), weight*idsf*pileup_reweight*recosf*weigtbytrig, -.3, .3, 300);

     FillHist(prefix+"DXY"+postfix1, electrons[0].dxy(), weight*idsf*pileup_reweight*recosf*weigtbytrig, -0.2, 0.2, 200);
     FillHist(prefix+"DXY"+postfix1, electrons[1].dxy(), weight*idsf*pileup_reweight*recosf*weigtbytrig, -0.2, 0.2, 200);

     FillHist(prefix+"Nvtx"+postfix1,eventbase->GetEvent().nVertices() ,weight*idsf*pileup_reweight*recosf*weigtbytrig, 0. , 50., 50);
   } // event selection for opposite sign electrons 


   }// exactly two pog wp electrons

   double idwipsf = mcdata_correction->ElectronScaleFactor("ELECTRON_POG_MEDIUM", electrons_ip, 0);
   double idwipsf_up = mcdata_correction->ElectronScaleFactor("ELECTRON_POG_MEDIUM", electrons_ip, 1);
   double idwipsf_down = mcdata_correction->ElectronScaleFactor("ELECTRON_POG_MEDIUM", electrons_ip, -1);

   double recowipsf = mcdata_correction->ElectronRecoScaleFactor(electrons_ip, 0);
   double recowipsf_up = mcdata_correction->ElectronRecoScaleFactor(electrons_ip, 1);
   double recowipsf_down = mcdata_correction->ElectronRecoScaleFactor(electrons_ip, -1);

   // pog id + ip cuts
   bool is_doubleelectronip = (electrons_ip.size() == 2);
   if(is_doubleelectronip){
  
   bool is_os = (electrons_ip.at(0).Charge() == (-electrons_ip.at(1).Charge()));

   double ptlep1 = electrons_ip[0].Pt();
   double ptlep2 = electrons_ip[1].Pt();
   double etalep1 = electrons_ip[0].Eta();
   double etalep2 = electrons_ip[1].Eta();
   double dipt = (electrons_ip[0]+electrons_ip[1]).Pt();
   double dimass = (electrons_ip[0]+electrons_ip[1]).M();

   bool TriggerMatch1 = (electrons_ip[0].TriggerMatched(dielectron_trig) && electrons_ip[1].TriggerMatched(dielectron_trig));

   bool mcfromtau = (electrons_ip[0].MCFromTau()||electrons_ip[1].MCFromTau()); // is there case only one of the lepton decaying from tau in DY sample?
   TString prefix="";
   if(mcfromtau&&k_sample_name.Contains("DY")) prefix="tau_";

   // set histogram range
   double dimass_min = 0., dimass_max = 500.;
   int ndimass = 1000;
   double dipt_min = 0., dipt_max = 100.;
   int ndipt = 100;
   double leppt_min = 0., leppt_max = 500.;
   int npt = 500;
   double lepeta_min = -2.5, lepeta_max = 2.5;
   int neta = 50;

   TString postfix="_ipcut";

   // temporary kinematic cut
   if(trig_pass && is_os && TriggerMatch1 && ptlep1 > 25. && ptlep2 > 15. && fabs(etalep1) < 2.4 && fabs(etalep2) < 2.4 && dipt < 100.){

     for(int i = 0; i < nBoundaries-1; i++){
        if(dimass > massBinBoundaries[i] && dimass < massBinBoundaries[i+1]){

          // dilepton pt for each mass bin
          FillHist(prefix+"dielectronpt_"+massbins[i]+postfix,dipt,weight*idwipsf*pileup_reweight*recowipsf*weigtbytrig, dipt_min, dipt_max, ndipt);
          FillHist(prefix+"dielectronmass_"+massbins[i]+postfix,dimass,weight*idwipsf*pileup_reweight*recowipsf*weigtbytrig, dimass_min, dimass_max, ndimass);
          // leading and subleading pt for each mass bin
          FillHist(prefix+"leadingpt_"+massbins[i]+postfix,ptlep1,weight*idwipsf*pileup_reweight*recowipsf*weigtbytrig, leppt_min, leppt_max, npt);
          FillHist(prefix+"subleadingpt_"+massbins[i]+postfix,ptlep2,weight*idwipsf*pileup_reweight*recowipsf*weigtbytrig, leppt_min, leppt_max, npt);

          FillHist(prefix+"leadingeta_"+massbins[i]+postfix,etalep1,weight*idwipsf*pileup_reweight*recowipsf*weigtbytrig, lepeta_min, lepeta_max, neta);
          FillHist(prefix+"subleadingeta_"+massbins[i]+postfix,etalep2,weight*idwipsf*pileup_reweight*recowipsf*weigtbytrig, lepeta_min, lepeta_max, neta);
        }
     }
     FillHist(prefix+"Mass"+postfix,dimass,weight*idwipsf*pileup_reweight*recowipsf*weigtbytrig, dimass_min, dimass_max, ndimass);
     FillHist(prefix+"Pt"+postfix,dipt,weight*idwipsf*pileup_reweight*recowipsf*weigtbytrig, dipt_min, dipt_max, ndipt);
     FillHist(prefix+"Met"+postfix, eventbase->GetEvent().MET(),weight*idwipsf*pileup_reweight*recowipsf*weigtbytrig, dimass_min, dimass_max, ndimass);

     FillHist(prefix+"DZ"+postfix, electrons_ip[0].dz(), weight*idwipsf*pileup_reweight*recowipsf*weigtbytrig, -.3, .3, 100);
     FillHist(prefix+"DZ"+postfix, electrons_ip[1].dz(), weight*idwipsf*pileup_reweight*recowipsf*weigtbytrig, -.3, .3, 100);

     FillHist(prefix+"DXY"+postfix, electrons_ip[0].dxy(), weight*idwipsf*pileup_reweight*recowipsf*weigtbytrig, -.2, .2, 200);
     FillHist(prefix+"DXY"+postfix, electrons_ip[1].dxy(), weight*idwipsf*pileup_reweight*recowipsf*weigtbytrig, -.2, .2, 200);

     FillHist(prefix+"Nvtx"+postfix,eventbase->GetEvent().nVertices() ,weight*idsf*pileup_reweight*recosf*weigtbytrig, 0. , 50., 50);

     // check mass and pt with MET cut
     if(eventbase->GetEvent().MET() < 35.){
       FillHist(prefix+"Mass_met"+postfix,dimass,weight*idwipsf*pileup_reweight*recowipsf*weigtbytrig, dimass_min, dimass_max, ndimass);
       FillHist(prefix+"Pt_met"+postfix,dipt,weight*idwipsf*pileup_reweight*recowipsf*weigtbytrig, dipt_min, dipt_max, ndipt);
     }
   } // event selection for opposite sign electrons 
  } 
   return;
}// End of execute event loop
  


void ISRee::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void ISRee::BeginCycle() throw( LQError ){
  
  Message("In begin Cycle", INFO);
  
  //
  //If you wish to output variables to output file use DeclareVariable
  // clear these variables in ::ClearOutputVectors function
  //DeclareVariable(obj, label, treename );
  //DeclareVariable(obj, label ); //-> will use default treename: LQTree
  //  DeclareVariable(out_electrons, "Signal_Electrons", "LQTree");
  //  DeclareVariable(out_muons, "Signal_Muons");

  
  return;
  
}

ISRee::~ISRee() {
  
  Message("In ISRee Destructor" , INFO);
  
}


void ISRee::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void ISRee::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  maphist3D.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this ISReeCore::MakeHistograms() to make new hists for your analysis
   **/

}


void ISRee::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}



