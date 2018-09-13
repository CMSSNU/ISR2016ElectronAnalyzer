// $Id: DYee.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQDYee Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "DYee.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (DYee);


 /**
  *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
  *
  */
DYee::DYee() :  AnalyzerCore(), out_electrons(0)  {
  
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("DYee");
  
  Message("In DYee constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();
  //  MakeCleverHistograms(sighist_mm,"DiElectron");


}


void DYee::InitialiseAnalysis() throw( LQError ) {
  
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


void DYee::ExecuteEvents()throw( LQError ){

  // Mass boundaries  
  double massBinBoundaries[] = {40., 60., 80., 100., 200., 350};
  const int nBoundaries = sizeof(massBinBoundaries)/sizeof(double); 

  vector<TString> massbins;
  massbins.push_back("m40to60");
  massbins.push_back("m60to80");
  massbins.push_back("m80to100");
  massbins.push_back("m100to200");
  massbins.push_back("m200to350");
  massbins.push_back("m350to2000");

  /// Apply the gen weight 
  if(!isData) weight*=MCweight;

  // Get generator information for gen/reco ratio correction (only for Drell-Yan MC).
  // For this, we needs full-phase space generator level information.
  // So, You should run this with FLATCAT.
  if(k_sample_name.Contains("DY")){  //only for Drell-Yan MC
    std::vector<snu::KTruth> truthcol =  eventbase->GetTruth();   //get truth particles
    TLorentzVector gendy,genl1,genl2;   //gendy: Z/gamma* fourvector, genl1: electron fourvector, genl2: anti-electron fourvector
    int truthsize=truthcol.size();

    //loop for collect dimuon Drell-Yan product
    for(int i=0;i<truthsize;i++){
      snu::KTruth truth=truthcol.at(i);
      if(truth.GenStatus()!=1) continue;  //stable-particle-requirement
      if(truth.PdgId()==11&&truth.StatusFlag(snu::KTruth::fromhardprocess)){
	genl1+=truth;
	gendy+=truth;
      }else if(truth.PdgId()==-11&&truth.StatusFlag(snu::KTruth::fromhardprocess)){
	genl2+=truth;
	gendy+=truth;
      }else if(truth.PdgId()==22){   //collect photons from DY muons by FSR
	int imother=truth.IndexMother();
	if(imother>=truthsize) continue;
	snu::KTruth mother=truthcol.at(truth.IndexMother());
	if(abs(mother.PdgId())==11&&mother.StatusFlag(snu::KTruth::fromhardprocess)){
	  gendy+=truth;
	}
      }
    }
   
    //Fill generator level hists, if we find both electron  and anti-electron
    if(genl1!=TLorentzVector(0,0,0,0)&&genl2!=TLorentzVector(0,0,0,0)){  
      std::cout << "found gen z" << std::endl;
      double genptlep1 = genl1.Pt();
      double genptlep2 = genl2.Pt();
      double genetalep1 = genl1.Eta();
      double genetalep2 = genl2.Eta();
      double gendipt = gendy.Pt();
      double gendimass = gendy.M();  
      
      FillHist("gendipt",gendipt,weight,0.,100.,100);
      FillHist("gendimass",gendimass,weight,0.,500.,1000);
      if(genptlep1>genptlep2){
	FillHist("genleadingpt",genptlep1,weight,0.,500.,500);
	FillHist("gensubleadingpt",genptlep2,weight,0.,500.,500);
	FillHist("genleadingeta",genetalep1,weight,-4.,4.,80);
	FillHist("gensubleadingeta",genetalep2,weight,-4.,4.,80);
      }else{
	FillHist("genleadingpt",genptlep2,weight,0.,500.,500);
	FillHist("gensubleadingpt",genptlep1,weight,0.,500.,50);
	FillHist("genleadingeta",genetalep2,weight,-4.,4.,80);
	FillHist("gensubleadingeta",genetalep1,weight,-4.,4.,80);
      }    

     for(int i = 0; i < nBoundaries-1; i++){
        if(gendimass > massBinBoundaries[i] && gendimass < massBinBoundaries[i+1]){
          FillHist("gendielectronpt_"+massbins[i],gendipt,weight,0.,100.,100);
          FillHist("gendielectronmass_"+massbins[i],gendimass,weight,0.,500.,1000);
        }
     }// loop for mass bins
    }

  } // gen information for DY MC
    
  m_logger << DEBUG << "RunNumber/Event Number = "  << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;
  m_logger << DEBUG << "isData = " << isData << LQLogger::endmsg;
   
  FillCutFlow("NoCut", weight);
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
  FillCutFlow("EventCut", weight); 
  FillHist("Basic_METFilter_Cut", 1, weight, 0.,2.,2);
  /// #### CAT::: triggers stored are all HLT_Ele/HLT_DoubleEle/HLT_Mu/HLT_TkMu/HLT_Photon/HLT_DoublePhoton
  
  if (!eventbase->GetEvent().HasGoodPrimaryVertex()) return; //// Make cut on event wrt vertex                                                                               
  FillHist("NonZero_Nvtx", 1, weight, 0.,2.,2);
  
  TString dielectron_trig="HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v";//
  
  // Now you should do an OR of 4 triggers 
  vector<TString> trignames;
  trignames.push_back(dielectron_trig);

  FillHist("TargetLumi", TargetLumi, weight, 0.,5000.,50000);

  std::vector<snu::KElectron> electrons =  GetElectrons(true,true, "ELECTRON_POG_MEDIUM");

   std::vector<snu::KMuon> muons =GetMuons("MUON_POG_TIGHT",true); 

   bool trig_pass= PassTriggerOR(trignames);

   /// CorrectMuonMomentum(muons);  will also work as Funcion in AnalyzerCore just calls mcdata_correction function
   
   double idsf = mcdata_correction->ElectronScaleFactor("ELECTRON_POG_MEDIUM", electrons, 0);
   double idsf_up = mcdata_correction->ElectronScaleFactor("ELECTRON_POG_MEDIUM", electrons, 1);
   double idsf_down = mcdata_correction->ElectronScaleFactor("ELECTRON_POG_MEDIUM", electrons, -1);
   
   double triggersf =mcdata_correction->GetDoubleEGTriggerEff(electrons); // FIXME need to validate

   double recosf = mcdata_correction->ElectronRecoScaleFactor(electrons, 0);
   double recosf_up = mcdata_correction->ElectronRecoScaleFactor(electrons, 1);
   double recosf_down = mcdata_correction->ElectronRecoScaleFactor(electrons, -1);

   bool is_doubleelectron = (electrons.size() == 2);
   if(!is_doubleelectron) return;
   
   bool is_os = (electrons.at(0).Charge() == (-electrons.at(1).Charge()));
   FillHist("os_cut", 1, weight, 0.,2.,2);

   double ptlep1 = electrons[0].Pt();
   double ptlep2 = electrons[1].Pt();
   double etalep1 = electrons[0].Eta();
   double etalep2 = electrons[1].Eta();
   double dipt = (electrons[0]+electrons[1]).Pt();
   double dimass = (electrons[0]+electrons[1]).M();

   bool TriggerMatch1 = (electrons[0].TriggerMatched(dielectron_trig) && electrons[1].TriggerMatched(dielectron_trig));
   double weigtbytrig = WeightByTrigger(dielectron_trig,TargetLumi);

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

   // temporary kinematic cut
   if(trig_pass && is_os && TriggerMatch1 && ptlep1 > 25. && ptlep2 > 25. && fabs(etalep1) < 2.4 && fabs(etalep2) < 2.4 && dipt < 100.){

     for(int i = 0; i < nBoundaries-1; i++){
        if(dimass > massBinBoundaries[i] && dimass < massBinBoundaries[i+1]){

          // dilepton pt for each mass bin
          FillHist(prefix+"dielectronpt_"+massbins[i],dipt,weight*idsf*pileup_reweight*triggersf*recosf*weigtbytrig, dipt_min, dipt_max, ndipt);
          FillHist(prefix+"dielectronmass_"+massbins[i],dimass,weight*idsf*pileup_reweight*triggersf*recosf*weigtbytrig, dimass_min, dimass_max, ndimass);
          // leading and subleading pt for each mass bin
          FillHist(prefix+"leadingpt_"+massbins[i],ptlep1,weight*idsf*pileup_reweight*triggersf*recosf*weigtbytrig, leppt_min, leppt_max, npt);
          FillHist(prefix+"subleadingpt_"+massbins[i],ptlep2,weight*idsf*pileup_reweight*triggersf*recosf*weigtbytrig, leppt_min, leppt_max, npt);

          FillHist(prefix+"leadingeta_"+massbins[i],etalep1,weight*idsf*pileup_reweight*triggersf*recosf*weigtbytrig, lepeta_min, lepeta_max, neta);
          FillHist(prefix+"subleadingeta_"+massbins[i],etalep2,weight*idsf*pileup_reweight*triggersf*recosf*weigtbytrig, lepeta_min, lepeta_max, neta);
        }
     }

     FillHist(prefix+"Mass",dimass,weight*idsf*pileup_reweight*triggersf*recosf*weigtbytrig, dimass_min, dimass_max, ndimass);
     FillHist(prefix+"Pt",dipt,weight*idsf*pileup_reweight*triggersf*recosf*weigtbytrig, dipt_min, dipt_max, ndipt);
     FillHist(prefix+"Met", eventbase->GetEvent().MET(),weight*idsf*pileup_reweight*triggersf*recosf*weigtbytrig, dimass_min, dimass_max, ndimass);

     // check mass and pt with MET cut
     if(eventbase->GetEvent().MET() < 35.){
       FillHist(prefix+"Mass_met",dimass,weight*idsf*pileup_reweight*triggersf*recosf*weigtbytrig, dimass_min, dimass_max, ndimass);
       FillHist(prefix+"Pt_met",dipt,weight*idsf*pileup_reweight*triggersf*recosf*weigtbytrig, dipt_min, dipt_max, ndipt);
     }
   } // event selection for opposite sign electrons 

   TString postfix="_nokincut";

   // histograms without kincut
   if(trig_pass && is_os && TriggerMatch1 && dipt < 100.){

     for(int i = 0; i < nBoundaries-1; i++){
        if(dimass > massBinBoundaries[i] && dimass < massBinBoundaries[i+1]){

          // dilepton pt for each mass bin
          FillHist(prefix+"dielectronpt_"+massbins[i]+postfix,dipt,weight*idsf*pileup_reweight*triggersf*recosf*weigtbytrig, dipt_min, dipt_max, ndipt);
          FillHist(prefix+"dielectronmass_"+massbins[i]+postfix,dimass,weight*idsf*pileup_reweight*triggersf*recosf*weigtbytrig, dimass_min, dimass_max, ndimass);
          // leading and subleading pt for each mass bin
          FillHist(prefix+"leadingpt_"+massbins[i]+postfix,ptlep1,weight*idsf*pileup_reweight*triggersf*recosf*weigtbytrig, leppt_min, leppt_max, npt);
          FillHist(prefix+"subleadingpt_"+massbins[i]+postfix,ptlep2,weight*idsf*pileup_reweight*triggersf*recosf*weigtbytrig, leppt_min, leppt_max, npt);

          FillHist(prefix+"leadingeta_"+massbins[i]+postfix,etalep1,weight*idsf*pileup_reweight*triggersf*recosf*weigtbytrig, lepeta_min, lepeta_max, neta);
          FillHist(prefix+"subleadingeta_"+massbins[i]+postfix,etalep2,weight*idsf*pileup_reweight*triggersf*recosf*weigtbytrig, lepeta_min, lepeta_max, neta);

        }
     }

     FillHist(prefix+"Mass"+postfix,dimass,weight*idsf*pileup_reweight*triggersf*recosf*weigtbytrig, dimass_min, dimass_max, ndimass);
     FillHist(prefix+"Pt"+postfix,dipt,weight*idsf*pileup_reweight*triggersf*recosf*weigtbytrig, dipt_min, dipt_max, ndipt);
     FillHist(prefix+"Met"+postfix, eventbase->GetEvent().MET(),weight*idsf*pileup_reweight*triggersf*recosf*weigtbytrig, dimass_min, dimass_max, ndimass);


   } // event selection for opposite sign electrons 


  
   return;
}// End of execute event loop
  


void DYee::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void DYee::BeginCycle() throw( LQError ){
  
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

DYee::~DYee() {
  
  Message("In DYee Destructor" , INFO);
  
}


void DYee::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void DYee::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this DYeeCore::MakeHistograms() to make new hists for your analysis
   **/

}


void DYee::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}



