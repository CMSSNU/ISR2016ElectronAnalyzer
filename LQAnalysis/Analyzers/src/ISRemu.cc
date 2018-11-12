// $Id: ISRemu.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQISRemu Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "ISRemu.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (ISRemu);


 /**
  *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
  *
  */
ISRemu::ISRemu() :  AnalyzerCore(), out_electrons(0)  {
  
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("ISRemu");
  
  Message("In ISRemu constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();
  //  MakeCleverHistograms(sighist_mm,"DiElectron");


}


void ISRemu::InitialiseAnalysis() throw( LQError ) {
  
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


void ISRemu::ExecuteEvents()throw( LQError ){

  vector<TString> trignames;

  bool event_m80to100 = false;
  double weigtbytrig = 1.;

  // use singlemuon trigger
  TString single_trig1 = "HLT_IsoMu24_v";
  TString single_trig2 = "HLT_IsoTkMu24_v";  
  trignames.push_back(single_trig1);
  trignames.push_back(single_trig2);
  if(!isData) weigtbytrig=WeightByTrigger(trignames,TargetLumi); // change MC luminocity from 1pb^-1(default) to data luminocity

  // Mass boundaries  
  double massBinBoundaries[] = {40., 60., 80., 100., 200., 350};
  const int nBoundaries = sizeof(massBinBoundaries)/sizeof(double); 

  vector<TString> massbins;
  massbins.push_back("m40to60");
  massbins.push_back("m60to80");
  massbins.push_back("m80to100");
  massbins.push_back("m100to200");
  massbins.push_back("m200to350");

  /// Apply the gen weight 
  if(!isData) weight*=MCweight;

  m_logger << DEBUG << "RunNumber/Event Number = "  << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;
  m_logger << DEBUG << "isData = " << isData << LQLogger::endmsg;
   
  FillHist("Nocut_event", 1., weight, 0.,2.,2);
  FillHist("isLumiMaskGold_nocut",eventbase->GetEvent().LumiMask(), weight,0.,2.,2);

  float pileup_reweight=(1.0);
  if (!k_isdata) {
    pileup_reweight=mcdata_correction->CatPileupWeight(eventbase->GetEvent(),0);
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
  std::vector<snu::KMuon> muons =GetMuons("MUON_POG_TIGHT",true); 

  bool trig_pass= PassTriggerOR(trignames);

  double idsf = mcdata_correction->ElectronScaleFactor("ELECTRON_POG_MEDIUM", electrons, 0);
  double recosf = mcdata_correction->ElectronRecoScaleFactor(electrons, 0);

  bool is_emuEvent = (electrons.size() > 0 && muons.size() > 0);

  if(trig_pass) {
    if(event_m80to100) FillCutFlow("Trigger",1.);
  }

    // for b veto
    std::vector<snu::KJet> jetsPreColl = GetJets("JET_NOLEPTONVETO", 20, 2.4);
    std::vector<snu::KJet> bjets_loose, bjets_medium, bjets_tight;

    for(int i = 0; i < jetsPreColl.size(); i++){
        if( jetsPreColl.at(i).IsBTagged(snu::KJet::CSVv2, snu::KJet::Loose) ) bjets_loose.push_back(jetsPreColl.at(i));
        if( jetsPreColl.at(i).IsBTagged(snu::KJet::CSVv2, snu::KJet::Medium) ) bjets_medium.push_back(jetsPreColl.at(i));
        if( jetsPreColl.at(i).IsBTagged(snu::KJet::CSVv2, snu::KJet::Tight) ) bjets_tight.push_back(jetsPreColl.at(i));
    }

 
  if(is_emuEvent && trig_pass){
    if ((electrons.at(0).Charge() == (-muons.at(0).Charge())) ){ // opposite sign emu
       if(electrons.at(0).Pt() > 26 && muons.at(0).Pt() > 26 && fabs(electrons.at(0).Eta()) < 2.4 && fabs(muons.at(0).Eta()) < 2.4){
         TLorentzVector emu = electrons.at(0) + muons.at(0);
         double dimass = emu.M();
         double dipt = emu.Pt();
         if(emu.Pt() < 100.){

           bool mcfromtau = (electrons[0].MCFromTau()||muons[0].MCFromTau()); // is there case only one of the lepton decaying from tau in DY sample?
           TString prefix="";
           if(mcfromtau&&k_sample_name.Contains("DY")) prefix="tau_";
           TString postfix1="_asymptcut";

           // set histogram range
           double dimass_min = 0., dimass_max = 500.;
           int ndimass = 1000;
           double dipt_min = 0., dipt_max = 100.;
           int ndipt = 100;
           double leppt_min = 0., leppt_max = 500.;
           int npt = 500;
           double lepeta_min = -2.5, lepeta_max = 2.5;
           int neta = 50;


           for(int i = 0; i < nBoundaries-1; i++){
              if(dimass > massBinBoundaries[i] && dimass < massBinBoundaries[i+1]){
    
                // dilepton pt for each mass bin
                FillHist(prefix+"dielectronpt_"+massbins[i]+postfix1,dipt,weight*weigtbytrig, dipt_min, dipt_max, ndipt);
                FillHist(prefix+"dielectronmass_"+massbins[i]+postfix1,dimass,weight*weigtbytrig, dimass_min, dimass_max, ndimass);

                // medium B
                if(bjets_medium.size()==0){
                FillHist(prefix+"dielectronpt_"+massbins[i]+postfix1+"_mediumBveto",dipt,weight*weigtbytrig, dipt_min, dipt_max, ndipt);
                FillHist(prefix+"dielectronmass_"+massbins[i]+postfix1+"_mediumBveto",dimass,weight*weigtbytrig, dimass_min, dimass_max, ndimass);
                }

                if(bjets_medium.size()>=1){
                FillHist(prefix+"dielectronpt_"+massbins[i]+postfix1+"_mediumBcontrol",dipt,weight*weigtbytrig, dipt_min, dipt_max, ndipt);
                FillHist(prefix+"dielectronmass_"+massbins[i]+postfix1+"_mediumBcontrol",dimass,weight*weigtbytrig, dimass_min, dimass_max, ndimass);
                }

              }
           }


         }// dilepton pt cut

       }
    }

  } // event selection for opposite sign electrons 

  return;
}// End of execute event loop
  


void ISRemu::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void ISRemu::BeginCycle() throw( LQError ){
  
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

ISRemu::~ISRemu() {
  
  Message("In ISRemu Destructor" , INFO);
  
}


void ISRemu::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void ISRemu::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  maphist3D.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this ISRemuCore::MakeHistograms() to make new hists for your analysis
   **/

}


void ISRemu::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}



