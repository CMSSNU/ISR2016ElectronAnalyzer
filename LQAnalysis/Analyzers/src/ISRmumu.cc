// $Id: ISRmumu.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQISRmumu Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "ISRmumu.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (ISRmumu);


 /**
  *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
  *
  */
ISRmumu::ISRmumu() :  AnalyzerCore(), out_electrons(0)  {
  
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("ISRmumu");
  
  Message("In ISRmumu constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();
  //  MakeCleverHistograms(sighist_mm,"DiElectron");


}


void ISRmumu::InitialiseAnalysis() throw( LQError ) {
  
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


void ISRmumu::ExecuteEvents()throw( LQError ){
  /// Apply the gen weight 
  if(!isData) weight*=MCweight;
  
  m_logger << DEBUG << "RunNumber/Event Number = "  << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;
  m_logger << DEBUG << "isData = " << isData << LQLogger::endmsg;   

  // for SingleMuon trigger (not default)
  TString single_trig1 = "HLT_IsoMu24_v";
  TString single_trig2 = "HLT_IsoTkMu24_v";  
  vector<TString> trignames_single;
  trignames_single.push_back(single_trig1);
  trignames_single.push_back(single_trig2);

  // for DoubleMuon trigger (default)
  vector<TString> trignames_double;
  trignames_double.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v");
  trignames_double.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v");
  trignames_double.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
  trignames_double.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");

  double lumiweight_single=1,lumiweight_double=1;
  if(!isData){
    lumiweight_single*=WeightByTrigger(trignames_single,TargetLumi); // change MC luminocity from 1pb^-1(default) to data luminocity
    lumiweight_double*=WeightByTrigger(trignames_double,TargetLumi);
  }
  double weightGen=weight*lumiweight_double;

  FillCutFlow("NoCut", weightGen);
  // Get generator information for gen/reco ratio correction (only for Drell-Yan MC).
  // For this, we needs full-phase space generator level information.
  // So, You should run this with FLATCAT.
  bool ishardFSR=false;
  if(k_sample_name.Contains("DY")){  //only for Drell-Yan MC
    std::vector<snu::KTruth> truthcol =  eventbase->GetTruth();   //get truth particles
    TLorentzVector gendy,genl1,genl2,genl1pre,genl2pre;   //gendy: Z/gamma* fourvector, genl1: muon fourvector, genl2: anti-muon fourvector, *pre: before FSR
    int truthsize=truthcol.size();

    //loop for collect dimuon Drell-Yan product
    for(int i=0;i<truthsize;i++){
      snu::KTruth truth=truthcol.at(i);
      if(truth.GenStatus()!=1) continue;  //stable-particle-requirement
      if(truth.PdgId()==13&&truth.StatusFlag(snu::KTruth::fromhardprocess)){
	genl1+=truth;
      }else if(truth.PdgId()==-13&&truth.StatusFlag(snu::KTruth::fromhardprocess)){
	genl2+=truth;
      }else if(truth.PdgId()==22){   //collect photons from DY muons by FSR
	int imother=truth.IndexMother();
	if(imother>=truthsize) continue;
	snu::KTruth mother=truthcol.at(truth.IndexMother());
	if(mother.PdgId()==13&&mother.StatusFlag(snu::KTruth::fromhardprocess)){
	  genl1pre+=truth;
	}else if(mother.PdgId()==-13&&mother.StatusFlag(snu::KTruth::fromhardprocess)){
	  genl2pre+=truth;
	}
      }
    }
    genl1pre+=genl1;
    genl2pre+=genl2;
    gendy=genl1pre+genl2pre;

    //Fill generator level hists, if we find both of muon and anti-muon
    if(genl1!=TLorentzVector(0,0,0,0)&&genl2!=TLorentzVector(0,0,0,0)){  
      double gendimass = gendy.M();  
      double gendipt = gendy.Pt();
      double genl1pt = genl1.Pt();
      double genl2pt = genl2.Pt();
      double genl1eta = genl1.Eta();
      double genl2eta = genl2.Eta();
      double genl1mass = genl1.M();
      double genl2mass = genl2.M();
      if(genl1pt<genl2pt){
	double temppt=genl2pt; 
	double tempeta=genl2eta;
	double tempmass=genl2mass;
	genl2pt=genl1pt;
	genl2eta=genl1eta;
	genl2mass=genl1mass;
	genl1pt=temppt;
	genl1eta=tempeta;
	genl1mass=tempmass;
      }
      if(gendy.M()-(genl1+genl2).M()>5){
	ishardFSR=true;
      } 
    }
   }

  if(!PassMETFilter()) return;     /// Initial event cuts : 
  if(!eventbase->GetEvent().HasGoodPrimaryVertex()) return; //// Make cut on event wrt vertex
  FillCutFlow("BasicEventCut", weightGen);    

  std::vector<snu::KMuon> muons =GetMuons("MUON_POG_TIGHT",true); //get muon collection
  std::vector<snu::KElectron> electrons =  GetElectrons(false,false, "ELECTRON_NOCUT"); //get electron collection

  std::vector<snu::KMuon> muon_lead;
  std::vector<snu::KMuon> muon_sublead;

  //select dimuon events
  if(muons.size()!=2) return;
  //if(electrons.size()>0) return;
  FillCutFlow("dimuonCut",weightGen);

  muon_lead.push_back(muons.at(0));
  muon_sublead.push_back(muons.at(1));  

  //check whether this event pass SingleMuon trigger
  bool passtrigger_single=true;
  if(!PassTriggerOR(trignames_single)) passtrigger_single=false;
  //check whether one of dimuon fire trigger
  if(!(muons[0].TriggerMatched(single_trig1)||muons[1].TriggerMatched(single_trig1)||muons[0].TriggerMatched(single_trig2)||muons[1].TriggerMatched(single_trig2))) passtrigger_single=false;
  
  //check whether this event pass DoubleMuon trigger
  bool passtrigger_double=PassTriggerOR(trignames_double);
  //check whether dimuon fire trigger
  bool matchtrigger_double=false;
  for(int i=0;i<4;i++){
    if(muons[0].TriggerMatched(trignames_double[i])&&muons[1].TriggerMatched(trignames_double[i])) matchtrigger_double=true;
  }
  passtrigger_double&=matchtrigger_double;

  //scale factors
  double PUreweight=1,PUreweight_up=1,PUreweight_down=1;
  double IDSF=1,IDSF_up=1,IDSF_down=1;
  double IDSF_lead = 1, IDSF_sublead = 1;
  double ISOSF_lead = 1, ISOSF_sublead = 1;
  double ISOSF=1,ISOSF_up=1,ISOSF_down=1;
  double triggerSF_double=1.;
  if(!isData){
    PUreweight=mcdata_correction->CatPileupWeight(eventbase->GetEvent(),0);  //pileup reweight
    PUreweight_up=mcdata_correction->CatPileupWeight(eventbase->GetEvent(),1);  //pileup reweight sys
    PUreweight_down=mcdata_correction->CatPileupWeight(eventbase->GetEvent(),-1);  //pileup reweight sys

    IDSF=mcdata_correction->MuonScaleFactor("MUON_POG_TIGHT", muons, 0);  //muon ID scale factor
    IDSF_lead=mcdata_correction->MuonScaleFactor("MUON_POG_TIGHT", muon_lead, 0);  //muon ID scale factor
    IDSF_sublead=mcdata_correction->MuonScaleFactor("MUON_POG_TIGHT", muon_sublead, 0);  //muon ID scale factor
    IDSF_up=mcdata_correction->MuonScaleFactor("MUON_POG_TIGHT", muons, 1);  //muon ID scale factor sys
    IDSF_down=mcdata_correction->MuonScaleFactor("MUON_POG_TIGHT", muons, -1);  //muon ID scale factor sys

    ISOSF=mcdata_correction->MuonISOScaleFactor("MUON_POG_TIGHT", muons, 0);  //muon isolation scale factor
    ISOSF_lead=mcdata_correction->MuonISOScaleFactor("MUON_POG_TIGHT", muon_lead, 0);  //muon isolation scale factor
    ISOSF_sublead=mcdata_correction->MuonISOScaleFactor("MUON_POG_TIGHT", muon_sublead, 0);  //muon isolation scale factor
    ISOSF_up=mcdata_correction->MuonISOScaleFactor("MUON_POG_TIGHT", muons, 1);  //muon isolation scale factor sys
    ISOSF_down=mcdata_correction->MuonISOScaleFactor("MUON_POG_TIGHT", muons, -1);  //muon isolation scale factor sys

    triggerSF_double=mcdata_correction->GetDoubleMUTriggerEffISR(muons, 0);  //no scale factor for double muon trigger yet
  }
  double weightTotal=weight*lumiweight_double*PUreweight*IDSF*ISOSF;
  FillCutFlow("AfterScaleFactor",weightTotal);

  //select opposite sign events
  if(muons.at(0).Charge()==muons.at(1).Charge()) return;
  FillCutFlow("osCut",weightTotal);
  
  //relative isolation cut
  if(muons.at(0).RelIso04() > 0.15) return;
  if(muons.at(1).RelIso04() > 0.15) return;
  FillCutFlow("RisoCut",weightTotal);

  CorrectMuonMomentum(muons);   //apply rochester correction
  CorrectedMETRochester(muons);  //update MET after rochester correction
  
  double dimass = (muons[0]+muons[1]).M();  
  double dipt = (muons[0]+muons[1]).Pt();
  double dieta = (muons[0]+muons[1]).Eta();
  double met=eventbase->GetEvent().MET();
  int nvtx=eventbase->GetEvent().nVertices();
  double l1pt = muons[0].Pt();
  double l2pt = muons[1].Pt();
  double l1eta = muons[0].Eta();
  double l2eta = muons[1].Eta();
  double l1mass = muons[0].M();
  double l2mass = muons[1].M();

  if(dimass<15) return;
  FillCutFlow("MassCut",weightTotal);

  //missing E_T cut
  //if(met>35) return;
  //FillCutFlow("METCut",weightTotal);

  //mark 'DY -> tau tau' events
  bool mcfromtau = (muons[0].MCFromTau()||muons[1].MCFromTau());
  TString prefix="";
  if(mcfromtau&&k_sample_name.Contains("DY")) prefix="tau_";

  if(l1pt<l2pt){
    double temppt=l2pt; 
    double tempeta=l2eta;
    double tempmass=l2mass;
    l2pt=l1pt;
    l2eta=l1eta;
    l2mass=l1mass;
    l1pt=temppt;
    l1eta=tempeta;
    l1mass=tempmass;
  }

   // set histogram range
   double dimass_min = 0., dimass_max = 500.;
   int ndimass = 1000;
   double dipt_min = 0., dipt_max = 100.;
   int ndipt = 100;
   double leppt_min = 0., leppt_max = 500.;
   int npt = 500;
   double lepeta_min = -2.5, lepeta_max = 2.5;
   int neta = 50;


  if(passtrigger_double){
    //////////////Default/////////////////////
    if(l1pt>20&&l2pt>10&&(fabs(l1eta)<2.4)&&(fabs(l2eta)<2.4) && dipt < 100.){
      FillCutFlow("PtEtaCut",weightTotal);
      TString postfix1="_asymptcut";

      if(dimass > 40. && dimass < 350.){
        // dilepton pt for each mass bin, without any scale factor
        FillHist(prefix+"dimuonptNoRecoSFNoTrgSFNoIDSF_"+postfix1,dipt,weight*lumiweight_double*PUreweight_up, dipt_min, dipt_max, ndipt);
        FillHist(prefix+"dimuonmassNoRecoSFNoTrgSFNoIDSF_"+postfix1,dimass,weight*lumiweight_double*PUreweight_up, dimass_min, dimass_max, ndimass);

        // dilepton pt for each mass bin
        FillHist(prefix+"dimuonptRecoSFNoTrgSFNoIDSF_"+postfix1,dipt,weight*lumiweight_double*PUreweight_up*IDSF, dipt_min, dipt_max, ndipt);
        FillHist(prefix+"dimuonmassRecoSFNoTrgSFNoIDSF_"+postfix1,dimass,weight*lumiweight_double*PUreweight_up*IDSF, dimass_min, dimass_max, ndimass);

        // dilepton pt for each mass bin, No trigger scale factor applied 
        FillHist(prefix+"dimuonptRecoSFNoTrgSFIDSF_"+postfix1,dipt,weight*lumiweight_double*PUreweight_up*IDSF*ISOSF, dipt_min, dipt_max, ndipt);
        FillHist(prefix+"dimuonmassRecoSFNoTrgSFIDSF_"+postfix1,dimass,weight*lumiweight_double*PUreweight_up*IDSF*ISOSF, dimass_min, dimass_max, ndimass);

        // dilepton pt for each mass bin, trigger and ID scale factor applied
        FillHist(prefix+"dimuonptRecoSFTrgSFIDSF_"+postfix1,dipt,weight*lumiweight_double*PUreweight_up*IDSF*ISOSF*triggerSF_double, dipt_min, dipt_max, ndipt);
        FillHist(prefix+"dimuonmassRecoSFTrgSFIDSF_"+postfix1,dimass,weight*lumiweight_double*PUreweight_up*IDSF*ISOSF*triggerSF_double, dimass_min, dimass_max, ndimass);


        // leading and subleading pt, without any SF 
        FillHist(prefix+"leadingptNoRecoSFNoTrgSFNoIDSF_"+postfix1,   l1pt,weight*lumiweight_double*PUreweight_up, leppt_min, leppt_max, npt);
        FillHist(prefix+"subleadingptNoRecoSFNoTrgSFNoIDSF_"+postfix1,l2pt,weight*lumiweight_double*PUreweight_up, leppt_min, leppt_max, npt);
        // leading and subleading pt, Reco SF
        FillHist(prefix+"leadingptRecoSFNoTrgSFNoIDSF_"+postfix1,   l1pt,weight*lumiweight_double*PUreweight_up*IDSF_lead, leppt_min, leppt_max, npt);
        FillHist(prefix+"subleadingptRecoSFNoTrgSFNoIDSF_"+postfix1,l2pt,weight*lumiweight_double*PUreweight_up*IDSF_sublead, leppt_min, leppt_max, npt);
        // leading and subleading pt, ID SF 
        FillHist(prefix+"leadingptRecoSFNoTrgSFIDSF_"+postfix1,   l1pt,weight*lumiweight_double*PUreweight_up*IDSF_lead*ISOSF_lead, leppt_min, leppt_max, npt);
        FillHist(prefix+"subleadingptRecoSFNoTrgSFIDSF_"+postfix1,l2pt,weight*lumiweight_double*PUreweight_up*IDSF_sublead*ISOSF_sublead, leppt_min, leppt_max, npt);
        // leading and subleading pt, Trig SF 
        FillHist(prefix+"leadingptRecoSFTrgSFIDSF_"+postfix1,   l1pt,weight*lumiweight_double*PUreweight_up*IDSF_lead*ISOSF_lead*triggerSF_double, leppt_min, leppt_max, npt);
        FillHist(prefix+"subleadingptRecoSFTrgSFIDSF_"+postfix1,l2pt,weight*lumiweight_double*PUreweight_up*IDSF_sublead*ISOSF_sublead*triggerSF_double, leppt_min, leppt_max, npt);

        // leading and subleading eta
        FillHist(prefix+"leadingetaNoRecoSFNoTrgSFNoIDSF_"+postfix1,   l1eta,weight*lumiweight_double*PUreweight_up, lepeta_min, lepeta_max, neta);
        FillHist(prefix+"subleadingetaNoRecoSFNoTrgSFNoIDSF_"+postfix1,l2eta,weight*lumiweight_double*PUreweight_up, lepeta_min, lepeta_max, neta);
        // leading and subleading eta
        FillHist(prefix+"leadingetaRecoSFNoTrgSFNoIDSF_"+postfix1,   l1eta,weight*lumiweight_double*PUreweight_up*IDSF_lead, lepeta_min, lepeta_max, neta);
        FillHist(prefix+"subleadingetaRecoSFNoTrgSFNoIDSF_"+postfix1,l2eta,weight*lumiweight_double*PUreweight_up*IDSF_sublead, lepeta_min, lepeta_max, neta);
        // leading and subleading eta
        FillHist(prefix+"leadingetaRecoSFNoTrgSFIDSF_"+postfix1,   l1eta,weight*lumiweight_double*PUreweight_up*IDSF_lead*ISOSF_lead, lepeta_min, lepeta_max, neta);
        FillHist(prefix+"subleadingetaRecoSFNoTrgSFIDSF_"+postfix1,l2eta,weight*lumiweight_double*PUreweight_up*IDSF_sublead*ISOSF_sublead, lepeta_min, lepeta_max, neta);
        // leading and subleading eta
        FillHist(prefix+"leadingetaRecoSFTrgSFIDSF_"+postfix1,   l1eta,weight*lumiweight_double*PUreweight_up*IDSF_lead*ISOSF_lead*triggerSF_double, lepeta_min, lepeta_max, neta);
        FillHist(prefix+"subleadingetaRecoSFTrgSFIDSF_"+postfix1,l2eta,weight*lumiweight_double*PUreweight_up*IDSF_sublead*ISOSF_sublead*triggerSF_double, lepeta_min, lepeta_max, neta);

      }
    }
  }


   return;
}// End of execute event loop
  


void ISRmumu::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void ISRmumu::BeginCycle() throw( LQError ){
  
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

ISRmumu::~ISRmumu() {
  
  Message("In ISRmumu Destructor" , INFO);
  
}


void ISRmumu::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void ISRmumu::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  maphist3D.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this ISRmumuCore::MakeHistograms() to make new hists for your analysis
   **/

}


void ISRmumu::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}



