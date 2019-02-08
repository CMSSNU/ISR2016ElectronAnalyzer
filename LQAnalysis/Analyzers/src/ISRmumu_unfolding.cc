// $Id: ISRmumu_unfolding.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQISRmumu_unfolding Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "ISRmumu_unfolding.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (ISRmumu_unfolding);


 /**
  *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
  *
  */
ISRmumu_unfolding::ISRmumu_unfolding() :  AnalyzerCore(), out_electrons(0)  {
  
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("ISRmumu_unfolding");
  
  Message("In ISRmumu_unfolding constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();
  //  MakeCleverHistograms(sighist_mm,"DiElectron");


}


void ISRmumu_unfolding::InitialiseAnalysis() throw( LQError ) {
  
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


void ISRmumu_unfolding::ExecuteEvents()throw( LQError ){
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
  weightGen=weight*lumiweight_double;

  FillCutFlow("NoCut", weightGen);
  // Get generator information for gen/reco ratio correction (only for Drell-Yan MC).
  // For this, we needs full-phase space generator level information.
  // So, You should run this with FLATCAT.

 if(k_sample_name.Contains("DY")){  //only for Drell-Yan MC
   std::vector<snu::KTruth> truthcol =  eventbase->GetTruth();   //get truth particles
   TLorentzVector gendy,genl1,genl2,genl1pre,genl2pre;   //gendy: Z/gamma* fourvector, genl1: electron fourvector, genl2: anti-electron fourvector
   int truthsize=truthcol.size();

   //loop for collect dielectron Drell-Yan product
   for(int i=0;i<truthsize;i++){
     snu::KTruth truth=truthcol.at(i);
     if(truth.GenStatus()!=1) continue;  //stable-particle-requirement
     if(truth.PdgId()==13&&truth.StatusFlag(snu::KTruth::fromhardprocess)){
       genl1+=truth;
     }
     else if(truth.PdgId()==-13&&truth.StatusFlag(snu::KTruth::fromhardprocess)){
       genl2+=truth;
     }
     else if(truth.PdgId()==22){   //collect photons from DY muons by FSR
       int imother=truth.IndexMother();
       if(imother>=truthsize) continue;
       snu::KTruth mother=truthcol.at(truth.IndexMother());
       if(mother.PdgId()==13&&mother.StatusFlag(snu::KTruth::fromhardprocess)){
         genl1pre+=truth; //collect photons first, then later add to post fsr lepton i.e., gen1 or gen2
       }else if(mother.PdgId()==-13&&mother.StatusFlag(snu::KTruth::fromhardprocess)){
         genl2pre+=truth;
       }
     }
   }
   genl1pre+=genl1;
   genl2pre+=genl2;
   gendy=genl1pre+genl2pre;

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

      ///////for unfolding///////
      ptPreFSR.push_back(genl1pre.Pt());
      ptPreFSR.push_back(genl2pre.Pt());
      ptPreFSR.push_back(gendipt);
      mPreFSR.push_back(genl1pre.M());
      mPreFSR.push_back(genl2pre.M());
      mPreFSR.push_back(gendimass);
      ptGen.push_back(genl1pt);
      ptGen.push_back(genl2pt);
      ptGen.push_back((genl1+genl2).Pt());
      mGen.push_back(0.1057);
      mGen.push_back(0.1057);
      mGen.push_back((genl1+genl2).M());

      if((genl1pre.Pt()>20||genl2pre.Pt()>20)&&(genl1pre.Pt()>10&&genl2pre.Pt()>10)&&fabs(genl1pre.Eta())<2.4&&fabs(genl2pre.Eta())<2.4){
        isfiducialPreFSR = 1;
      }
      if(genl1pt>25&&genl2pt>10&&fabs(genl1eta)<2.4&&fabs(genl2eta)<2.4){
        isfiducialGen = 1;
      }
    }
  }
  float pileup_reweight=(1.0);
  float pileup_reweight_up=(1.0);
  float pileup_reweight_down=(1.0);

  if (!k_isdata) {
    pileup_reweight=mcdata_correction->CatPileupWeight(eventbase->GetEvent(),0);
    pileup_reweight_up=mcdata_correction->CatPileupWeight(eventbase->GetEvent(),1);
    pileup_reweight_down=mcdata_correction->CatPileupWeight(eventbase->GetEvent(),-1);

    weightGen *= pileup_reweight;

    // Pileup systematic
    weightGenPileUp = weightGen*pileup_reweight_up/pileup_reweight;
    weightGenPileDown = weightGen*pileup_reweight_down/pileup_reweight;

    // Scale variation
    vector<Float_t> scaleweight=eventbase->GetEvent().ScaleWeights();
    int imax=scaleweight.size();
    for(int i=0;i<imax;i++){
       weightGenScale.push_back(weightGen*scaleweight.at(i));
       //std::cout << i << " th scale: " << scaleweight.at(i) << std::endl;
    }
    vector<Float_t> pdfweight=eventbase->GetEvent().PdfWeights();
    imax=pdfweight.size();
    for(int i=0; i<imax; i++){
       weightGenPdf.push_back(weightGen*pdfweight.at(i));
    }

  }

  //
  bool passMETFilter = false;
  if(PassMETFilter()) passMETFilter = true;     

  bool hasGoodPrimaryVertex = false;
  if(eventbase->GetEvent().HasGoodPrimaryVertex()) hasGoodPrimaryVertex = true; //// Make cut on event wrt vertex

  std::vector<snu::KMuon> muons =GetMuons("MUON_POG_TIGHT",true); //get muon collection
  std::vector<snu::KElectron> electrons =  GetElectrons(false,false, "ELECTRON_NOCUT"); //get electron collection

  //select dimuon events
  bool is_doubleelectron = (muons.size()==2);

  if(is_doubleelectron){

    //check whether this event pass DoubleMuon trigger
    bool passtrigger_double=PassTriggerOR(trignames_double);
    //check whether dimuon fire trigger
    bool matchtrigger_double=false;
    for(int i=0;i<4;i++){
      if(muons[0].TriggerMatched(trignames_double[i])&&muons[1].TriggerMatched(trignames_double[i])) matchtrigger_double=true;
    }
    passtrigger_double&=matchtrigger_double;

    std::vector<snu::KMuon> muon_lead;
    std::vector<snu::KMuon> muon_sublead;

    //select opposite sign events
    bool is_os = (muons.at(0).Charge()!=muons.at(1).Charge());

    //relative isolation cut
    bool is_isolated = ((muons.at(0).RelIso04() < 0.15) && (muons.at(1).RelIso04() < 0.15));

    CorrectMuonMomentum(muons);   //apply rochester correction
    CorrectedMETRochester(muons);  //update MET after rochester correction
 
    double dimass = (muons[0]+muons[1]).M();
    double dipt = (muons[0]+muons[1]).Pt();
    double dieta = (muons[0]+muons[1]).Eta();
    double diphi = (muons[0]+muons[1]).Phi();
    double l1pt = muons[0].Pt();
    double l2pt = muons[1].Pt();
    double l1eta = muons[0].Eta();
    double l1phi = muons[0].Phi();
    double l2eta = muons[1].Eta();
    double l2phi = muons[1].Phi();
    double l1mass = muons[0].M();
    double l2mass = muons[1].M();


    //mark 'DY -> tau tau' events
    bool mcfromtau = (muons[0].MCFromTau()||muons[1].MCFromTau());
    if(mcfromtau&&k_sample_name.Contains("DY")) DYtautau = 1;;

    if(l1pt<l2pt){
      double temppt=l2pt;
      double tempeta=l2eta;
      double tempphi=l2phi;
      double tempmass=l2mass;
      l2pt=l1pt;
      l2eta=l1eta;
      l2phi=l1phi;
      l2mass=l1mass;

      l1pt=temppt;
      l1eta=tempeta;
      l1phi=tempphi;
      l1mass=tempmass;

      muon_lead.push_back(muons.at(1));
      muon_sublead.push_back(muons.at(0));
    }
    else{
         muon_lead.push_back(muons.at(0));
         muon_sublead.push_back(muons.at(1));

    }

    if(passtrigger_double&&matchtrigger_double&&passMETFilter&&hasGoodPrimaryVertex&&is_doubleelectron && is_os && is_isolated){
       double IDSF=1,IDSF_up=1,IDSF_down=1;
       double IDSF_lead = 1, IDSF_sublead = 1;
       double ISOSF_lead = 1, ISOSF_sublead = 1;
       double ISOSF=1,ISOSF_up=1,ISOSF_down=1;
       double triggerSF_double=1., triggerSF_double_up=1.,triggerSF_double_down=1.;

       if(!isData){

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
         triggerSF_double_up=mcdata_correction->GetDoubleMUTriggerEffISR(muons, 1);  //no scale factor for double muon trigger yet
         triggerSF_double_down=mcdata_correction->GetDoubleMUTriggerEffISR(muons, -1);  //no scale factor for double muon trigger yet

          for(int i=0; i<3; i++){
             int sys = 0;
             if(i==1) sys = 1; // up
             if(i==2) sys = -1; // down
             Id1SF.push_back(mcdata_correction->MuonScaleFactor("MUON_POG_TIGHT", muon_lead, sys));
             Id2SF.push_back(mcdata_correction->MuonScaleFactor("MUON_POG_TIGHT", muon_sublead, sys));

             Iso1SF.push_back(mcdata_correction->MuonISOScaleFactor("MUON_POG_TIGHT", muon_lead, sys));
             Iso2SF.push_back(mcdata_correction->MuonISOScaleFactor("MUON_POG_TIGHT", muon_sublead, sys));

             TrigSF.push_back(mcdata_correction->GetDoubleMUTriggerEffISR(muons, sys));
          }

       }

       // for b veto
       std::vector<snu::KJet> jetsPreColl = GetJets("JET_NOLEPTONVETO", 20, 2.4);
       std::vector<snu::KJet> bjets_loose, bjets_medium, bjets_tight;

       for(int i = 0; i < jetsPreColl.size(); i++){
           if( jetsPreColl.at(i).IsBTagged(snu::KJet::CSVv2, snu::KJet::Loose) ) bjets_loose.push_back(jetsPreColl.at(i));
           if( jetsPreColl.at(i).IsBTagged(snu::KJet::CSVv2, snu::KJet::Medium) ) bjets_medium.push_back(jetsPreColl.at(i));
           if( jetsPreColl.at(i).IsBTagged(snu::KJet::CSVv2, snu::KJet::Tight) ) bjets_tight.push_back(jetsPreColl.at(i));
       }

      //////////////Default/////////////////////
      if(l1pt>20&&l2pt>10&&(fabs(l1eta)<2.4)&&(fabs(l2eta)<2.4) && dipt < 100.){
           weightRec=IDSF*ISOSF*triggerSF_double;

           if(!isData){
               // ID systematic
               weightRecIdUp = IDSF_up*ISOSF*triggerSF_double;
               weightRecIdDown = IDSF_down*ISOSF*triggerSF_double;

               // Trigger systematic
               weightRecTriUp = IDSF*ISOSF*triggerSF_double_up;
               weightRecTriDown = IDSF*ISOSF*triggerSF_double_down;

               // Reconstruction systematic
               weightRecRecoUp = IDSF*ISOSF_up*triggerSF_double;
               weightRecRecoDown = IDSF*ISOSF_up*triggerSF_double;

            }
            if(bjets_medium.size() == 0) isBveto=1;
            ispassRec=1;
            nVtx = eventbase->GetEvent().nVertices();
            ptRec.push_back(l1pt);
            ptRec.push_back(l2pt);
            ptRec.push_back(dipt);
            etaRec.push_back(l1eta);
            etaRec.push_back(l2eta);
            etaRec.push_back(dieta);
            phiRec.push_back(l1phi);
            phiRec.push_back(l2phi);
            phiRec.push_back(diphi);
            mRec.push_back(l1mass);
            mRec.push_back(l2mass);
            mRec.push_back(dimass);

      }
    }

  } // two muons


   return;
}// End of execute event loop
  


void ISRmumu_unfolding::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void ISRmumu_unfolding::BeginCycle() throw( LQError ){
  
  Message("In begin Cycle", INFO);
  
  //
  //If you wish to output variables to output file use DeclareVariable
  // clear these variables in ::ClearOutputVectors function
  //DeclareVariable(obj, label, treename );
  //DeclareVariable(obj, label ); //-> will use default treename: LQTree
  //  DeclareVariable(out_electrons, "Signal_Electrons", "LQTree");
  //  DeclareVariable(out_muons, "Signal_Muons");

  DeclareVariable(nVtx,"nVtx","tree");
  DeclareVariable(ptPreFSR,"ptPreFSR","tree");
  DeclareVariable(mPreFSR,"mPreFSR","tree");
  DeclareVariable(phiRec,"phiRec","tree");
  DeclareVariable(etaRec,"etaRec","tree");
  DeclareVariable(ptRec,"ptRec","tree");
  DeclareVariable(mRec,"mRec","tree");
  DeclareVariable(ptGen,"ptGen","tree");
  DeclareVariable(etaGen,"etaGen","tree");
  DeclareVariable(mGen,"mGen","tree");

  DeclareVariable(TrigSF,"TrigSF","tree");
  DeclareVariable(Iso1SF,"Iso1SF","tree");
  DeclareVariable(Iso2SF,"Iso2SF","tree");
  DeclareVariable(Id1SF,"Id1SF","tree");
  DeclareVariable(Id2SF,"Id2SF","tree");


  DeclareVariable(ispassRec,"ispassRec","tree");
  DeclareVariable(isBveto,"isBveto","tree");
  DeclareVariable(isfiducialGen,"isfiducialGen","tree");
  DeclareVariable(isfiducialPreFSR,"isfiducialPreFSR","tree");

  // systematic weight
  DeclareVariable(weightGenScale,"weightGenScale","tree");
  DeclareVariable(weightGenPdf,"weightGenPdf","tree");
  DeclareVariable(weightGen,"weightGen","tree");
  DeclareVariable(weightRec,"weightRec","tree");
  DeclareVariable(weightRecIdUp,"weightRecIdUp","tree");
  DeclareVariable(weightRecIdDown,"weightRecIdDown","tree");
  DeclareVariable(weightRecTriUp,"weightRecTriUp","tree");
  DeclareVariable(weightRecTriDown,"weightRecTriDown","tree");
  DeclareVariable(weightRecRecoUp,"weightRecRecoUp","tree");
  DeclareVariable(weightRecRecoDown,"weightRecRecoDown","tree");
  DeclareVariable(weightGenPileUp,"weightGenPileUp","tree");
  DeclareVariable(weightGenPileDown,"weightGenPileDown","tree");

  DeclareVariable(DYtautau,"DYtautau","tree");

  
  return;
  
}

ISRmumu_unfolding::~ISRmumu_unfolding() {
  
  Message("In ISRmumu_unfolding Destructor" , INFO);
  
}


void ISRmumu_unfolding::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void ISRmumu_unfolding::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  maphist3D.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this ISRmumu_unfoldingCore::MakeHistograms() to make new hists for your analysis
   **/

}


void ISRmumu_unfolding::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  //
  isfiducialGen = 0;
  isfiducialPreFSR = 0;
  ispassRec = 0;
  isBveto = 0;
  DYtautau = 0;
  nVtx = 0;

  weightGen = 1.;
  weightRec = 1.;
  weightRecIdUp = 1.;
  weightRecIdDown = 1.;
  weightRecTriUp = 1.;
  weightRecTriDown = 1.;
  weightRecRecoUp = 1.;
  weightRecRecoDown = 1.;
  weightGenPileUp = 1.;
  weightGenPileDown = 1.;
  weightGenScale.clear();
  weightGenPdf.clear();

  TrigSF.clear();
  Iso1SF.clear();
  Iso2SF.clear();
  Id1SF.clear();
  Id2SF.clear();

  phiRec.clear();
  etaRec.clear();
  ptRec.clear();
  mRec.clear();

  etaGen.clear();
  ptGen.clear();
  mGen.clear();

  ptPreFSR.clear();
  mPreFSR.clear();
}



