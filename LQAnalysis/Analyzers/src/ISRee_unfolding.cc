// $Id: ISRee_unfolding.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQISRee_unfolding Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "ISRee_unfolding.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (ISRee_unfolding);


 /**
  *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
  *
  */
ISRee_unfolding::ISRee_unfolding() :  AnalyzerCore(), out_electrons(0)  {
  
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("ISRee_unfolding");
  
  Message("In ISRee_unfolding constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();
  //  MakeCleverHistograms(sighist_mm,"DiElectron");


}


void ISRee_unfolding::InitialiseAnalysis() throw( LQError ) {
  
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

void ISRee_unfolding::ExecuteEvents()throw( LQError ){

 double totWeight = 1.;
 // double electron trigger
 TString dielectron_trig="HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v";//

 vector<TString> trignames;
 trignames.push_back(dielectron_trig);

 double weigtbytrig = WeightByTrigger(dielectron_trig,TargetLumi);

 /// Apply the gen weight 
 if(!isData) weight*=MCweight;
 weightGen=weight*weigtbytrig;
 
 totWeight = weight*weigtbytrig;

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
     if(truth.PdgId()==11&&truth.StatusFlag(snu::KTruth::fromhardprocess)){
       genl1+=truth;
     }
     else if(truth.PdgId()==-11&&truth.StatusFlag(snu::KTruth::fromhardprocess)){
       genl2+=truth;
     }
     else if(truth.PdgId()==22){   //collect photons from DY muons by FSR
       int imother=truth.IndexMother();
       if(imother>=truthsize) continue;
       snu::KTruth mother=truthcol.at(truth.IndexMother());
       if(mother.PdgId()==11&&mother.StatusFlag(snu::KTruth::fromhardprocess)){
         genl1pre+=truth; //collect photons first, then later add to post fsr lepton i.e., gen1 or gen2
       }else if(mother.PdgId()==-11&&mother.StatusFlag(snu::KTruth::fromhardprocess)){
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
      mGen.push_back(0.000511);
      mGen.push_back(0.000511);
      mGen.push_back((genl1+genl2).M());

      if((genl1pre.Pt()>25||genl2pre.Pt()>25)&&(genl1pre.Pt()>15&&genl2pre.Pt()>15)&&fabs(genl1pre.Eta())<2.5&&fabs(genl2pre.Eta())<2.5){
	isfiducialPreFSR = 1;
      }
      if(genl1pt>25&&genl2pt>15&&fabs(genl1eta)<2.5&&fabs(genl2eta)<2.5){
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

  std::vector<snu::KElectron> electrons =  GetElectrons(true,true, "ELECTRON_POG_MEDIUM"); // Cut Based POG Medium WP 
  bool trig_pass= PassTriggerOR(trignames);

  double idsf = mcdata_correction->ElectronScaleFactor("ELECTRONISR_POG_MEDIUM", electrons, 0);
  double idsf_up = mcdata_correction->ElectronScaleFactor("ELECTRONISR_POG_MEDIUM", electrons, 1);
  double idsf_down = mcdata_correction->ElectronScaleFactor("ELECTRONISR_POG_MEDIUM", electrons, -1);

  double recosf = mcdata_correction->ElectronRecoScaleFactor(electrons, 0);
  double recosf_up = mcdata_correction->ElectronRecoScaleFactor(electrons, 1);
  double recosf_down = mcdata_correction->ElectronRecoScaleFactor(electrons, -1);
 
  double trigsf = 1.; 
  trigsf = mcdata_correction->GetDoubleEGTriggerEffISR(electrons,0);

  std::vector<snu::KPhoton> photons = eventbase->GetPhotons();
  std::vector<snu::KJet> jets = eventbase->GetJets();

  double l1prefire = 1.;
  if(!isData) l1prefire = mcdata_correction->GetL1ECALPrefiringWeight(photons, jets, 0);

  bool is_doubleelectron = (electrons.size() == 2);

  if(PassMETFilter()){     /// Initial event cuts : 
     if(eventbase->GetEvent().HasGoodPrimaryVertex()){ //// Make cut on event wrt vertex                                                                               
 
        // medium working point pog id without ip cuts
        if(is_doubleelectron){
        
          bool is_os = (electrons.at(0).Charge() == (-electrons.at(1).Charge()));

          double trig_sf = 1.;
          double trig_sf_up = 1.;
          double trig_sf_down = 1.;
          if(!isData){
            trig_sf = mcdata_correction->GetDoubleEGTriggerEffISR(electrons,0);
            trig_sf_up = mcdata_correction->GetDoubleEGTriggerEffISR(electrons,1);
            trig_sf_down = mcdata_correction->GetDoubleEGTriggerEffISR(electrons,-1);
          }

          double l1pt = electrons[0].Pt();
          double l2pt = electrons[1].Pt();
          double l1eta = electrons[0].Eta();
          double l2eta = electrons[1].Eta();
          double l1phi = electrons[0].Phi();
          double l2phi = electrons[1].Phi();
          double l1mass = electrons[0].M();
          double l2mass = electrons[1].M();
          double dipt = (electrons[0]+electrons[1]).Pt();
          double dieta = (electrons[0]+electrons[1]).Eta();
          double diphi = (electrons[0]+electrons[1]).Phi();
          double dimass = (electrons[0]+electrons[1]).M();
          //bool TriggerMatch1 = (electrons[0].TriggerMatched(dielectron_trig) && electrons[1].TriggerMatched(dielectron_trig));
          bool TriggerMatch1 = true; 

          std::vector<snu::KElectron> electron1; 
          std::vector<snu::KElectron> electron2; 
          electron1.push_back(electrons[0]);
          electron2.push_back(electrons[1]);

          for(int i=0; i<3; i++){
             int sys = 0;
             if(i==1) sys = 1; // up
             if(i==2) sys = -1; // down
             Id1SF.push_back(mcdata_correction->ElectronScaleFactor("ELECTRONISR_POG_MEDIUM", electron1, sys));
             Id2SF.push_back(mcdata_correction->ElectronScaleFactor("ELECTRONISR_POG_MEDIUM", electron2, sys));
             
             Reco1SF.push_back(mcdata_correction->ElectronRecoScaleFactor(electron1, sys)); 
             Reco2SF.push_back(mcdata_correction->ElectronRecoScaleFactor(electron2, sys)); 

             TrigSF.push_back(mcdata_correction->GetDoubleEGTriggerEffISR(electrons,sys));
          }

          bool mcfromtau = (electrons[0].MCFromTau()||electrons[1].MCFromTau()); // is there case only one of the lepton decaying from tau in DY sample?

          // for b veto
          std::vector<snu::KJet> jetsPreColl = GetJets("JET_NOLEPTONVETO", 20, 2.4);
          std::vector<snu::KJet> bjets_loose, bjets_medium, bjets_tight;

          for(int i = 0; i < jetsPreColl.size(); i++){
              if( jetsPreColl.at(i).IsBTagged(snu::KJet::CSVv2, snu::KJet::Loose) ) bjets_loose.push_back(jetsPreColl.at(i));
              if( jetsPreColl.at(i).IsBTagged(snu::KJet::CSVv2, snu::KJet::Medium) ) bjets_medium.push_back(jetsPreColl.at(i));
              if( jetsPreColl.at(i).IsBTagged(snu::KJet::CSVv2, snu::KJet::Tight) ) bjets_tight.push_back(jetsPreColl.at(i));
          }

          if(l1pt < l2pt) std::cout << "l1pt: " << l1pt << " l2pt: " << l2pt << std::endl;
          if(trig_pass && is_os && TriggerMatch1 && l1pt > 25. && l2pt > 15. && fabs(l1eta) < 2.5 && fabs(l2eta) < 2.5 && dipt < 100.){


            TString prefix="";
            if(mcfromtau && k_sample_name.Contains("DY") ){ 
              DYtautau = 1;
              prefix="tau_";
            }

            weightRec = idsf*recosf*trig_sf;
 
            if(!isData){ 
               // ID systematic
               weightRecIdUp = idsf_up*recosf*trig_sf;
               weightRecIdDown = idsf_down*recosf*trig_sf;

               // Trigger systematic
               weightRecTriUp = idsf*recosf*trig_sf_up;
               weightRecTriDown = idsf*recosf*trig_sf_down;

               // Reconstruction systematic
               weightRecRecoUp = idsf*recosf_up*trig_sf;
               weightRecRecoDown = idsf*recosf_down*trig_sf;

               l1PreFire.push_back(mcdata_correction->GetL1ECALPrefiringWeight(photons, jets, 0));
               l1PreFire.push_back(mcdata_correction->GetL1ECALPrefiringWeight(photons, jets, 1));
               l1PreFire.push_back(mcdata_correction->GetL1ECALPrefiringWeight(photons, jets, 2));

               //std::cout << "photon size: " << photons.size() << " jet size: " << jets.size() << std::endl;
               //std::cout << "l1PreFire up: " << mcdata_correction->GetL1ECALPrefiringWeight(photons, jets, 1) << " norm: " << mcdata_correction->GetL1ECALPrefiringWeight(photons, jets, 0) << " down: " << mcdata_correction->GetL1ECALPrefiringWeight(photons, jets, 2) << std::endl;
            }

            FillHist(prefix+"dielectronmassRecoSFTrgSFIDSF",dimass,weightGen*idsf*recosf*trig_sf, 40., 350, 620);
            if(bjets_medium.size() == 0) isBveto=1;
            lepSmearFactor.push_back(electrons[0].SmearFactor());
            lepSmearFactor.push_back(electrons[1].SmearFactor());
            ispassRec=1;
            nVtx = eventbase->GetEvent().nVertices();
            qLep.push_back(electrons[0].Charge());
            qLep.push_back(electrons[1].Charge());
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

          } // event selection for opposite sign electrons 

        }// exactly two pog wp electrons
     } // Good PV filter
  } //MET filter

  TString prefix = "";
  if(DYtautau == 1) 
    prefix="tau_";

  // cutflow
  if(eventbase->GetEvent().HasGoodPrimaryVertex() && PassMETFilter() && PassTriggerOR(trignames)){
    FillHist(prefix+"cutflow",0.5,totWeight,0.,20.,20); 
    FillHist(prefix+"cutflow",1.5,totWeight*pileup_reweight,0.,20.,20); 

    if(electrons.size() > 1){ 
      FillHist(prefix+"cutflow",2.5,totWeight*pileup_reweight,0.,20.,20);
      if(electrons.at(0).Pt()>25.&&electrons.at(1).Pt()>15 && fabs(electrons.at(0).Eta())<2.5&&fabs(electrons.at(1).Eta())<2.5){
         FillHist(prefix+"cutflow",3.5,totWeight*pileup_reweight,0.,20.,20);
         if(electrons.size() == 2){
            FillHist(prefix+"cutflow",4.5,totWeight*pileup_reweight,0.,20.,20);
            FillHist(prefix+"cutflow",5.5,totWeight*pileup_reweight*recosf,0.,20.,20);
            FillHist(prefix+"cutflow",6.5,totWeight*pileup_reweight*recosf*idsf,0.,20.,20);
            FillHist(prefix+"cutflow",7.5,totWeight*pileup_reweight*recosf*idsf,0.,20.,20);
            FillHist(prefix+"cutflow",8.5,totWeight*pileup_reweight*recosf*idsf*trigsf,0.,20.,20);
            if((electrons.at(0)+electrons.at(1)).M()>60){
               FillHist(prefix+"cutflow",9.5,totWeight*pileup_reweight*recosf*idsf*trigsf,0.,20.,20);
               FillHist(prefix+"cutflow",10.5,totWeight*pileup_reweight*recosf*idsf*trigsf*l1prefire,0.,20.,20);
               if((electrons.at(0).Charge() == (-electrons.at(1).Charge()))){
                  FillHist(prefix+"cutflow",11.5,totWeight*pileup_reweight*recosf*idsf*trigsf*l1prefire,0.,20.,20);
                  if((electrons.at(0)+electrons.at(1)).M()<120){
                     FillHist(prefix+"cutflow",12.5,totWeight*pileup_reweight*recosf*idsf*trigsf*l1prefire,0.,20.,20);
                  }
               } // os
            }// mass > 60 GeV
         } // exactly two electrons
      }// kinematic cut 
    }// at least 2 electron
  }// basic fillter & trigger 
   return;
}// End of execute event loop
  


void ISRee_unfolding::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void ISRee_unfolding::BeginCycle() throw( LQError ){
  
  Message("In begin Cycle", INFO);
  
  //
  //If you wish to output variables to output file use DeclareVariable
  // clear these variables in ::ClearOutputVectors function
  //DeclareVariable(obj, label, treename );
  //DeclareVariable(obj, label ); //-> will use default treename: LQTree
  //  DeclareVariable(out_electrons, "Signal_Electrons", "LQTree");
  //  DeclareVariable(out_muons, "Signal_Muons");

//  DeclareVariable(nVtx,"nVtx","tree"); 
//  DeclareVariable(ptPreFSR,"ptPreFSR","tree"); 
//  DeclareVariable(mPreFSR,"mPreFSR","tree");
//  DeclareVariable(lepSmearFactor,"lepSmearFactor","tree"); 
//  DeclareVariable(qLep,"qLep","tree"); 
//  DeclareVariable(etaRec,"etaRec","tree"); 
//  DeclareVariable(phiRec,"phiRec","tree"); 
//  DeclareVariable(ptRec,"ptRec","tree"); 
//  DeclareVariable(mRec,"mRec","tree");
//  DeclareVariable(ptGen,"ptGen","tree"); 
//  DeclareVariable(etaGen,"etaGen","tree"); 
//  DeclareVariable(mGen,"mGen","tree"); 
//
//  DeclareVariable(TrigSF,"TrigSF","tree"); 
//  DeclareVariable(Reco1SF,"Reco1SF","tree"); 
//  DeclareVariable(Reco2SF,"Reco2SF","tree"); 
//  DeclareVariable(Id1SF,"Id1SF","tree"); 
//  DeclareVariable(Id2SF,"Id2SF","tree"); 
//  DeclareVariable(l1PreFire,"l1PreFire","tree"); 
// 
//  DeclareVariable(ispassRec,"ispassRec","tree"); 
//  DeclareVariable(isBveto,"isBveto","tree");
//  DeclareVariable(isfiducialGen,"isfiducialGen","tree"); 
//  DeclareVariable(isfiducialPreFSR,"isfiducialPreFSR","tree"); 
//
//  // systematic weight
//  DeclareVariable(weightGenScale,"weightGenScale","tree"); 
//  DeclareVariable(weightGenPdf,"weightGenPdf","tree");
//  DeclareVariable(weightGen,"weightGen","tree"); 
//  DeclareVariable(weightRec,"weightRec","tree"); 
//  DeclareVariable(weightRecIdUp,"weightRecIdUp","tree"); 
//  DeclareVariable(weightRecIdDown,"weightRecIdDown","tree"); 
//  DeclareVariable(weightRecTriUp,"weightRecTriUp","tree"); 
//  DeclareVariable(weightRecTriDown,"weightRecTriDown","tree"); 
//  DeclareVariable(weightRecRecoUp,"weightRecRecoUp","tree"); 
//  DeclareVariable(weightRecRecoDown,"weightRecRecoDown","tree"); 
//  DeclareVariable(weightGenPileUp,"weightGenPileUp","tree"); 
//  DeclareVariable(weightGenPileDown,"weightGenPileDown","tree"); 
//
//  DeclareVariable(DYtautau,"DYtautau","tree"); 

  return;
}

ISRee_unfolding::~ISRee_unfolding() {
  
  Message("In ISRee_unfolding Destructor" , INFO);
  
}


void ISRee_unfolding::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void ISRee_unfolding::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this ISRee_unfoldingCore::MakeHistograms() to make new hists for your analysis
   **/

}


void ISRee_unfolding::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
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

  qLep.clear();
  lepSmearFactor.clear();
  l1PreFire.clear();
  TrigSF.clear();
  Reco1SF.clear();
  Reco2SF.clear();
  Id1SF.clear();
  Id2SF.clear();

  weightGenScale.clear();
  weightGenPdf.clear();

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

