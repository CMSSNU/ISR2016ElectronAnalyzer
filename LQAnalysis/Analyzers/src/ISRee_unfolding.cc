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

  // double electron trigger
  TString dielectron_trig="HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v";//

  vector<TString> trignames;
  trignames.push_back(dielectron_trig);

  double weigtbytrig = WeightByTrigger(dielectron_trig,TargetLumi);

  // Mass boundaries  
  double massBinBoundaries[] = {40., 50., 60., 80., 100., 200., 350};
  const int nBoundaries = sizeof(massBinBoundaries)/sizeof(double); 

  vector<TString> massbins;
  massbins.push_back("m40to50");
  massbins.push_back("m50to60");
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
    TLorentzVector gendy_postFSR,genl1_postFSR,genl2_postFSR;   
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
      }
      else if(truth.PdgId()==-11&&truth.StatusFlag(snu::KTruth::fromhardprocess)){
	genl2+=truth;
	gendy+=truth;
	genl2_postFSR+=truth;
	gendy_postFSR+=truth;
      }
      else if(truth.PdgId()==22){   //collect photons from DY muons by FSR
	int imother=truth.IndexMother();
	if(imother>=truthsize) continue;
	snu::KTruth mother=truthcol.at(truth.IndexMother());
	if(abs(mother.PdgId())==11&&mother.StatusFlag(snu::KTruth::fromhardprocess)){
	  gendy+=truth;
	}
      }
    }
    
    if(genl1!=TLorentzVector(0,0,0,0)&&genl2!=TLorentzVector(0,0,0,0)){

      if(fabs(genl1_postFSR.Eta()) < 2.4 && fabs(genl2_postFSR.Eta()) < 2.4){
         if((genl1_postFSR.Pt() > 25 && genl2_postFSR.Pt() > 15)||(genl2_postFSR.Pt() > 25 && genl1_postFSR.Pt() > 15)){
            issignal = 1;

            weight_ = weight*weigtbytrig;
            ptGen.push_back(genl1.Pt()); // FIXME
            ptGen.push_back(genl2.Pt());
            ptGen.push_back((gendy_postFSR).Pt());

            mGen.push_back(0.000511);
            mGen.push_back(0.000511);
            mGen.push_back((gendy_postFSR).M());
         }
      }
    }
  }
   
  float pileup_reweight=(1.0);
  if (!k_isdata) {
    pileup_reweight=mcdata_correction->CatPileupWeight(eventbase->GetEvent(),0);
  }
  
  std::vector<snu::KElectron> electrons =  GetElectrons(true,true, "ELECTRON_POG_MEDIUM"); // Cut Based POG Medium WP 

  bool trig_pass= PassTriggerOR(trignames);

  double idsf = mcdata_correction->ElectronScaleFactor("ELECTRON_POG_MEDIUM", electrons, 0);
  double recosf = mcdata_correction->ElectronRecoScaleFactor(electrons, 0);

  bool is_doubleelectron = (electrons.size() == 2);

  double dimass_min = 0., dimass_max = 500.;
  int ndimass = 1000;

  if(PassMETFilter()){     /// Initial event cuts : 
     if(eventbase->GetEvent().HasGoodPrimaryVertex()){ //// Make cut on event wrt vertex                                                                               
 
        // medium working point pog id without ip cuts
        if(is_doubleelectron){
        
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

          if(trig_pass && is_os && TriggerMatch1 && ptlep1 > 25. && ptlep2 > 15. && fabs(etalep1) < 2.4 && fabs(etalep2) < 2.4 && dipt < 100.){

            istriggered = 1;
            if(mcfromtau && k_sample_name.Contains("DY") ) DYtautau = 1;
            if(DYtautau){
              std::cout << "name: "<< k_sample_name << std::endl;
              std::cout << "from tau?: " << mcfromtau << std::endl;
            }

               weightTotal = weight*weigtbytrig*idsf*recosf*pileup_reweight;

               ptRec.push_back(ptlep1);
               ptRec.push_back(ptlep2);
               ptRec.push_back(dipt);

               etaRec.push_back(etalep1);
               etaRec.push_back(etalep2);
               etaRec.push_back((electrons[0]+electrons[1]).Eta());

               mRec.push_back(electrons[0].M());
               mRec.push_back(electrons[1].M());
               mRec.push_back((electrons[0]+electrons[1]).M());

               for(int i = 0; i < nBoundaries-1; i++){
                  if(dimass > massBinBoundaries[i] && dimass < massBinBoundaries[i+1]){
                    // dilepton pt for each mass bin
                    FillHist(prefix+"dielectronmass_"+massbins[i],dimass,weight*idsf*pileup_reweight*recosf*weigtbytrig, dimass_min, dimass_max, ndimass);
                    FillHist(prefix+"dielectronmass_"+massbins[i]+"_weightTotal",dimass,weightTotal, dimass_min, dimass_max, ndimass);

                  }
               }

          } // event selection for opposite sign electrons 

        }// exactly two pog wp electrons
     } // Good PV filter
  } //MET filter
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

  DeclareVariable(etaRec,"etarec","tree"); 
  DeclareVariable(ptRec,"ptrec","tree"); 
  DeclareVariable(mRec,"mrec","tree"); 
  DeclareVariable(istriggered,"istriggered","tree"); 
  DeclareVariable(weight_,"weight","tree"); 
  DeclareVariable(weightTotal,"weightTotal","tree"); 
  DeclareVariable(ptGen,"ptgen","tree"); 
  DeclareVariable(etaGen,"etagen","tree"); 
  DeclareVariable(mGen,"mgen","tree"); 
  DeclareVariable(issignal,"issignal","tree"); 
  DeclareVariable(DYtautau,"DYtautau","tree"); 

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

  issignal = 0;
  istriggered = 0;
  DYtautau = 0;
  weight_ = 1.;
  weightTotal = 1.;

  etaRec.clear();
  ptRec.clear();
  mRec.clear();

  etaGen.clear();
  ptGen.clear();
  mGen.clear();
}



