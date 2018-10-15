// $Id: ISRee_FSRunfolding.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQISRee_FSRunfolding Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "ISRee_FSRunfolding.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (ISRee_FSRunfolding);


 /**
  *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
  *
  */
ISRee_FSRunfolding::ISRee_FSRunfolding() :  AnalyzerCore(), out_electrons(0)  {
  
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("ISRee_FSRunfolding");
  
  Message("In ISRee_FSRunfolding constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();
  //  MakeCleverHistograms(sighist_mm,"DiElectron");


}


void ISRee_FSRunfolding::InitialiseAnalysis() throw( LQError ) {
  
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


void ISRee_FSRunfolding::ExecuteEvents()throw( LQError ){

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
        genl1_postFSR+=truth;
      }
      else if(truth.PdgId()==-11&&truth.StatusFlag(snu::KTruth::fromhardprocess)){
	genl2+=truth;
	genl2_postFSR+=truth;
      }
      else if(truth.PdgId()==22){   //collect photons from DY muons by FSR
	int imother=truth.IndexMother();
	if(imother>=truthsize) continue;
	snu::KTruth mother=truthcol.at(truth.IndexMother());

        if(mother.PdgId() == 11 && mother.StatusFlag(snu::KTruth::fromhardprocess)){
          genl1+=truth;
        }
        if(mother.PdgId() == -11 && mother.StatusFlag(snu::KTruth::fromhardprocess)){
          genl2+=truth;
        }
      }
    }
    
    if(genl1!=TLorentzVector(0,0,0,0)&&genl2!=TLorentzVector(0,0,0,0)){

      // full phase space
      weight_ = weight*weigtbytrig;
      ptPreFSR.push_back((genl1+genl2).Pt());
      mPreFSR.push_back((genl1+genl2).M());

      if(fabs(genl1.Eta()) < 2.4 && fabs(genl2.Eta()) < 2.4){
         if((genl1.Pt() > 25 && genl2.Pt() > 15)||(genl2.Pt() > 25 && genl1.Pt() > 15)){
            issignal = 1;

            //weight_ = weight*weigtbytrig;
            //ptPreFSR.push_back((genl1+genl2).Pt());
            //mPreFSR.push_back((genl1+genl2).M());
         }
      }

      if(fabs(genl1_postFSR.Eta()) < 2.4 && fabs(genl2_postFSR.Eta()) < 2.4){
         if((genl1_postFSR.Pt() > 25 && genl2_postFSR.Pt() > 15)||(genl2_postFSR.Pt() > 25 && genl1_postFSR.Pt() > 15)){
            istriggered = 1;

            weight_ = weight*weigtbytrig;
            ptPostFSR.push_back((genl1_postFSR+genl2_postFSR).Pt());
            mPostFSR.push_back((genl1_postFSR+genl2_postFSR).M());
         }
      }
    }
  }
  return;
}// End of execute event loop
  


void ISRee_FSRunfolding::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void ISRee_FSRunfolding::BeginCycle() throw( LQError ){
  
  Message("In begin Cycle", INFO);
  
  //
  //If you wish to output variables to output file use DeclareVariable
  // clear these variables in ::ClearOutputVectors function
  //DeclareVariable(obj, label, treename );
  //DeclareVariable(obj, label ); //-> will use default treename: LQTree
  //  DeclareVariable(out_electrons, "Signal_Electrons", "LQTree");
  //  DeclareVariable(out_muons, "Signal_Muons");

  DeclareVariable(ptPreFSR,"ptPreFSR","tree"); 
  DeclareVariable(mPreFSR,"mPreFSR","tree"); 
  DeclareVariable(istriggered,"istriggered","tree"); 
  DeclareVariable(weight_,"weight","tree"); 
  DeclareVariable(ptPostFSR,"ptPostFSR","tree"); 
  DeclareVariable(mPostFSR,"mPostFSR","tree"); 
  DeclareVariable(issignal,"issignal","tree"); 

  return;
  
}

ISRee_FSRunfolding::~ISRee_FSRunfolding() {
  
  Message("In ISRee_FSRunfolding Destructor" , INFO);
  
}


void ISRee_FSRunfolding::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void ISRee_FSRunfolding::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this ISRee_FSRunfoldingCore::MakeHistograms() to make new hists for your analysis
   **/

}


void ISRee_FSRunfolding::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //

  issignal = 0;
  istriggered = 0;
  weight_ = 1.;
  weightTotal = 1.;

  ptPreFSR.clear();
  mPreFSR.clear();

  ptPostFSR.clear();
  mPostFSR.clear();
}



