#ifndef AnalyzerCore_H
#define AnalyzerCore_H

//forward declarations                                                                                                                                            
class Reweight;
class EventBase;
class MuonPlots;
class ElectronPlots;
class JetPlots;
class SignalPlots;
class EventBase;

#include "LQCycleBase.h"
#include "HNCommonLeptonFakes/HNCommonLeptonFakes/HNCommonLeptonFakes.h"
#include "rochcor2012/rochcor2012/rochcor2012jan22.h"

class AnalyzerCore : public LQCycleBase {
  
 public:
  
  // Default constructor
  AnalyzerCore();

  //destructor
  virtual ~AnalyzerCore();

  // SetUpEvent CORE function: accesses event in ntuple
  virtual void SetUpEvent(Long64_t entry, float ev_weight)throw( LQError );
  virtual void EndEvent()throw( LQError );
  virtual void WriteHistograms()throw( LQError );


  TDirectory*   getTemporaryDirectory(void) const;

  double ElectronScaleFactor( double eta, double pt);
  double MuonScaleFactor(double eta, double pt);
  float  JetResCorr(snu::KJet jet, std::vector<snu::KGenJet> genjets);
  float SumPt( std::vector<snu::KJet> particles);
  bool isPrompt(long pdgid);
  bool IsTight(snu::KElectron el, std::vector<snu::KJet> jets, double jetrho , double dxy, double biso, double eiso, bool usedr3, bool usetrkiso, bool usetight);
  bool IsTight(snu::KElectron electron, std::vector<snu::KJet> jets, double rho);
  bool IsTight(snu::KMuon muon);
  std::vector<snu::KElectron> GetTruePrompt(vector<snu::KElectron> electrons,  bool keep_chargeflip, bool keepfake);
  int NBJet(std::vector<snu::KJet> jets);
  bool Zcandidate(vector<snu::KElectron> electrons, float interval, bool require_os=true);
  bool SameCharge(std::vector<snu::KElectron> electrons);
  float CFRate(snu::KElectron el);
  void CorrectMuonMomentum(vector<snu::KMuon>& k_muons);
  std::vector<snu::KElectron>  ShiftElectronEnergy(std::vector<snu::KElectron> el, bool applyshift);
  float Get_DataDrivenWeight_E(vector<snu::KElectron> k_electrons, vector<snu::KJet> jets ,std::vector<snu::KJet> alljets , int njets, int nbjets, double rho, double dxy, double biso, double eiso, bool usedr3, bool usetrkiso, bool    usetight,TString cut, bool applypucorr);

  float Get_DataDrivenWeight_EE(vector<snu::KElectron> k_electrons, std::vector<snu::KJet> alljets ,int njets,  double rho);
  float Get_DataDrivenWeight_EE(vector<snu::KElectron> k_electrons, std::vector<snu::KJet> alljets ,int njets, double rho, double dxy, double biso, double eiso, bool usedr3, bool usetrkiso, bool usetight,TString cut ,int nbjet, float ht, bool user, bool useht, bool usenjet);
  float Get_DataDrivenWeight_EE(vector<snu::KElectron> k_electrons,std::vector<snu::KJet> alljets , int njets, double rho, double dxy, double biso, double eiso, bool usedr3, bool usetrkiso, bool usetight,TString cut ,int nbjet, float ht, bool close1, bool close2);
  float Get_DataDrivenWeight_EE(vector<snu::KElectron> k_electrons, std::vector<snu::KJet> alljets ,int njets, double rho, double dxy, double biso, double eiso, bool usedr3, bool usetrkiso, bool usetight,TString cut );
  float Get_DataDrivenWeight_r1_EE(vector<snu::KElectron> k_electrons,std::vector<snu::KJet> alljets , int njets, double rho, double dxy, double biso, double eiso, bool usedr3, bool usetrkiso, bool  usetight,TString cut, int eventtype, bool setr1);
  float Get_DataDrivenWeight_MM(vector<snu::KMuon> k_muons);
  float Get_DataDrivenWeight_EM(vector<snu::KMuon> k_muons, vector<snu::KElectron> k_electrons, std::vector<snu::KJet> alljets ,int njets, double rho);
  
  double MuonDYMassCorrection(std::vector<snu::KMuon> mu, double w);

  
  vector<TLorentzVector> MakeTLorentz( vector<snu::KElectron> el);
  vector<TLorentzVector> MakeTLorentz( vector<snu::KMuon> mu);
  vector<TLorentzVector> MakeTLorentz( vector<snu::KJet> jet);
  // enum for plotting functions/classes
  enum histtype {muhist, elhist, jethist, sighist};
  
  
  //
  // Useful message function 
  //
  void Message(TString message, LQMsgType type=INFO);


  //
  //  Specify which triggers will be avaiable in KTrigger
  //
  void AddTriggerToList(TString triggername);
  
  /// Pileup Reweighting class
  static const Bool_t MC_pu = true;
  Reweight *reweightPU;

  //// Event base pointer. Used to get all objects for analysis
  EventBase* eventbase;
  
  UInt_t numberVertices;
  Bool_t *goodVerticiesB;

  TDirectory *Dir;
  map<TString, TH1*> maphist;
  map<TString, TH2*> maphist2D;
  TH2F* FRHist;
  TH2F* MuonSF;
  HNCommonLeptonFakes* m_fakeobj;
  rochcor2012 *rmcor;
  
  /// Event weights
  Double_t MCweight, weight;

  // used to get trigger prescale
  Int_t prescale;
  
  std::vector<TString> triggerlist;

  //// Making cleaver hist maps
  map<TString, SignalPlots*> mapCLhistSig;
  map<TString, ElectronPlots*> mapCLhistEl;
  map<TString, MuonPlots*> mapCLhistMu;
  map<TString, JetPlots*> mapCLhistJet;
  
  //
  // Function that closes rootfile
  //
  void CloseFiles();
  

  //
  // Make Histograms and fill maphist
  //
  void MakeHistograms();
  void MakeHistograms(TString hname, int nbins, float xmin, float xmax);
  void MakeHistograms(TString hname, int nbins, float xbins[]);
  void MakeHistograms2D(TString hname, int nbinsx, float xbins[], int nbinsy, float ybins[]);
  void MakeHistograms2D(TString hname, int nbinsx, float xmin, float xmax, int nbinsy, float ymin, float ymax);
    //
    // Makes temporary dir
    //
    TDirectory* GetTemporaryDirectory(void) const;                                                                                                                                 
  //
  // Checks if a file exists
  //
  void CheckFile(TFile* file) throw( LQError );
  
  //// Plotting 
  TH1* GetHist(TString hname);
  TH2* GetHist2D(TString hname);

  /// Fills hist in maphist
  void FillHist(TString histname, float value, float w );
  void FillHist(TString histname, float value, float w , float xmin, float xmax, int nbins=0);
  void FillHist(TString histname, float value, float w , float xmin[], int nbins=0);
  void FillHist(TString histname, float value1, float value2, float w , float x[], int nbinsx, float y[], int nbinsy);
  void FillHist(TString histname, float value1,  float value2, float w , float xmin, float xmax, int nbinsx,  float ymin, float ymax, int nbinsy);
  /// Fills clever hists
  void FillCLHist(histtype type, TString hist, snu::KEvent ev,vector<snu::KMuon> muons, vector<snu::KElectron> electrons, vector<snu::KJet> jets,double weight);
  void FillCLHist(histtype type, TString hist, snu::KEvent ev,vector<snu::KMuon> muons, vector<snu::KJet> jets,double weight);
  void FillCLHist(histtype type, TString hist, snu::KEvent ev, vector<snu::KElectron> electrons, vector<snu::KJet> jets,double weight);
  void FillCLHist(histtype type, TString hist, vector<snu::KMuon> muons , double weight);
  void FillCLHist(histtype type, TString hist, vector<snu::KElectron> electrons , double rho, double weight);
  void FillCLHist(histtype type, TString hist, vector<snu::KJet> jets , double weight);

  // Makes clever histograms
  void MakeCleverHistograms(histtype type, TString clhistname );

  /// File related                                                                                                                                                
  void OpenPutputFile();
  void WriteHists();
  void WriteCLHists();

  //// Event related                                                                                                                                              
  bool PassTrigger(std::vector<TString> list, int& prescale);
  bool PassBasicEventCuts();

};
#endif
