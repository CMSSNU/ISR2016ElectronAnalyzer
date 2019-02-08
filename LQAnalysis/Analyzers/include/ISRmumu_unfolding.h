#ifndef ISRmumu_unfolding_h
#define ISRmumu_unfolding_h

#include "AnalyzerCore.h"
class ISRmumu_unfolding : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  ISRmumu_unfolding();
  ~ISRmumu_unfolding();

  /// Functions from core
  virtual void BeginCycle() throw( LQError );
  virtual void BeginEvent()throw( LQError );
  virtual void ExecuteEvents()throw( LQError );
  virtual void EndCycle()throw( LQError );
  virtual void ClearOutputVectors()throw( LQError );
  
  void InitialiseAnalysis() throw( LQError );
  void MakeHistograms();
 private:
  
  //
  // The output variables 
  //
  /// Vectors for output objetcs
  std::vector<snu::KMuon> out_muons;
  std::vector<snu::KElectron> out_electrons;

  // https://github.com/CMSSNU/ISR2016MuonAnalyzer/blob/master/LQAnalysis/Analyzers/include/ISR2016MuonAnalyzer.h#L30
  std::vector<Double_t> phiRec, etaRec,ptRec,mRec;
  std::vector<Double_t> etaGen,ptGen,mGen;
  std::vector<Double_t> ptPreFSR,mPreFSR;
  std::vector<Double_t> TrigSF, Iso1SF, Iso2SF, Id1SF, Id2SF;
  Double_t weightGen, weightRec;
  Double_t weightRecIdUp, weightRecIdDown, weightRecTriUp, weightRecTriDown, weightRecRecoUp, weightRecRecoDown;
  Double_t weightGenPileUp, weightGenPileDown;
  std::vector<Double_t> weightGenScale, weightGenPdf;
  Int_t isfiducialGen,ispassRec,isfiducialPreFSR,DYtautau, isBveto, nVtx;

  ClassDef ( ISRmumu_unfolding, 1);
};
#endif
