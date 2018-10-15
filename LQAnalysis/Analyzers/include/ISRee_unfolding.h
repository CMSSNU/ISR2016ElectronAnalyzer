#ifndef ISRee_unfolding_h
#define ISRee_unfolding_h

#include "AnalyzerCore.h"
class ISRee_unfolding : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  ISRee_unfolding();
  ~ISRee_unfolding();

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

  std::vector<Double_t> etaRec,ptRec,mRec;
  std::vector<Double_t> etaGen,ptGen,mGen;
  Double_t weight_, weightTotal;
  Int_t istriggered,issignal, DYtautau;

  ClassDef ( ISRee_unfolding, 1);
};
#endif
