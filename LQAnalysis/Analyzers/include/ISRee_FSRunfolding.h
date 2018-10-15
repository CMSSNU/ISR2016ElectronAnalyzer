#ifndef ISRee_FSRunfolding_h
#define ISRee_FSRunfolding_h

#include "AnalyzerCore.h"
class ISRee_FSRunfolding : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  ISRee_FSRunfolding();
  ~ISRee_FSRunfolding();

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

  std::vector<Double_t> ptPreFSR,mPreFSR;
  std::vector<Double_t> ptPostFSR,mPostFSR;
  Double_t weight_, weightTotal;
  Int_t istriggered,issignal;

  ClassDef ( ISRee_FSRunfolding, 1);
};
#endif
