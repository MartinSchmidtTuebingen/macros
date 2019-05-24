#include "TFile.h"
#include "AliCFContainer.h"
#include "AliPID.h"

#include "THnSparseDefinitions.h"

#include <iostream>
using namespace std;

Int_t iPt     = 0;
Int_t iMCid   = 0;
Int_t iEta    = 0;
Int_t iCharge = 0;
Int_t iMult   = 0;
Int_t iJetPt  = 0;
Int_t iZ = 0;
Int_t iXi = 0;
Int_t iDistance = 0;
Int_t ijT = 0;

Int_t iObsAxis = 0;

Int_t writeOutCorrectionFiles(TString effFile, TString pathNameData, TString outfilepath, TString addoutFileName) {
  TFile* fileEff = new TFile(effFile.Data());
  if (!fileEff) {
    printf("Failed to open efficiency file \"%s\"\n", effFile.Data());
    return -1;
  }
  
  AliCFContainer* data = (AliCFContainer*)(fileEff->Get("containerEff"));
  if (!data) {
    printf("Failed to load efficiency container!\n");
    return -1;
  }  
  
  Int_t nOfJetBins = 4;
  Double_t jetBinLimits[2*nOfJetBins] = {5.0,10.0,10.0,15.0,15.0,20.0,20.0,30.0};
  
  Int_t nOfEffSteps = 2;
  EffSteps usedEffSteps[nOfEffSteps] = {kStepGenWithGenCuts, kStepRecWithRecCutsPrimaries};
  TString dirNameEffSteps[nOfEffSteps] = {"Gen", "RecCuts"};
  
  // For backward compatibility:
  // Check whether "P_{T}" or "p_{T}" is used
  TString momentumString = "p";
  for (Int_t i = 0; i < data->GetNVar(); i++) {
    TString temp = data->GetVarTitle(i);
    if (temp.Contains("P_{")) {
      momentumString = "P";
      break;
    }
    else if (temp.Contains("p_{")) {
      momentumString = "p";
      break;
    }
  }
  
//   AliPID::ParticleName(species)
  
  iPt     = data->GetVar(Form("%s_{T} (GeV/c)", momentumString.Data()));
  iMCid   = data->GetVar("MC ID");
  iEta    = data->GetVar("#eta");
  iCharge = data->GetVar("Charge (e_{0})");
  iMult   = data->GetVar("Centrality Percentile");
  // Will be set later, if jet pT is restricted
  iJetPt  = 0;
  iZ = 0;
  iXi = 0;  
  
  iJetPt = data->GetVar(Form("%s_{T}^{jet} (GeV/c)", momentumString.Data()));
  iZ     = data->GetVar(Form("z = %s_{T}^{track} / %s_{T}^{jet}", momentumString.Data(), momentumString.Data()));
  iXi    = data->GetVar(Form("#xi = ln(%s_{T}^{jet} / %s_{T}^{track})", momentumString.Data(), momentumString.Data()));
  iDistance = data->GetVar("R");
  ijT = data->GetVar("j_{T} (GeV/c)");  
 
  
  TFile* outFile = new TFile((outfilepath + TString("/outCorrections_PythiaFastJet_") + addoutFileName + TString(".root")).Data());
//   outFile->mkdir(
//   
//   Int_t lowerJetPtBinLimit = -1;
//   Int_t upperJetPtBinLimit = -2;
//   Bool_t restrictJetPtAxis = kFALSE;
//   Double_t actualLowerJetPt = -1.;
//   Double_t actualUpperJetPt = -1.;
//   
//   if (lowerJetPt >= 0 && upperJetPt >= 0) {
//     // Add subtract a very small number to avoid problems with values right on the border between to bins
//     lowerJetPtBinLimit = dataRebinned->GetAxis(iJetPt, 0)->FindFixBin(lowerJetPt + 0.001);
//     upperJetPtBinLimit = dataRebinned->GetAxis(iJetPt, 0)->FindFixBin(upperJetPt - 0.001);
//   }
  
  
  
}