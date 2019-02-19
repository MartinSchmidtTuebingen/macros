#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TGraph.h"
#include "TF1.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TMath.h"
#include "TString.h"
#include "TStyle.h"

// #include "AliPID.h"

#include <iostream>
#include <iomanip>

void CalculateMiddleValueUnderlyingEvent(TString path) {
  
  const Int_t numJetPtBins = 4;
  const Int_t jetPt[numJetPtBins + 1] = { 5, 10, 15, 20, 30};
  
  const Int_t maxMode = 2;
  TString mode[maxMode] = {"Pt", "Z"};
  
  const Int_t nOfYieldHistograms = 5;
  
  const TString histNames[nOfYieldHistograms] = {"hFractionElectrons", "hFractionMuons", "hFractionPions", "hFractionKaons", "hFractionProtons" };
  const TString histNamesYields[nOfYieldHistograms] = {"hYieldElectrons", "hYieldMuons", "hYieldPions", "hYieldKaons", "hYieldProtons" };
  const TString histNamesToPiRatios[nOfYieldHistograms] = {"hRatioToPiElectrons", "hRatioToPiMuons", "", "hRatioToPiKaons", "hRatioToPiProtons" }; 
  
  const Int_t nOfPionHistogram = 2;
  
  TH1* yieldTotal = 0x0;
  
  TH1::SetDefaultSumw2();
  
  const Int_t nOfUEFiles = 2;
  TString methods[nOfUEFiles] = {"PC", "RC"};
  TString inputNames[nOfUEFiles];
  TFile* inputFiles[nOfUEFiles];
  
  TH1 *yieldHistogramsUE[nOfYieldHistograms][nOfUEFiles]; 
  TH1* yieldHistogramAverage[nOfYieldHistograms];
  TH1* yieldHistogramFractions[nOfYieldHistograms];
  TH1* ratioToPiHistograms[nOfYieldHistograms];
  
  TString begin = TString::Format("%s/maschmid_PID_Jets_UE_",path.Data());
  TString middle = "_results_LLFit__";
  TString end = "_2_reg1_regFac1.00_noMuons_idSpectra_jetPt";
  
  TString outputbegin = TString::Format("%s/maschmid_PID_Jets_UE_Average",path.Data());

  for (Int_t iMode = 0; iMode < maxMode; iMode++) {
    for (Int_t iJetPt = 0; iJetPt < numJetPtBins; iJetPt++) {
      for (Int_t iFileNumber = 0 ; iFileNumber < nOfUEFiles; ++iFileNumber) {
        
        //Open UE Files

        inputFiles[iFileNumber] = new TFile(TString::Format("%s%s%s%s%s%d.0_%d.0.root", begin.Data(), methods[iFileNumber].Data(), middle.Data(),mode[iMode].Data(),end.Data(),jetPt[iJetPt], jetPt[iJetPt + 1] ),"READ");
        
        if (!inputFiles[iFileNumber]) {
          Printf("File %s not found", Form("%s%s%s%s%s%d.0_%d.0.root", begin.Data(), methods[iFileNumber].Data(), middle.Data(),mode[iMode].Data(),end.Data(),jetPt[iJetPt], jetPt[iJetPt + 1] ));
          return;
        }
        
        //Retrieving Yield Histograms
        for (Int_t iSpecies = 0; iSpecies < nOfYieldHistograms; ++iSpecies) {
          inputFiles[iFileNumber]->GetObject(histNamesYields[iSpecies],yieldHistogramsUE[iSpecies][iFileNumber]);
        }
      }

      for (Int_t iSpecies = 0; iSpecies < nOfYieldHistograms; ++iSpecies) {
        yieldHistogramAverage[iSpecies] = dynamic_cast<TH1*>(yieldHistogramsUE[iSpecies][0]->Clone(histNamesYields[iSpecies]));
        for (Int_t iFileNumber = 1 ; iFileNumber < nOfUEFiles; ++iFileNumber) {
          yieldHistogramAverage[iSpecies]->Add(yieldHistogramsUE[iSpecies][iFileNumber]);
        }
        yieldHistogramAverage[iSpecies]->Scale(1.0/(Double_t)(nOfUEFiles));
      }
      yieldTotal = dynamic_cast<TH1*>(yieldHistogramAverage[0]->Clone("hYieldTotal"));
      for (Int_t iSpecies = 1; iSpecies < nOfYieldHistograms; ++iSpecies) {
        yieldTotal->Add(yieldHistogramAverage[iSpecies]);
      }
      
      for (Int_t iSpecies = 0; iSpecies < nOfYieldHistograms; ++iSpecies) {
        yieldHistogramFractions[iSpecies] = dynamic_cast<TH1*>(yieldHistogramAverage[iSpecies]->Clone(histNames[iSpecies]));
        yieldHistogramFractions[iSpecies]->Divide(yieldTotal);
        if (iSpecies != nOfPionHistogram) {
          ratioToPiHistograms[iSpecies] = dynamic_cast<TH1*>(yieldHistogramAverage[iSpecies]->Clone(histNamesToPiRatios[iSpecies]));
          ratioToPiHistograms[iSpecies]->Divide(yieldHistogramAverage[nOfPionHistogram]);
        }
      }
      
      //Open output File, write and Close

      TFile* output = new TFile(TString::Format("%s%s%s%s%d.0_%d.0.root", outputbegin.Data(), middle.Data(),mode[iMode].Data(),end.Data(),jetPt[iJetPt], jetPt[iJetPt + 1] ),"RECREATE");
      for (Int_t iSpecies = 0; iSpecies < nOfYieldHistograms; ++iSpecies) {
        yieldHistogramAverage[iSpecies]->Write();
        yieldHistogramFractions[iSpecies]->Write();
        if (iSpecies != nOfPionHistogram)
          ratioToPiHistograms[iSpecies]->Write();
      }
      yieldTotal->Write();
      output->Close();

      for (Int_t iFileNumber = 0 ; iFileNumber < nOfUEFiles; ++iFileNumber) {
        //Close UE Files
        inputFiles[iFileNumber]->Close();
      }
    }
  }
}