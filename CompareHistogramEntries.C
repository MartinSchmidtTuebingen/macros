#include "TH1.h"
#include "TFile.h"
#include "TMath.h"
#include <iostream>
using namespace std;

void CompareHistogramEntries(TH1D *hbase, TH1* hcomp, TString output) {
  TH1D* result = new TH1D(*hbase);
  result->SetNameTitle("hRelDev", "Relative Deviation");
  TH1D* resultScaled = new TH1D(*hbase);
  resultScaled->SetNameTitle("hRelDevScaled","Relative Deviation scaled with hcomp"); 
  
  Int_t binhcomp = 0;
  
  for (Int_t iBin=1;iBin<=hbase->GetNbinsX()+1;++iBin) {
    result->SetBinContent(iBin,0.0);
    result->SetBinError(iBin,0.0);
    resultScaled->SetBinContent(iBin,0.0);
    resultScaled->SetBinError(iBin,0.0);
  }
  
  result->GetYaxis()->SetRangeUser(-1.0,1.0);
  resultScaled->GetYaxis()->SetRangeUser(-1.0,1.0);
  
  for (Int_t iBin=1;iBin<=hbase->GetNbinsX()+1;++iBin) {
    binhcomp = hcomp->GetXaxis()->FindFixBin(hbase->GetXaxis()->GetBinCenter(iBin));
    if (binhcomp == 0 || TMath::Abs(hcomp->GetBinContent(binhcomp)) <= 10e-3)
      continue;
    if (binhcomp == hcomp->GetNbinsX() + 1)
      break;
    

    if (hbase->GetBinContent(iBin) > 10e-6) {
      Double_t D = hcomp->GetBinContent(binhcomp);
      Double_t B = hbase->GetBinContent(iBin);
      Double_t DeltaD = hcomp->GetBinError(binhcomp);
      Double_t DeltaB = hbase->GetBinError(iBin);
      Double_t relDeviation = D/B - 1.0;
      result->SetBinContent(iBin,relDeviation); 
      result->SetBinError(iBin,DeltaD/B - DeltaB * (D/(B*B)));
      resultScaled->SetBinContent(iBin,relDeviation * (1.0-B));
      resultScaled->SetBinError(iBin,- DeltaB * (D/(B*B)) + DeltaD * (1.0/B - 1.0));
    }
    
  }
  
  TFile* f = new TFile(output.Data(),"RECREATE");
  hbase->Write();
  hcomp->Write();
  result->Write();
  resultScaled->Write();
  f->Close();
  return;
  
}