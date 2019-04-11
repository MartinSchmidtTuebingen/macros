#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"

#include "./PieceWisePoly.h"

#include <iostream>

using namespace std;

void FitCorrFactors(TString effFile = "") {
  TFile* f = new TFile(effFile.Data());
  if (f->IsZombie())
    return;
  
  TH1* h = (TH1*)f->FindObjectAny("hSingleTrackEfficiency_GF_e");
  const Int_t parts = 4;
  Double_t cuts[parts-1] = {0.6,3.2,8.0};
  Int_t nparameters[parts] = {7,5,3,2};
    
  PieceWisePoly* pwpelectrons = new PieceWisePoly(parts,cuts,nparameters,0,50,0x0,2);
  TF1* funcElectrons = new TF1("func",pwpelectrons,0,50,pwpelectrons->GetNOfParam());

  for (Int_t i=1;i<=h->GetNbinsX();++i) {
    if (h->GetBinContent(i) == 1 && h->GetBinError(i) == 0) {
      h->SetBinContent(i,0);
    }
  }
  
  h->Fit(funcElectrons);
  TCanvas* c = new TCanvas();
  c->SetLogx(kTRUE);
  h->Draw("same");
}

