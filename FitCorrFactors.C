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
  
  TH1* h = (TH1*)f->FindObjectAny("hSingleTrackEfficiency_GF_p");
  //For electrons
//   const Int_t parts = 4;
//   Double_t cuts[parts-1] = {0.6,3.2,8.0};
//   Int_t nparameters[parts] = {7,5,3,2};
  
  
  //For protons
//   const Int_t parts = 6;
//   Double_t cuts[parts-1] = {0.4,1.6,2.5,8.0,12.0};
//   Int_t nparameters[parts] = {5,4,4,2,5,2};
    
  //For kaons
//   const Int_t parts = 5;
//   Double_t cuts[parts-1] = {0.4,1.2,6,15};
//   Int_t nparameters[parts] = {3,3,5,4,2};  
  
  //For pions
//   const Int_t parts = 6;
//   Double_t cuts[parts-1] = {0.8,1.6,3.0,10.0,12.0};
//   Int_t nparameters[parts] = {9,4,4,3,3,2};    
  
  PieceWisePoly* pwp = new PieceWisePoly(parts,cuts,nparameters,0,50,0x0,2);
  TF1* func = new TF1("func",pwp,0,50,pwp->GetNOfParam());  

  for (Int_t i=1;i<=h->GetNbinsX();++i) {
    if (h->GetBinContent(i) == 1 && h->GetBinError(i) == 0) {
      h->SetBinContent(i,0);
    }
  }
  
  TCanvas* c = new TCanvas("c");
  c->SetLogx(kTRUE);
  h->Fit(func,"I","same");
}

