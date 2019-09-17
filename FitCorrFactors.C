#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TString.h"

#include "TF1.h"
#include "./PieceWisePoly.cxx"

#include "AliCFContainer.h"
#include "AliCFEffGrid.h"

#include "THnSparseDefinitions.h"

#include <iostream>

using namespace std;

const Int_t speciescolors[5] = {kMagenta, kYellow, kRed, kGreen, kBlue};

Int_t FitCorrFactors(TString effFile, TString outputfile) {
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
  
  TF1* effFunctions[AliPID::kSPECIES] = {0x0};
  PieceWisePoly* polynoms[AliPID::kSPECIES] = {0x0};
  
  //For electrons
    const Int_t parts_e = 4;
    Double_t cuts_e[parts_e-1] = {0.6,3.2,8.0};
    Int_t nparameters_e[parts_e] = {7,5,3,2};
    PieceWisePoly* pwp_e = new PieceWisePoly(parts_e,cuts_e,nparameters_e,0,50,0x0,2);
    TF1* func_e = new TF1("func_e",pwp_e,0,50,pwp_e->GetNOfParam());
    polynoms[AliPID::kElectron] = pwp_e;
//     Double_t parameters_e[11] = {-1.22427e+00, 2.25003e+01, -9.00154e+01, 1.42536e+02, 1.98605e+00, -2.33708e+02, 1.74505e+02, -4.40750e-01, 1.43504e-01, -1.59226e-02, -6.14939e-04};
//     func_e->SetParameters(parameters_e);
    effFunctions[AliPID::kElectron] = func_e;
    
  //   For muons
    const Int_t parts_mu = 6;
    Double_t cuts_mu[parts_mu-1] = {0.8,1.6,3.0,10.0,12.0};
    Int_t nparameters_mu[parts_mu] = {9,4,4,3,3,2};   
    PieceWisePoly* pwp_mu = new PieceWisePoly(parts_mu,cuts_mu,nparameters_mu,0,50,0x0,2);
    TF1* func_mu = new TF1("func_mu",pwp_mu,0,50,pwp_mu->GetNOfParam());
    polynoms[AliPID::kMuon] = pwp_mu;
    effFunctions[AliPID::kMuon] = func_mu;     
    
    //For protons
    const Int_t parts_p = 6;
    Double_t cuts_p[parts_p-1] = {0.4,1.6,2.5,8.0,12.0};
    Int_t nparameters_p[parts_p] = {6,4,4,2,5,2};
    PieceWisePoly* pwp_p = new PieceWisePoly(parts_p,cuts_p,nparameters_p,0,50,0x0,2);
    TF1* func_p = new TF1("func_p",pwp_p,0,50,pwp_p->GetNOfParam());
    polynoms[AliPID::kProton] = pwp_p;
    effFunctions[AliPID::kProton] = func_p;   

    //For kaons
    const Int_t parts_k = 5;
    Double_t cuts_k[parts_k-1] = {0.4,1.2,6,15};
    Int_t nparameters_k[parts_k] = {3,3,5,4,2};  
    PieceWisePoly* pwp_k = new PieceWisePoly(parts_k,cuts_k,nparameters_k,0,50,0x0,2);
    TF1* func_k = new TF1("func_k",pwp_k,0,50,pwp_k->GetNOfParam());
    polynoms[AliPID::kKaon] = pwp_k;
    effFunctions[AliPID::kKaon] = func_k;

  //   For pions
    const Int_t parts_pi = 6;
    Double_t cuts_pi[parts_pi-1] = {0.8,1.6,3.0,10.0,12.0};
    Int_t nparameters_pi[parts_pi] = {9,4,4,3,3,2};   
    PieceWisePoly* pwp_pi = new PieceWisePoly(parts_pi,cuts_pi,nparameters_pi,0,50,0x0,2);
    TF1* func_pi = new TF1("func_pi",pwp_pi,0,50,pwp_pi->GetNOfParam());
    polynoms[AliPID::kPion] = pwp_pi;
    effFunctions[AliPID::kPion] = func_pi;                                                              

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
  
  Int_t iPt     = data->GetVar(Form("%s_{T} (GeV/c)", momentumString.Data()));
  Int_t iMCid   = data->GetVar("MC ID");
  Int_t iEta    = data->GetVar("#eta");
  Int_t iCharge = data->GetVar("Charge (e_{0})");
  Int_t iMult   = data->GetVar("Centrality Percentile");
  
  TFile* output = new TFile(outputfile.Data(),"RECREATE");
  TCanvas* c = new TCanvas("cEfficiencies","Efficiency x Acceptance");
  c->SetLogx(kTRUE);
  
  for (Int_t species=0;species<AliPID::kSPECIES;++species) {
    data->SetRangeUser(iMCid,species+1,species+1,kTRUE);
    for (Int_t chargeBin=1;chargeBin<=2;chargeBin++) {
      
      data->SetRangeUser(iCharge,chargeBin,chargeBin,kTRUE);
//       data->SetRangeUser(iCharge,-1,-1,kTRUE);
      AliCFEffGrid* effSingleTrack = new AliCFEffGrid("effSingleTrack", "Efficiency x Acceptance", *data);
      effSingleTrack->CalculateEfficiency(kStepRecWithRecCutsMeasuredObsPrimaries, kStepGenWithGenCuts);      
      TH1* h = (TH1D*)effSingleTrack->Project(iPt);
      h->SetNameTitle(Form("hSingleTrackEfficiency_%s_%d", AliPID::ParticleShortName(species), chargeBin),Form("%s for %d", AliPID::ParticleLatexName(species), chargeBin));
      h->SetLineColor(speciescolors[species] + 2*chargeBin - 2);
      h->SetMarkerColor(speciescolors[species] + 2*chargeBin -2);
      h->GetXaxis()->SetRangeUser(0.15,50);
      for (Int_t i=1;i<=h->GetNbinsX();++i) {
        if (h->GetBinContent(i) == 1 && h->GetBinError(i) == 0) {
          h->SetBinContent(i,0);
        }
      }      
//       TF1* func = (TF1*)effFunctions[species]->Clone();
//       
//       TF1* func = new TF1("func",pwp_pi,0,50,pwp_pi->GetNOfParam()); 
      TF1* func = new TF1("func",polynoms[species],0,50,polynoms[species]->GetNOfParam()); 
      func->SetLineColor(speciescolors[species] + 2*chargeBin - 2);
      h->Fit(func,"I","same");
      
      Int_t nOfParts = polynoms[species]->GetNParts();
      TString parameters = TString::Format("%d",nOfParts);
      Double_t* cuts = polynoms[species]->GetCuts();
      Int_t* nOfParameters = polynoms[species]->GetNParameters();
      for (Int_t part = 0;part<nOfParts - 1;++part) {
        parameters += TString::Format(";%f",cuts[part]);
      }
      for (Int_t part = 0;part<nOfParts;++part) {
        parameters += TString::Format(";%d",nOfParameters[part]);
      }      
      for (Int_t i=0;i<func->GetNpar();++i) {
        parameters += TString::Format(";%f",func->GetParameter(i));
      }
      TNamed* paramSave = new TNamed(TString::Format("fastSimulationParameters_%s_%d",AliPID::ParticleShortName(species),chargeBin-1).Data(),parameters.Data());
      paramSave->Write();
      h->Draw("same");
      h->Write();
//       delete h;
      h = 0x0;
      delete effSingleTrack;
      effSingleTrack = 0x0;
    }
  }
  c->Write();
  output->Close();
  

return 0;
}

