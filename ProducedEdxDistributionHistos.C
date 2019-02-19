#include "TString.h"
#include "TObjArray.h"
#include "TFile.h"
#include "TH2.h"
#include "TH1.h"
#include "TAxis.h"
#include "TObjString.h"
#include "TDirectoryFile.h"
#include "TH3.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TFitResultPtr.h"
#include "FitAsymmetricShape.C"
#include "TCanvas.h"

#include <iostream>

//     Double_t binsYMult[nBinsY[2] + 1] = {0, 100, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800 }; - y-Axis
//     Double_t binsYclusters[nBinsY[3] + 1] = {60, 80, 100, 120, 140, 160}; z-Axis

const Int_t nOfArrays = 4;
const TString arrayNames[nOfArrays] = {"Pi first", "Pi second", "Ka", "Pr"};
const TString legendNames[nOfArrays] = {"#pi first", "#pi second", "Kaons", "Protons"};
// const TString legendNames[3*nOfArrays] = {"#pi first Mult 0-30", "#pi first Mult 30-80", "#pi first Mult 80-200", "#pi second Mult 0-30", "#pi second Mult 30-80", "#pi second Mult 80-200", "Kaons Mult 0-30", "Kaons Mult 30-80", "Kaons Mult 80-200", "Protons Mult 0-30", "Protons Mult 30-80", "Protons Mult 80-200"};

const Int_t nOfHistos = 3;
const TString histNames[nOfHistos] = {"fh2Eta", "fh2Phi", "fh3MultClusters"};

const Int_t speciescolor[4] = {kBlack, kRed, kGreen, kBlue};

const Int_t nMultProj = 9;
const Int_t binsForMultProjection[2*nMultProj] = {1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9};
// const Int_t binsForMultProjection[2*nMultProj] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 23, 24, 26, 27, 29, 30, 32, 33, 35, 36, 40};
// const Int_t binsForMultProjection[2*nMultProj] = {1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14, 14, 15, 15, 16, 16, 17, 17, 18, 18, 19, 19, 20, 20, 21, 21, 22, 22, 23, 23, 24, 24, 25, 25, 26, 26, 27, 27, 28, 29, 30, 30, 31, 31, 32, 32, 33, 33, 34, 34, 35, 35, 36, 36, 37, 37, 38, 38, 39, 39, 40, 40};
// const Int_t binsForMultProjection[2*nMultProj] = {1, 40};
const Bool_t kSaveMultHistos = kTRUE;

const Int_t nClustProj = 8;
Int_t binsForClustProjection[2*nClustProj] = {1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8/*, 9, 9, 10, 11, 20 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14, 14, 15, 15, 16, 16, 17, 17, 18, 18, 19, 19, 20, 20*/ };
// Int_t binsForClustProjection[2*nClustProj] = {7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14, 14, 15, 15 };
// Int_t binsForClustProjection[2*nClustProj] = {14,15};
const Bool_t kSaveClustHistos = kTRUE;

void ProcessArray(TObjArray* histarray, TFile* output);
void ProjectHistos(TDirectoryFile* df, TH2* hist, TString namearray);
void ProjectHistos(TDirectoryFile* df, TH3* hist, TString namearray);

TCanvas* lambdacanvas_mult = 0x0;
TCanvas* meancanvas_mult = 0x0;

void ProducedEdxDistributionHistos(TString filename, TString outfile) {
  
  gStyle->SetOptStat("ne");
  gROOT->LoadMacro("./FitAsymmetricShape.C+");
  
  TFile* f = new TFile(filename.Data(),"READ");
  
  if (!f)
    return;
  
  TObjArray* arr = 0x0;
  
  TObjArray* tokens = filename.Tokenize("/");
  
  TObjString* objString = (TObjString*)tokens->Last();
  
//   TString arrayName = objString->GetString().ReplaceAll(".root","");
  
  TString arrayName = "maschmid_PID_Jets_Inclusive";
  
  f->GetObject(arrayName.Data(),arr);
  
  TFile* output = new TFile(outfile.Data(),"RECREATE");
  
  for (Int_t i=0;i<nOfArrays;++i) {
    TObjArray* histarray = (TObjArray*)arr->FindObject(arrayNames[i]);
    ProcessArray(histarray,output);
  }
  
  output->cd();
  meancanvas_mult = new TCanvas("cmean","cmean");
  meancanvas_mult->SetGridx(kTRUE);
  meancanvas_mult->SetGridy(kTRUE);
  TLegend* legend = new TLegend(0.1, 0.9, 0.4, 0.7);
  TH1F* h = 0x0;
  for (Int_t i=0;i<nOfArrays;++i) {
    h = (TH1F*)output->FindObjectAny(TString::Format("h1meanMult_140_160_" + arrayNames[i]).Data());
    if (h) {
      if (i == 0) {
        h->SetTitle("Mean d#it{E}/d#it{x}"); 
        h->SetLineColor(kBlack);
        h->SetMarkerColor(kBlack);
      }
      h->GetYaxis()->SetRangeUser(0.97,1.01);
      legend->AddEntry(h, legendNames[i]);
      h->SetStats(0);
      h->Draw(i == 0 ? "" : "same");
    }
/*    h = (TH1F*)output->FindObjectAny(TString::Format("h1mean_Mult_30_80_" + arrayNames[i]).Data());
    if (h) {
      h->SetLineColor(i == 0 ? kGray+2 : h->GetLineColor() + 1);
      h->SetMarkerColor(i == 0 ? kGray+2 : h->GetMarkerColor() + 1);
      h->GetYaxis()->SetRangeUser(0.97,1.01);
      legend->AddEntry(h, legendNames[3*i+1]);
      h->SetStats(0);
      h->Draw("same");
    }
    h = (TH1F*)output->FindObjectAny(TString::Format("h1mean_Mult_80_200_" + arrayNames[i]).Data());
    if (h) {
      h->SetLineColor(i == 0 ? kBlack : h->GetLineColor() + 2);
      h->SetMarkerColor(i == 0 ? kBlack : h->GetMarkerColor() + 2);
      h->GetYaxis()->SetRangeUser(0.97,1.01);
      legend->AddEntry(h, legendNames[3*i+2]);
      h->SetStats(0);
      h->Draw("same");
    }  */      
  }
  legend->Draw("same");
  meancanvas_mult->Update();
  TPaveText* paveTextTitle = (TPaveText*)meancanvas_mult->FindObject("title");
  if (paveTextTitle) {
    paveTextTitle->Clear();
    paveTextTitle->AddText("dE/dx/dEdx_expected");
  }
  meancanvas_mult->Modified();  
  meancanvas_mult->Write();
  meancanvas_mult->Close();
  delete legend;
//   
//   lambdacanvas_mult = new TCanvas("clambda","clambda");
//   lambdacanvas_mult->SetGridx(kTRUE);
//   lambdacanvas_mult->SetGridy(kTRUE);
//   legend = new TLegend(0.1, 0.9, 0.4, 0.7);
//   for (Int_t i=0;i<nOfArrays;++i) {
//     h = (TH1F*)output->FindObjectAny(TString::Format("h1lambda_nearTracks_120_160_" + arrayNames[i]).Data());
//     if (h) {
//       h->GetYaxis()->SetRangeUser(1.5,3.5);
//       legend->AddEntry(h, legendNames[i]);
//       h->SetStats(0);
//       h->Draw(i == 0 ? "" : "same");
//     }
//   }  
//   legend->Draw("same");
//   lambdacanvas_mult->Update();
//   paveTextTitle = (TPaveText*)lambdacanvas_mult->FindObject("title");
//   if (paveTextTitle) {
//     paveTextTitle->Clear();
//     paveTextTitle->AddText("Mean dE/dx for TPC cluster number 120 - 150");
//   }
//   
//   lambdacanvas_mult->Modified();    
//   lambdacanvas_mult->Write();
//   lambdacanvas_mult->Close();  
//   delete meancanvas_mult;
//   meancanvas_mult = 0x0;
//   delete lambdacanvas_mult;
//   lambdacanvas_mult = 0x0;
  f->Close();
  output->Close();
  
}

void ProcessArray(TObjArray* histarray, TFile* output) {
  output->cd();
  TDirectoryFile* df = new TDirectoryFile(histarray->GetName(),histarray->GetName());
  df->Write();
  df->cd();
  for (Int_t i=0;i<nOfHistos;++i) {
    if (i != 2) {
      TH2* hist = (TH2*)histarray->At(i);
//       ProjectHistos(df, hist, TString(histarray->GetName()));
    }
    else {
      TH3* hist = (TH3*)histarray->At(i);
      ProjectHistos(df, hist, TString(histarray->GetName()));
    }
  }
  return;
}


void ProjectHistos(TDirectoryFile* df, TH3* hist, TString namearray) {
  
  const Int_t lambdainitialvalue = 2.0;
  const Int_t lambdalowerbound = 0.2;
  const Int_t lambdaupperbound = 10.0;
  
  static Int_t species = 0;
    
  df->cd();
  TString histName = hist->GetName();
  TH1* projHisto = 0x0;
  TString name;
  
  Double_t* multBins = new Double_t[nMultProj+1];
  Double_t* clustBins = new Double_t[nClustProj+1];
  
  for (Int_t i=0;i<nClustProj;i++) {
    clustBins[i] = hist->GetZaxis()->GetBinLowEdge(binsForClustProjection[2*i]);
  }
  clustBins[nClustProj] = hist->GetZaxis()->GetBinUpEdge(binsForClustProjection[2*nClustProj-1]);
  
  for (Int_t i=0;i<nMultProj;i++) {
    multBins[i] = hist->GetYaxis()->GetBinLowEdge(binsForMultProjection[2*i]);
  }
  multBins[nMultProj] = hist->GetYaxis()->GetBinUpEdge(binsForMultProjection[2*nMultProj-1]);
  
  TH1F** h1lambda_forMultClass = new TH1F*[nMultProj];
  TH1F** h1mean_forMultClass = new TH1F*[nMultProj];  
  
  for (Int_t i=0;i<nMultProj;i++) {
    TString multName = TString::Format("Mult_%i_%i_", (Int_t)multBins[i],(Int_t)multBins[i+1]);
    h1lambda_forMultClass[i] = new TH1F(("h1lambda_" + multName + namearray).Data(),("#lambda_" + multName + namearray).Data(),nClustProj,clustBins);
    h1lambda_forMultClass[i]->SetXTitle("Near Tracks");
    h1mean_forMultClass[i] = new TH1F(("h1mean_" + multName + namearray).Data(),("Mean_" + multName + namearray).Data(),nClustProj,clustBins);
    h1mean_forMultClass[i]->SetXTitle("Multiplicity"); 
    
    h1lambda_forMultClass[i]->SetMarkerStyle(20);
    h1lambda_forMultClass[i]->SetDrawOption("e");
    h1lambda_forMultClass[i]->SetLineColor(speciescolor[species]);
    h1lambda_forMultClass[i]->SetMarkerColor(speciescolor[species]);
    
    h1mean_forMultClass[i]->SetMarkerStyle(20);
    h1mean_forMultClass[i]->SetDrawOption("e");
    h1mean_forMultClass[i]->SetLineColor(speciescolor[species]);
    h1mean_forMultClass[i]->SetMarkerColor(speciescolor[species]);        
  }

  for (Int_t j=0;j<nClustProj;++j) {
    
    Double_t clustlowedge = hist->GetZaxis()->GetBinLowEdge(binsForClustProjection[2*j]);
    Double_t clustupedge = hist->GetZaxis()->GetBinUpEdge(binsForClustProjection[2*j+1]);
    TString clustName = TString::Format("Eta_%.1f_%.1f_", clustlowedge,clustupedge);
    
    TDirectoryFile* dfclust = 0x0;
    
    if (nClustProj > 1) {
      dfclust = new TDirectoryFile(clustName.Data(), clustName.Data());
      dfclust->cd();
    }
    
    TH1F* h1lambda = new TH1F(("h1lambda_" + clustName + namearray).Data(),("#lambda_" + clustName + namearray).Data(),nMultProj,multBins);
    h1lambda->SetXTitle("Multiplicity");
    
    TH1F* h1mean = new TH1F(("h1mean" + clustName + namearray).Data(),("Mean_" + clustName + namearray).Data(),nMultProj,multBins);
    h1mean->SetXTitle("Multiplicity");
    
    h1lambda->SetMarkerStyle(20);
    h1lambda->SetDrawOption("e");
    h1lambda->SetLineColor(speciescolor[species]);
    h1lambda->SetMarkerColor(speciescolor[species]);
    
    h1mean->SetMarkerStyle(20);
    h1mean->SetDrawOption("e");
    h1mean->SetLineColor(speciescolor[species]);
    h1mean->SetMarkerColor(speciescolor[species]);    
    
    for (Int_t i=0;i<nMultProj;++i) {
      name = TString::Format("fh1_Mult_%i_%i_%s",(Int_t)hist->GetYaxis()->GetBinLowEdge(binsForMultProjection[2*i]),(Int_t)hist->GetYaxis()->GetBinUpEdge(binsForMultProjection[2*i+1]),clustName.Data());
      projHisto = hist->ProjectionX(name.Data(),binsForMultProjection[2*i], binsForMultProjection[2*i+1], binsForClustProjection[2*j], binsForClustProjection[2*j+1]);
      Double_t nOfEntries = projHisto->Integral();
      if (nOfEntries < 20)
        continue;
      
      projHisto->SetXTitle("dE/dx / dE/dx_expected");
      projHisto->SetYTitle("Entries");   
      projHisto->SetTitle((namearray + " " + name.ReplaceAll("fh1_","")).Data()); 
//       TFitResultPtr fitResult = 0;
      TFitResultPtr fitResult = FitAsymmetricShapeHisto(projHisto, lambdainitialvalue, lambdalowerbound, lambdaupperbound, kFALSE);
      TF1* gausFunction = 0x0;
      TF1* asymFunction = 0x0;
      if (fitResult && fitResult != 4) {    
        TCanvas* c = new TCanvas(("c_" + namearray + "_" + name).Data(), name.Data(),  100,10,1380,800);
        c->Divide(1,1);
        c->GetPad(1)->SetGridx(kTRUE);
        c->GetPad(1)->SetGridy(kTRUE);
        c->GetPad(1)->SetLogy(kTRUE);
        c->cd(1);
        projHisto->Draw("e");
        c->Write();
        c->Close();
        gausFunction = (TF1*)projHisto->GetListOfFunctions()->At(0);
        asymFunction = (TF1*)projHisto->GetListOfFunctions()->At(1);
        delete c;
        c = 0;        
      }
      
      projHisto->Write();
      
      if (asymFunction) {
        h1lambda->SetBinContent(i+1,asymFunction->GetParameter(3));
        h1lambda->SetBinError(i+1,asymFunction->GetParError(3));
        h1mean->SetBinContent(i+1,asymFunction->GetParameter(1));
        h1mean->SetBinError(i+1,asymFunction->GetParError(1));
        h1lambda_forMultClass[i]->SetBinContent(j+1,asymFunction->GetParameter(3));
        h1lambda_forMultClass[i]->SetBinError(j+1,asymFunction->GetParError(3));
        h1mean_forMultClass[i]->SetBinContent(j+1,asymFunction->GetParameter(1));
        h1mean_forMultClass[i]->SetBinError(j+1,asymFunction->GetParError(1));        
      }
      
      projHisto->Scale(1.0/projHisto->Integral());
      projHisto->GetListOfFunctions()->Clear();
      for (Int_t m=0;m<=projHisto->GetXaxis()->GetNbins();++m)
        projHisto->SetBinError(m,projHisto->GetBinError(m)/nOfEntries);
      fitResult = FitAsymmetricShapeHisto(projHisto, lambdainitialvalue, lambdalowerbound, lambdaupperbound, kFALSE);
      projHisto->SetName(("fh1"+name+"_scaled").Data());
      projHisto->SetYTitle("Entries scaled");
      if (fitResult && fitResult != 4) {
        TCanvas* c = new TCanvas(("c_" + namearray + "_" + name + "_scaled").Data(), name.Data(),  100,10,1380,800);
        c->Divide(1,1);
        c->GetPad(1)->SetGridx(kTRUE);
        c->GetPad(1)->SetGridy(kTRUE);
        c->GetPad(1)->SetLogy(kTRUE);
        c->cd(1);
        projHisto->Draw("e");
        c->Write(); 
        c->Close();
        delete c;
        c = 0;
      }
      projHisto->Write();
      delete projHisto;
      projHisto = 0x0;
    }
    if (kSaveClustHistos) {
      h1lambda->Write();
    }
    delete h1lambda;
    h1lambda = 0x0;
    if (kSaveClustHistos) {
      h1mean->Write();
    }
    df->cd();
    delete h1mean;
    h1mean = 0x0;
    if (dfclust)
      dfclust->Write();
    df->cd();
  }
  for (Int_t i=0;i<nMultProj;++i) {
    if (kSaveMultHistos) {
      h1lambda_forMultClass[i]->Write();
    }
    delete h1lambda_forMultClass[i];
    h1lambda_forMultClass[i] = 0x0;
    if (kSaveMultHistos) {
      h1mean_forMultClass[i]->Write();
    }
    delete h1mean_forMultClass[i];
    h1mean_forMultClass[i] = 0x0;
  }
  species++;
//   Int_t nEntries = hist->Integral();
//   TH1F* clusterFractions = new TH1F("fh1clustFractions", "Fractions of Entries in Cluster Bin on Total Entries", hist->GetZaxis()->GetNbins(),hist->GetZaxis()->GetXbins()->GetArray());
//   Int_t xlast = hist->GetXaxis()->GetLast();
//   Int_t ylast = hist->GetXaxis()->GetLast();
//   for (Int_t i=1;i<=hist->GetZaxis()->GetNbins();++i) {
//     clusterFractions->SetBinContent(i,hist->Integral(1,xlast,1,ylast,i,i)/nEntries);
//   }
//   clusterFractions->Write();
//   delete clusterFractions;
//   clusterFractions = 0x0;
  return;
}

void ProjectHistos(TDirectoryFile* df, TH2* hist, TString namearray) {
  df->cd();
  TString histName = hist->GetName();
  TAxis* axe = hist->GetYaxis();
  Int_t nBins = axe->GetNbins();
  TH1* projHisto = 0x0;
  TString name;
  for (Int_t binNumber = 1;binNumber<=nBins;++binNumber) {
    name = TString::Format("%s_%1.1f_%1.1f",(histName.ReplaceAll("fh2","fh1")).Data(),(Double_t)axe->GetBinLowEdge(binNumber),(Double_t)axe->GetBinUpEdge(binNumber));
    projHisto = hist->ProjectionX(name.Data(),binNumber, binNumber);
    projHisto->SetXTitle("dE/dx / dE/dx_expected");
    projHisto->SetYTitle("Entries");
    projHisto->SetTitle((namearray + " " + name.ReplaceAll("fh1","")).Data()); 
    projHisto->Write();
    Double_t nOfEntries = projHisto->Integral();
    if (nOfEntries > 0) {
      projHisto->Scale(1.0/projHisto->Integral());
      projHisto->SetName((name+"_scaled").Data());
      projHisto->SetYTitle("Entries scaled");
      projHisto->Write();
    }
  }
  return;
}     