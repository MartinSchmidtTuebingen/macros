#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TString.h"
#include "TLegend.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TGraph.h"
#include "TMath.h"

#include <iostream>
#include <iomanip>

const Double_t nSigma = 2;

TCanvas* calculateSystematics(TString canvName, TString canvTitle, TH1D* histos[4], TH1D* finalHistos[4],
                              TString* systematicsHistosName, Int_t reference, TH1D* hSystematics[4], TGraphAsymmErrors* gr[4])
{
  Double_t ymax = 0;
  Double_t ymin = 0;
  
  for (Int_t j = 0; j < 4; j++) {
    hSystematics[j] = (TH1D*)histos[j]->Clone(systematicsHistosName[j].Data());
    hSystematics[j]->Reset(); 
    hSystematics[j]->GetXaxis()->SetRange(0, -1);
    
    for (Int_t bin = 1; bin <= histos[j]->GetNbinsX(); bin++) {
      hSystematics[j]->SetBinContent(bin, histos[reference]->GetBinContent(bin) - histos[j]->GetBinContent(bin));
      hSystematics[j]->SetBinError(bin, TMath::Sqrt(TMath::Abs(TMath::Power(histos[reference]->GetBinError(bin), 2) - 
                                                               TMath::Power(histos[j]->GetBinError(bin), 2))));
  
      if (hSystematics[j]->GetBinError(bin) == 0)
        hSystematics[j]->SetBinError(bin, 1e-10);
      Double_t temp = hSystematics[j]->GetBinContent(bin) + hSystematics[j]->GetBinError(bin);
      if (temp > ymax)
        ymax = temp;
        
      temp = hSystematics[j]->GetBinContent(bin) - hSystematics[j]->GetBinError(bin);
      if (temp < ymin)
        ymin = temp;
    }
  }
  
  TCanvas* canv = new TCanvas(canvName.Data(), canvTitle.Data(),100,10,1200,800);
  canv->SetGridy(1);
  
  hSystematics[reference]->Draw("e p");
  hSystematics[reference]->GetYaxis()->SetRangeUser(ymin, ymax);
  hSystematics[0]->Draw("e p same");
  hSystematics[1]->Draw("e p same");
  hSystematics[3]->Draw("e p same");
  
  TLegend* legend = new TLegend(0.622126, 0.605932, 0.862069, 0.855932);    
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->AddEntry(hSystematics[2], "#Delta_{#pi}", "p");
  legend->AddEntry(hSystematics[0], "#Delta_{e}", "p");
  legend->AddEntry(hSystematics[1], "#Delta_{K}", "p");
  legend->AddEntry(hSystematics[3], "#Delta_{p}", "p");
  legend->Draw();
			       
			       
  const Int_t nBins = 68; //histos[reference]->GetNbinsX();
  Double_t x[nBins];
  Double_t y[nBins];
  Double_t xerr[nBins];
  Double_t yerrl[nBins];
  Double_t yerrh[nBins];
  
  for (Int_t bin = 0; bin < nBins; bin++) {
    x[bin] = histos[reference]->GetBinCenter(bin + 1);
    xerr[bin] = histos[reference]->GetBinWidth(bin + 1) / 2.;
    y[bin] = finalHistos[reference]->GetBinContent(bin + 1);
    
    // Take all points that are more than nSigma sigma away from 0.
    // If there are at least 2 such points, take the difference between
    // the extreme values (i.e. maximum and minimum) as a measure of
    // the systematics
    Int_t count = 0;
    Double_t deltaMin = 0;
    Double_t deltaMax = 0;
    
    for (Int_t j = 0; j < 4; j++) {
      if (hSystematics[j]->GetBinError(bin + 1) == 0)
        continue;
        
      Double_t delta = hSystematics[j]->GetBinContent(bin + 1);
      if (TMath::Abs(delta / hSystematics[j]->GetBinError(bin + 1)) > nSigma)  {
        if (count == 0) {
          deltaMin = delta;
          deltaMax = delta;
        }
        else {
          if (delta < deltaMin)
            deltaMin = delta;
          if (delta > deltaMax)
            deltaMax = delta;
        }
        count++;
      }
    }
    
    if (deltaMax > 0.) 
      yerrh[bin] = deltaMax;
    else
      yerrh[bin] = 0.;
      
    if (deltaMin < 0.)
      yerrl[bin] = -deltaMin;
    else
      yerrl[bin] = 0.;
    /*
    if (count < 2)
      yerr[bin] = 0.;
    else
      yerr[bin] = (deltaMax - deltaMin) / TMath::Sqrt(2);
    */
  }
  
  gr[reference] = new TGraphAsymmErrors(nBins, x, y, xerr, xerr, yerrl, yerrh);
  gr[reference]->SetFillColor(kGray);
  gr[reference]->SetFillStyle(3004 + reference);
  
  return canv;
}

TCanvas* DrawFractionHistos(TString canvName, TString canvTitle, Double_t pLow, Double_t pHigh, TH1D* hist[4],
                            TGraphAsymmErrors* gr[4])
{
  TCanvas* canv = new TCanvas(canvName.Data(), canvTitle.Data(),100,10,1200,800);
  canv->SetGridy(1);
  hist[2]->GetYaxis()->SetRangeUser(0.0, 1.0);
  hist[2]->GetYaxis()->SetTitle(canvTitle.Data());
  hist[2]->GetXaxis()->SetRangeUser(pLow - 0.1, pHigh + 0.1);
  hist[2]->SetMarkerStyle(20);
  hist[2]->SetFillStyle(3004 + 2);
  hist[2]->SetFillColor(kGray);
  hist[2]->Draw("e p");
  hist[0]->GetYaxis()->SetTitle(canvTitle.Data());
  hist[0]->GetXaxis()->SetRangeUser(pLow - 0.1, pHigh + 0.1);
  hist[0]->SetMarkerStyle(21);
  hist[0]->SetFillStyle(3004 + 0);
  hist[0]->SetFillColor(kGray);
  hist[0]->Draw("e p same");
  hist[1]->GetYaxis()->SetTitle(canvTitle.Data());
  hist[1]->GetXaxis()->SetRangeUser(pLow - 0.1, pHigh + 0.1);
  hist[1]->SetMarkerStyle(22);
  hist[1]->SetFillStyle(3004 + 1);
  hist[1]->SetFillColor(kGray);
  hist[1]->Draw("e p same");
  hist[3]->GetYaxis()->SetTitle(canvTitle.Data());
  hist[3]->GetXaxis()->SetRangeUser(pLow - 0.1, pHigh + 0.1);
  hist[3]->SetMarkerStyle(29);
  hist[3]->SetFillStyle(3004 + 3);
  hist[3]->SetFillColor(kGray);
  hist[3]->Draw("e p same");
  
  gr[0]->Draw("2 same");
  gr[1]->Draw("2 same");
  gr[2]->Draw("2 same");
  gr[3]->Draw("2 same");
  
  TLegend* legend = new TLegend(0.622126, 0.605932, 0.862069, 0.855932);    
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->AddEntry(hist[2], "#pi", "fp");
  legend->AddEntry(hist[0], "e", "fp");
  legend->AddEntry(hist[1], "K", "fp");
  legend->AddEntry(hist[3], "p", "fp");
  legend->Draw();
  
  return canv;
}

TCanvas* loadCanvas(TFile* f, TString canvName) 
{
  TCanvas* c = 0x0;
  c = dynamic_cast<TCanvas*>(f->Get(canvName.Data()));
  if (!c)
    std::cout << "Failed to load canvas \"" << canvName.Data() << "\"!" << std::endl;
    
  return c;
}

Bool_t loadHistos(TH1D* histos[4], TString* histNames, TCanvas* canv)
{
  for (Int_t i = 0; i < 4; i++) {
    histos[i] = dynamic_cast<TH1D*>(canv->FindObject(histNames[i].Data()));
    if (!histos[i]) {
      std::cout << "Failed to load histo \"" << histNames[i].Data() << "\" in canvas \"" << canv->GetName() << "\"!" << std::endl;
      return kFALSE;
    } 
  }
  
  return kTRUE;
}

Int_t Systematics(TString fileName) { 

	TFile* f = TFile::Open(fileName.Data());
  if (!f)  {
    std::cout << "Failed to open file \"" << fileName.Data() << "\"!" << std::endl;
    return -1;
  }
  
  // Extract the data histogram
  TCanvas* cFractions = loadCanvas(f, "cFractions");
  if (!cFractions)
    return -1;
    
  TCanvas* cPions = loadCanvas(f, "cFractionsPions");
  if (!cPions)
    return -1;
    
  TCanvas* cElectrons = loadCanvas(f, "cFractionsElectrons");
  if (!cElectrons)
    return -1;
  
  TCanvas* cKaons = loadCanvas(f, "cFractionsKaons");
  if (!cKaons)
    return -1;
  
  TCanvas* cProtons = loadCanvas(f, "cFractionsProtons");
  if (!cProtons)
    return -1;
    
  
  TH1D* hFractionElectrons[4] = {0, 0, 0, 0};
  TString histNameElectrons[4] = {"hFractionElectronsDeltaElectron", "hFractionElectronsDeltaKaon", "hFractionElectronsDeltaPion",
                                  "hFractionElectronsDeltaProton" };
  if (!loadHistos(hFractionElectrons, histNameElectrons, cElectrons))
    return -1;
    
  TH1D* hFractionKaons[4] = {0, 0, 0, 0};
  TString histNameKaons[4] = {"hFractionKaonsDeltaElectron", "hFractionKaonsDeltaKaon", "hFractionKaonsDeltaPion",
                              "hFractionKaonsDeltaProton" };
  if (!loadHistos(hFractionKaons, histNameKaons, cKaons))
    return -1;
  
  TH1D* hFractionPions[4] = {0, 0, 0, 0};
  TString histNamePions[4] = {"hFractionPionsDeltaElectron", "hFractionPionsDeltaKaon", "hFractionPionsDeltaPion",
                              "hFractionPionsDeltaProton" };
  if (!loadHistos(hFractionPions, histNamePions, cPions))
    return -1;
  
  TH1D* hFractionProtons[4] = {0, 0, 0, 0};
  TString histNameProtons[4] = {"hFractionProtonsDeltaElectron", "hFractionProtonsDeltaKaon", "hFractionProtonsDeltaPion",
                                "hFractionProtonsDeltaProton" };
  if (!loadHistos(hFractionProtons, histNameProtons, cProtons))
    return -1;
    
    
  TH1D* hFractions[4] = {0, 0, 0, 0};
  TString histNameAll[4] = {"hFractionElectrons", "hFractionKaons", "hFractionPions", "hFractionProtons" };
  if (!loadHistos(hFractions, histNameAll, cFractions))
    return -1;
  
  
  
  TH1D* hSystematicsPions[4] = {0, 0, 0, 0};
  TGraphAsymmErrors* grSysErrors[4] = {0, 0, 0, 0 };
  TString systematicsHistosNamePions[4] = {"hSystematicsPionsDeltaElectron", "hSystematicsPionsDeltaKaon", "hSystematicsPionsDeltaPion", 
                                           "hSystematicsPionsDeltaProton" };
  
  TCanvas* cSystematicsPions = calculateSystematics("cSystematicsPions", "Systematics Pions", hFractionPions, hFractions, 
                                                    systematicsHistosNamePions, 2, hSystematicsPions, grSysErrors);
  
  TH1D* hSystematicsElectrons[4] = {0, 0, 0, 0};
  TString systematicsHistosNameElectrons[4] = {"hSystematicsElectronsDeltaElectron", "hSystematicsElectronsDeltaKaon", 
                                               "hSystematicsElectronsDeltaPion", "hSystematicsElectronsDeltaProton" };
  
  TCanvas* cSystematicsElectrons = calculateSystematics("cSystematicsElectrons", "Systematics Electrons", hFractionElectrons, hFractions,
                                                        systematicsHistosNameElectrons, 0, hSystematicsElectrons,
                                                        grSysErrors);
  
  
  TH1D* hSystematicsKaons[4] = {0, 0, 0, 0};
  TString systematicsHistosNameKaons[4] = {"hSystematicsKaonsDeltaElectron", "hSystematicsKaonsDeltaKaon", "hSystematicsKaonsDeltaPion", 
                                           "hSystematicsKaonsDeltaProton" };

  TCanvas* cSystematicsKaons = calculateSystematics("cSystematicsKaons", "Systematics Kaons", hFractionKaons, hFractions,
                                                    systematicsHistosNameKaons, 1, hSystematicsKaons, grSysErrors);
  
  
  TH1D* hSystematicsProtons[4] = {0, 0, 0, 0};
  TString systematicsHistosNameProtons[4] = {"hSystematicsProtonsDeltaElectron", "hSystematicsProtonsDeltaKaon", 
                                             "hSystematicsProtonsDeltaPion", "hSystematicsProtonsDeltaProton" };
  
  TCanvas* cSystematicsProtons = calculateSystematics("cSystematicsProtons", "Systematics Protons", hFractionProtons, hFractions,
                                                      systematicsHistosNameProtons, 3, hSystematicsProtons, grSysErrors);
  
  
  
  
  Double_t pLow = 0.0;
  Double_t pHigh = 50.;
  DrawFractionHistos("cFractionsWithSystematicError", "Particle fractions", pLow, pHigh, hFractions, grSysErrors);
  
  
  return 0;
  /*
  //TODO
  // Save final results
  saveF->cd();
  
  if (fElectronFraction)
    fElectronFraction->Write();
  
  if (hFractionElectrons)
    hFractionElectrons->Write(0, TObject::kWriteDelete);
    
  if (hFractionKaons)
    hFractionKaons->Write(0, TObject::kWriteDelete);
  
  if (hFractionPions)
    hFractionPions->Write(0, TObject::kWriteDelete);
  
  if (hFractionProtons)
    hFractionProtons->Write(0, TObject::kWriteDelete);
    
  if (cFractions)
    cFractions->Write();
  if (cFractionsPions)
    cFractionsPions->Write();
  if (cFractionsElectrons)
    cFractionsElectrons->Write();
  if (cFractionsKaons)
    cFractionsKaons->Write();
  if (cFractionsProtons)
    cFractionsProtons->Write();
  
  saveF->Close();
   
  return 0; */
}
