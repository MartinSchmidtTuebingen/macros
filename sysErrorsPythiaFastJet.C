#include "TROOT.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLine.h"
#include "TF1.h"
#include "TProfile.h"
#include "TPaveText.h"
#include "TMath.h"
#include "TPaletteAxis.h"

#include <iostream>
#include <fstream>

Bool_t  savePlots = kTRUE;
Bool_t  writeSysErrors = kTRUE;
TString saveDir = "files/";


// --------------------------------------------------

void setHistoStyleColor(TH1* hist, Int_t color=2){

  hist->SetLineColor(1);
  hist->SetMarkerStyle(20);
  hist->SetMarkerSize(0.7);
  hist->SetMarkerColor(color);
  hist->SetTitleSize(0.08);
  
  hist->GetXaxis()->SetTitleSize(0.06);
  hist->GetYaxis()->SetTitleSize(0.06); 
  hist->GetXaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetLabelSize(0.05);
  hist->GetXaxis()->SetTitleOffset(1.0); 
  hist->GetYaxis()->SetTitleOffset(0.7); 

  hist->GetXaxis()->SetTitleFont(62); 
  hist->GetYaxis()->SetTitleFont(62); 
  hist->GetXaxis()->SetLabelFont(62); 
  hist->GetYaxis()->SetLabelFont(62); 


  /*
  hist->GetXaxis()->SetTitleSize(0.04);
  hist->GetYaxis()->SetTitleSize(0.04); 
  hist->GetYaxis()->SetTitleOffset(1.0); 
  */
}


// ---------------------------------------------------


TLegend* createLegend(Int_t type)
{
  TLegend* leg = 0x0;
  if (type == 0)
    leg = new TLegend(0.55,0.20,0.96,0.41);
  else if (type == 1)
    leg = new TLegend(0.49,0.63,0.83,0.84);
  else
    leg = new TLegend(0.55,0.61,0.96,0.83);
  
  leg->SetTextSize(0.05);
  leg->SetFillColor(0);
  leg->SetBorderSize(1);
  
  return leg;
}

// ---------------------------------------------------


TCanvas* createCanvas(TString name)
{
  TCanvas* c = new TCanvas(name.Data(),"",760,420); 
  c->Divide(2,2);
  
  for (Int_t i = 1; i <= 4; i++) {
    c->GetPad(i)->SetLeftMargin(0.081);
    c->GetPad(i)->SetRightMargin(0.01);
    c->GetPad(i)->SetBottomMargin(0.13);
    c->GetPad(i)->SetTopMargin(0.12);
  }
  
  return c;
}

// ---------------------------------------------------

Int_t getColor(Int_t col){
  
  Int_t color = kBlack;

  switch(col)
    {
    case 0 : color = kBlue;    break;
    case 1 : color = kCyan;    break;
    case 2 : color = kGreen;   break;
    case 3 : color = kOrange;  break;
    case 4 : color = kRed;     break;
    case 5 : color = kPink;    break;
    case 6 : color = kMagenta; break;
    case 7 : color = kViolet;  break;
    case 8 : color = kBlue-2;  break;
    case 9 : color = kGreen-2;  break;
    case 10: color = kRed+2;  break;
    }

  return color;
  
}

// -------------------------------------------------------------------------
TList* getList(TString name = "jets_noGenJets_trackTypeUndef_jetTypeUndef"){

  TList* list = 0x0; 

  if(gDirectory->Get("obusch_list")){
    list = (TList*) gDirectory->Get("obusch_list");
  }
  else{
    
    TString dirName  = "PWGJE_FragmentationFunction_" + name;
    TString listName = "fracfunc_" + name;
    
    gDirectory->cd(dirName);
    list = (TList*) gDirectory->Get(listName);
        //gDirectory->cd("PWG4_FragmentationFunction_jets_noGenJets_KINEb_KINEb");
    //list = (TList*) gDirectory->Get("fracfunc_jets_noGenJets_KINEb_KINEb");
  }

  return list;
}

// ----------------------------------------------------------------------------

void sysErrorsPythiaFastJet(){
  
  gStyle->SetTitleStyle(0);
  gStyle->SetTitleX(0.5);
  gStyle->SetTitleAlign(23);
  gStyle->SetTitleY(1.0);
  gStyle->SetTitleH(0.06);//0.054);

  Double_t jetPtLim[] = {5,10,15,20,30,80}; // nBins+1 entries
  const Int_t nJetPtBins = 5;

  TString strSp[] = {"","_pi","_K","_p","_e","_mu"}; 
  const Int_t nSpecies   = 6;

  TString strTitSp[] = {"h^{+} + h^{-}","#pi^{+} + #pi^{-}","K^{+} + K^{-}","p + #bar{p}","e^{+} + e^{-}","#mu^{+} + #mu^{-}"}; 

  TString strSp_10f6a[] = {"","_piPlusPiMinus","_kPlusKMinus","_ppbar","_ePlusEMinus","_muPlusMuMinus"}; 
   
  TH1F* hSysErrEffPt[nSpecies][nJetPtBins];
  TH1F* hSysErrEffZ[nSpecies][nJetPtBins];
  TH1F* hSysErrEffXi[nSpecies][nJetPtBins];

  TH1F* hSysErrResPt[nSpecies][nJetPtBins];
  TH1F* hSysErrResZ[nSpecies][nJetPtBins];
  TH1F* hSysErrResXi[nSpecies][nJetPtBins];

  // ---

  TH1F* fh1FFTrackPtGenPrim[nSpecies][nJetPtBins];
  TH1F* fh1FFZGenPrim[nSpecies][nJetPtBins];
  TH1F* fh1FFXiGenPrim[nSpecies][nJetPtBins];
    
  TH1F* fh1FFTrackPtRecPrim_recPt_100_100[nSpecies][nJetPtBins];
  TH1F* fh1FFZRecPrim_recPt_100_100[nSpecies][nJetPtBins];
  TH1F* fh1FFXiRecPrim_recPt_100_100[nSpecies][nJetPtBins];

  TH1F* fh1FFTrackPtRecPrim_recPt_105_100[nSpecies][nJetPtBins];
  TH1F* fh1FFZRecPrim_recPt_105_100[nSpecies][nJetPtBins];
  TH1F* fh1FFXiRecPrim_recPt_105_100[nSpecies][nJetPtBins];

  TH1F* fh1FFTrackPtRecPrim_recPt_095_100[nSpecies][nJetPtBins];
  TH1F* fh1FFZRecPrim_recPt_095_100[nSpecies][nJetPtBins];
  TH1F* fh1FFXiRecPrim_recPt_095_100[nSpecies][nJetPtBins];

  TH1F* fh1FFTrackPtRecPrim_recPt_100_120[nSpecies][nJetPtBins];
  TH1F* fh1FFZRecPrim_recPt_100_120[nSpecies][nJetPtBins];
  TH1F* fh1FFXiRecPrim_recPt_100_120[nSpecies][nJetPtBins];

  TH1F* fh1FFTrackPtRecPrim_recPt_100_080[nSpecies][nJetPtBins];
  TH1F* fh1FFZRecPrim_recPt_100_080[nSpecies][nJetPtBins];
  TH1F* fh1FFXiRecPrim_recPt_100_080[nSpecies][nJetPtBins];

  TH1F* corrFacPt_10f6a[nSpecies][nJetPtBins];
  TH1F* corrFacZ_10f6a[nSpecies][nJetPtBins];
  TH1F* corrFacXi_10f6a[nSpecies][nJetPtBins];

  TH1F* corrFacPt_100_100[nSpecies][nJetPtBins];
  TH1F* corrFacZ_100_100[nSpecies][nJetPtBins];
  TH1F* corrFacXi_100_100[nSpecies][nJetPtBins];

  TH1F* corrFacPtSysEff[nSpecies][nJetPtBins];
  TH1F* corrFacZSysEff[nSpecies][nJetPtBins];
  TH1F* corrFacXiSysEff[nSpecies][nJetPtBins];

  TH1F* corrFacPtSysRes[nSpecies][nJetPtBins];
  TH1F* corrFacZSysRes[nSpecies][nJetPtBins];
  TH1F* corrFacXiSysRes[nSpecies][nJetPtBins];

  TH1F* corrFacPt_095_100[nSpecies][nJetPtBins];
  TH1F* corrFacZ_095_100[nSpecies][nJetPtBins];
  TH1F* corrFacXi_095_100[nSpecies][nJetPtBins];

  TH1F* corrFacPt_105_100[nSpecies][nJetPtBins];
  TH1F* corrFacZ_105_100[nSpecies][nJetPtBins];
  TH1F* corrFacXi_105_100[nSpecies][nJetPtBins];

  TH1F* corrFacPt_100_120[nSpecies][nJetPtBins];
  TH1F* corrFacZ_100_120[nSpecies][nJetPtBins];
  TH1F* corrFacXi_100_120[nSpecies][nJetPtBins];

  TH1F* corrFacPt_100_080[nSpecies][nJetPtBins];
  TH1F* corrFacZ_100_080[nSpecies][nJetPtBins];
  TH1F* corrFacXi_100_080[nSpecies][nJetPtBins];

  //TString strInFile10f6a = "files/outCorrections_10f6a_tpcCut.root";
  TString strInFile10f6a = "files/outCorrections_10f6a.root";
  
  TString strResDir = "files";

  TString strInFileGen(Form("%s/outCorrections_PythiaFastJet_eff100_res100.root",strResDir.Data()));
  TString strInFile_100_100(Form("%s/outCorrections_PythiaFastJet_eff100_res100.root",strResDir.Data()));
  TString strInFile_095_100(Form("%s/outCorrections_PythiaFastJet_eff095_res100.root",strResDir.Data()));
  TString strInFile_105_100(Form("%s/outCorrections_PythiaFastJet_eff105_res100.root",strResDir.Data()));
  TString strInFile_100_080(Form("%s/outCorrections_PythiaFastJet_eff100_res080.root",strResDir.Data()));
  TString strInFile_100_120(Form("%s/outCorrections_PythiaFastJet_eff100_res120.root",strResDir.Data()));
 
  // ---

  TFile f1(strInFile10f6a,"READ");
  
  //gDirectory->cd("tpcCut1");

  for(Int_t sp=0; sp<nSpecies; sp++){

    for(Int_t i=0; i<nJetPtBins; i++){
    
      TString strTitlePt(Form("hBbBCorrPt_%02d_%02d%s",(int)jetPtLim[i],(int)jetPtLim[i+1],strSp_10f6a[sp].Data()));
      TString strTitleZ(Form("hBbBCorrZ_%02d_%02d%s",(int)jetPtLim[i],(int)jetPtLim[i+1],strSp_10f6a[sp].Data()));
      TString strTitleXi(Form("hBbBCorrXi_%02d_%02d%s",(int)jetPtLim[i],(int)jetPtLim[i+1],strSp_10f6a[sp].Data()));

      corrFacPt_10f6a[sp][i] = (TH1F*) gDirectory->Get(strTitlePt);
      corrFacZ_10f6a[sp][i]  = (TH1F*) gDirectory->Get(strTitleZ);
      corrFacXi_10f6a[sp][i] = (TH1F*) gDirectory->Get(strTitleXi);

      corrFacPt_10f6a[sp][i]->SetDirectory(0);
      corrFacZ_10f6a[sp][i]->SetDirectory(0);
      corrFacXi_10f6a[sp][i]->SetDirectory(0);
    }
  }
  
  f1.Close();  

  // ---

  for(Int_t sp=0; sp<nSpecies; sp++){

    TFile f2(strInFileGen,"READ");     
    gDirectory->cd(Form("Gen%s",strSp[sp].Data()));

    for(Int_t i=0; i<nJetPtBins; i++){
    
      // for genLevel doesn't matter which histos we take - eff == 1 
      TString strTitlePtGen(Form("fh1FFTrackPtGen%s_%02d_%02d",strSp[sp].Data(),(int)jetPtLim[i],(int)jetPtLim[i+1]));
      TString strTitleZGen(Form("fh1FFZGen%s_%02d_%02d",strSp[sp].Data(),(int)jetPtLim[i],(int)jetPtLim[i+1]));
      TString strTitleXiGen(Form("fh1FFXiGen%s_%02d_%02d",strSp[sp].Data(),(int)jetPtLim[i],(int)jetPtLim[i+1]));
      
      fh1FFTrackPtGenPrim[sp][i] = (TH1F*) gDirectory->Get(strTitlePtGen);
      fh1FFZGenPrim[sp][i]       = (TH1F*) gDirectory->Get(strTitleZGen);
      fh1FFXiGenPrim[sp][i]      = (TH1F*) gDirectory->Get(strTitleXiGen);
      
      fh1FFTrackPtGenPrim[sp][i]->SetDirectory(0);
      fh1FFZGenPrim[sp][i]->SetDirectory(0);
      fh1FFXiGenPrim[sp][i]->SetDirectory(0);
    }

    f2.Close();
  }


  // ---

  for(Int_t sp=0; sp<nSpecies; sp++){
    
    TFile f3(strInFile_100_100,"READ");
    
    gDirectory->cd(Form("RecCuts%s",strSp[sp].Data()));
    
    for(Int_t i=0; i<nJetPtBins; i++){
    
      TString strTitlePtRec(Form("fh1FFTrackPtRecCuts%s_%02d_%02d",strSp[sp].Data(),(int)jetPtLim[i],(int)jetPtLim[i+1]));
      TString strTitleZRec(Form("fh1FFZRecCuts%s_%02d_%02d",strSp[sp].Data(),(int)jetPtLim[i],(int)jetPtLim[i+1]));
      TString strTitleXiRec(Form("fh1FFXiRecCuts%s_%02d_%02d",strSp[sp].Data(),(int)jetPtLim[i],(int)jetPtLim[i+1]));
      
      fh1FFTrackPtRecPrim_recPt_100_100[sp][i] = (TH1F*) gDirectory->Get(strTitlePtRec);
      fh1FFZRecPrim_recPt_100_100[sp][i]       = (TH1F*) gDirectory->Get(strTitleZRec);
      fh1FFXiRecPrim_recPt_100_100[sp][i]      = (TH1F*) gDirectory->Get(strTitleXiRec);
      
      // -- corr factors
      
      TString strTitPtCorr_100_100(Form("corrFacPt%s_100_100_%d",strSp[sp].Data(),i));
      TString strTitZCorr_100_100(Form("corrFacZ%s_100_100_%d",strSp[sp].Data(),i));
      TString strTitXiCorr_100_100(Form("corrFacXi%s_100_100_%d",strSp[sp].Data(),i));
      
      corrFacPt_100_100[sp][i] = (TH1F*) fh1FFTrackPtGenPrim[sp][i]->Clone(strTitPtCorr_100_100);
      corrFacPt_100_100[sp][i]->Divide(fh1FFTrackPtGenPrim[sp][i],fh1FFTrackPtRecPrim_recPt_100_100[sp][i],1,1,"B");
      
      corrFacZ_100_100[sp][i] = (TH1F*) fh1FFZGenPrim[sp][i]->Clone(strTitZCorr_100_100);
      corrFacZ_100_100[sp][i]->Divide(fh1FFZGenPrim[sp][i],fh1FFZRecPrim_recPt_100_100[sp][i],1,1,"B");
      
      corrFacXi_100_100[sp][i] = (TH1F*) fh1FFXiGenPrim[sp][i]->Clone(strTitXiCorr_100_100);
      corrFacXi_100_100[sp][i]->Divide(fh1FFXiGenPrim[sp][i],fh1FFXiRecPrim_recPt_100_100[sp][i],1,1,"B");
      
      // --
      
      fh1FFTrackPtRecPrim_recPt_100_100[sp][i]->SetDirectory(0);
      fh1FFZRecPrim_recPt_100_100[sp][i]->SetDirectory(0);
      fh1FFXiRecPrim_recPt_100_100[sp][i]->SetDirectory(0);
      
      corrFacPt_100_100[sp][i]->SetDirectory(0);
      corrFacZ_100_100[sp][i]->SetDirectory(0);
      corrFacXi_100_100[sp][i]->SetDirectory(0); 
    }
    
    f3.Close();
  }

  // ---

  for(Int_t sp=0; sp<nSpecies; sp++){
    
    TFile f4(strInFile_095_100,"READ");
    
    gDirectory->cd(Form("RecCuts%s",strSp[sp].Data()));
    
    for(Int_t i=0; i<nJetPtBins; i++){
      
      TString strTitlePtRec(Form("fh1FFTrackPtRecCuts%s_%02d_%02d",strSp[sp].Data(),(int)jetPtLim[i],(int)jetPtLim[i+1]));
      TString strTitleZRec(Form("fh1FFZRecCuts%s_%02d_%02d",strSp[sp].Data(),(int)jetPtLim[i],(int)jetPtLim[i+1]));
      TString strTitleXiRec(Form("fh1FFXiRecCuts%s_%02d_%02d",strSp[sp].Data(),(int)jetPtLim[i],(int)jetPtLim[i+1]));
      
      fh1FFTrackPtRecPrim_recPt_095_100[sp][i] = (TH1F*) gDirectory->Get(strTitlePtRec);
      fh1FFZRecPrim_recPt_095_100[sp][i]       = (TH1F*) gDirectory->Get(strTitleZRec);
      fh1FFXiRecPrim_recPt_095_100[sp][i]      = (TH1F*) gDirectory->Get(strTitleXiRec);
      
      // -- corr factors
      
      TString strTitPtCorr_095_100(Form("corrFacPt%s_095_100_%d",strSp[sp].Data(),i));
      TString strTitZCorr_095_100(Form("corrFacZ%s_095_100_%d",strSp[sp].Data(),i));
      TString strTitXiCorr_095_100(Form("corrFacXi%s_095_100_%d",strSp[sp].Data(),i));
      
      corrFacPt_095_100[sp][i] = (TH1F*) fh1FFTrackPtGenPrim[sp][i]->Clone(strTitPtCorr_095_100);
      corrFacPt_095_100[sp][i]->Divide(fh1FFTrackPtGenPrim[sp][i],fh1FFTrackPtRecPrim_recPt_095_100[sp][i],1,1,"B");
      
      corrFacZ_095_100[sp][i] = (TH1F*) fh1FFZGenPrim[sp][i]->Clone(strTitZCorr_095_100);
      corrFacZ_095_100[sp][i]->Divide(fh1FFZGenPrim[sp][i],fh1FFZRecPrim_recPt_095_100[sp][i],1,1,"B");
      
      corrFacXi_095_100[sp][i] = (TH1F*) fh1FFXiGenPrim[sp][i]->Clone(strTitXiCorr_095_100);
      corrFacXi_095_100[sp][i]->Divide(fh1FFXiGenPrim[sp][i],fh1FFXiRecPrim_recPt_095_100[sp][i],1,1,"B");
      
      // --
      
      fh1FFTrackPtRecPrim_recPt_095_100[sp][i]->SetDirectory(0);
      fh1FFZRecPrim_recPt_095_100[sp][i]->SetDirectory(0);
      fh1FFXiRecPrim_recPt_095_100[sp][i]->SetDirectory(0);
      
      corrFacPt_095_100[sp][i]->SetDirectory(0);
      corrFacZ_095_100[sp][i]->SetDirectory(0);
      corrFacXi_095_100[sp][i]->SetDirectory(0); 
    }
    
    f4.Close();
  }

  // ---

  for(Int_t sp=0; sp<nSpecies; sp++){
    
    TFile f5(strInFile_105_100,"READ");
    
    gDirectory->cd(Form("RecCuts%s",strSp[sp].Data()));
    
    for(Int_t i=0; i<nJetPtBins; i++){
      
      TString strTitlePtRec(Form("fh1FFTrackPtRecCuts%s_%02d_%02d",strSp[sp].Data(),(int)jetPtLim[i],(int)jetPtLim[i+1]));
      TString strTitleZRec(Form("fh1FFZRecCuts%s_%02d_%02d",strSp[sp].Data(),(int)jetPtLim[i],(int)jetPtLim[i+1]));
      TString strTitleXiRec(Form("fh1FFXiRecCuts%s_%02d_%02d",strSp[sp].Data(),(int)jetPtLim[i],(int)jetPtLim[i+1]));
      
      fh1FFTrackPtRecPrim_recPt_105_100[sp][i] = (TH1F*) gDirectory->Get(strTitlePtRec);
      fh1FFZRecPrim_recPt_105_100[sp][i]       = (TH1F*) gDirectory->Get(strTitleZRec);
      fh1FFXiRecPrim_recPt_105_100[sp][i]      = (TH1F*) gDirectory->Get(strTitleXiRec);
      
      // -- corr factors
      
      TString strTitPtCorr_105_100(Form("corrFacPt%s_105_100_%d",strSp[sp].Data(),i));
      TString strTitZCorr_105_100(Form("corrFacZ%s_105_100_%d",strSp[sp].Data(),i));
      TString strTitXiCorr_105_100(Form("corrFacXi%s_105_100_%d",strSp[sp].Data(),i));
      
      corrFacPt_105_100[sp][i] = (TH1F*) fh1FFTrackPtGenPrim[sp][i]->Clone(strTitPtCorr_105_100);
      corrFacPt_105_100[sp][i]->Divide(fh1FFTrackPtGenPrim[sp][i],fh1FFTrackPtRecPrim_recPt_105_100[sp][i],1,1,"B");
      
      corrFacZ_105_100[sp][i] = (TH1F*) fh1FFZGenPrim[sp][i]->Clone(strTitZCorr_105_100);
      corrFacZ_105_100[sp][i]->Divide(fh1FFZGenPrim[sp][i],fh1FFZRecPrim_recPt_105_100[sp][i],1,1,"B");
      
      corrFacXi_105_100[sp][i] = (TH1F*) fh1FFXiGenPrim[sp][i]->Clone(strTitXiCorr_105_100);
      corrFacXi_105_100[sp][i]->Divide(fh1FFXiGenPrim[sp][i],fh1FFXiRecPrim_recPt_105_100[sp][i],1,1,"B");
      
      // --
      
      fh1FFTrackPtRecPrim_recPt_105_100[sp][i]->SetDirectory(0);
      fh1FFZRecPrim_recPt_105_100[sp][i]->SetDirectory(0);
      fh1FFXiRecPrim_recPt_105_100[sp][i]->SetDirectory(0);
      
      corrFacPt_105_100[sp][i]->SetDirectory(0);
      corrFacZ_105_100[sp][i]->SetDirectory(0);
      corrFacXi_105_100[sp][i]->SetDirectory(0); 
    }
    
    f5.Close();
  }

  // ---

  for(Int_t sp=0; sp<nSpecies; sp++){
    
    TFile f5(strInFile_100_080,"READ");
    
    gDirectory->cd(Form("RecCuts%s",strSp[sp].Data()));
    
    for(Int_t i=0; i<nJetPtBins; i++){
      
      TString strTitlePtRec(Form("fh1FFTrackPtRecCuts%s_%02d_%02d",strSp[sp].Data(),(int)jetPtLim[i],(int)jetPtLim[i+1]));
      TString strTitleZRec(Form("fh1FFZRecCuts%s_%02d_%02d",strSp[sp].Data(),(int)jetPtLim[i],(int)jetPtLim[i+1]));
      TString strTitleXiRec(Form("fh1FFXiRecCuts%s_%02d_%02d",strSp[sp].Data(),(int)jetPtLim[i],(int)jetPtLim[i+1]));
      
      fh1FFTrackPtRecPrim_recPt_100_080[sp][i] = (TH1F*) gDirectory->Get(strTitlePtRec);
      fh1FFZRecPrim_recPt_100_080[sp][i]       = (TH1F*) gDirectory->Get(strTitleZRec);
      fh1FFXiRecPrim_recPt_100_080[sp][i]      = (TH1F*) gDirectory->Get(strTitleXiRec);
      
      // -- corr factors
      
      TString strTitPtCorr_100_080(Form("corrFacPt%s_100_080_%d",strSp[sp].Data(),i));
      TString strTitZCorr_100_080(Form("corrFacZ%s_100_080_%d",strSp[sp].Data(),i));
      TString strTitXiCorr_100_080(Form("corrFacXi%s_100_080_%d",strSp[sp].Data(),i));
      
      corrFacPt_100_080[sp][i] = (TH1F*) fh1FFTrackPtGenPrim[sp][i]->Clone(strTitPtCorr_100_080);
      corrFacPt_100_080[sp][i]->Divide(fh1FFTrackPtGenPrim[sp][i],fh1FFTrackPtRecPrim_recPt_100_080[sp][i],1,1,"B");
      
      corrFacZ_100_080[sp][i] = (TH1F*) fh1FFZGenPrim[sp][i]->Clone(strTitZCorr_100_080);
      corrFacZ_100_080[sp][i]->Divide(fh1FFZGenPrim[sp][i],fh1FFZRecPrim_recPt_100_080[sp][i],1,1,"B");
      
      corrFacXi_100_080[sp][i] = (TH1F*) fh1FFXiGenPrim[sp][i]->Clone(strTitXiCorr_100_080);
      corrFacXi_100_080[sp][i]->Divide(fh1FFXiGenPrim[sp][i],fh1FFXiRecPrim_recPt_100_080[sp][i],1,1,"B");
      
      // --
      
      fh1FFTrackPtRecPrim_recPt_100_080[sp][i]->SetDirectory(0);
      fh1FFZRecPrim_recPt_100_080[sp][i]->SetDirectory(0);
      fh1FFXiRecPrim_recPt_100_080[sp][i]->SetDirectory(0);
      
      corrFacPt_100_080[sp][i]->SetDirectory(0);
      corrFacZ_100_080[sp][i]->SetDirectory(0);
      corrFacXi_100_080[sp][i]->SetDirectory(0); 
    }
    
    f5.Close();
  }

  // ---

  for(Int_t sp=0; sp<nSpecies; sp++){
    
    TFile f6(strInFile_100_120,"READ");
    
    gDirectory->cd(Form("RecCuts%s",strSp[sp].Data()));
    
    for(Int_t i=0; i<nJetPtBins; i++){
      
      TString strTitlePtRec(Form("fh1FFTrackPtRecCuts%s_%02d_%02d",strSp[sp].Data(),(int)jetPtLim[i],(int)jetPtLim[i+1]));
      TString strTitleZRec(Form("fh1FFZRecCuts%s_%02d_%02d",strSp[sp].Data(),(int)jetPtLim[i],(int)jetPtLim[i+1]));
      TString strTitleXiRec(Form("fh1FFXiRecCuts%s_%02d_%02d",strSp[sp].Data(),(int)jetPtLim[i],(int)jetPtLim[i+1]));
      
      fh1FFTrackPtRecPrim_recPt_100_120[sp][i] = (TH1F*) gDirectory->Get(strTitlePtRec);
      fh1FFZRecPrim_recPt_100_120[sp][i]       = (TH1F*) gDirectory->Get(strTitleZRec);
      fh1FFXiRecPrim_recPt_100_120[sp][i]      = (TH1F*) gDirectory->Get(strTitleXiRec);
      
      // -- corr factors
      
      TString strTitPtCorr_100_120(Form("corrFacPt%s_100_120_%d",strSp[sp].Data(),i));
      TString strTitZCorr_100_120(Form("corrFacZ%s_100_120_%d",strSp[sp].Data(),i));
      TString strTitXiCorr_100_120(Form("corrFacXi%s_100_120_%d",strSp[sp].Data(),i));
      
      corrFacPt_100_120[sp][i] = (TH1F*) fh1FFTrackPtGenPrim[sp][i]->Clone(strTitPtCorr_100_120);
      corrFacPt_100_120[sp][i]->Divide(fh1FFTrackPtGenPrim[sp][i],fh1FFTrackPtRecPrim_recPt_100_120[sp][i],1,1,"B");
      
      corrFacZ_100_120[sp][i] = (TH1F*) fh1FFZGenPrim[sp][i]->Clone(strTitZCorr_100_120);
      corrFacZ_100_120[sp][i]->Divide(fh1FFZGenPrim[sp][i],fh1FFZRecPrim_recPt_100_120[sp][i],1,1,"B");
      
      corrFacXi_100_120[sp][i] = (TH1F*) fh1FFXiGenPrim[sp][i]->Clone(strTitXiCorr_100_120);
      corrFacXi_100_120[sp][i]->Divide(fh1FFXiGenPrim[sp][i],fh1FFXiRecPrim_recPt_100_120[sp][i],1,1,"B");
      
      // --
      
      fh1FFTrackPtRecPrim_recPt_100_120[sp][i]->SetDirectory(0);
      fh1FFZRecPrim_recPt_100_120[sp][i]->SetDirectory(0);
      fh1FFXiRecPrim_recPt_100_120[sp][i]->SetDirectory(0);
      
      corrFacPt_100_120[sp][i]->SetDirectory(0);
      corrFacZ_100_120[sp][i]->SetDirectory(0);
      corrFacXi_100_120[sp][i]->SetDirectory(0); 
    }
    
    f6.Close();
  }


  // --------------------------------------------------------

  // calc syst. error as difference between corr factors 
  // store them as bin error in histos corrFacPt(Xi)SysEff / corrFacPt(Xi)SysRes . 

  for(Int_t sp=0; sp<nSpecies; sp++){
    for(Int_t i=0; i<nJetPtBins; i++){
      
      TString strTitPtSysEff(Form("corrFacPtSysEff%s_%d",strSp[sp].Data(),i));
      TString strTitZSysEff(Form("corrFacZSysEff%s_%d",strSp[sp].Data(),i));
      TString strTitXiSysEff(Form("corrFacXiSysEff%s_%d",strSp[sp].Data(),i));
      
      TString strTitPtSysRes(Form("corrFacPtSysRes%s_%d",strSp[sp].Data(),i));
      TString strTitZSysRes(Form("corrFacZSysRes%s_%d",strSp[sp].Data(),i));
      TString strTitXiSysRes(Form("corrFacXiSysRes%s_%d",strSp[sp].Data(),i));
      
      corrFacPtSysEff[sp][i] = (TH1F*) corrFacPt_100_100[sp][i]->Clone(strTitPtSysEff);
      corrFacZSysEff[sp][i]  = (TH1F*) corrFacZ_100_100[sp][i]->Clone(strTitZSysEff);
      corrFacXiSysEff[sp][i] = (TH1F*) corrFacXi_100_100[sp][i]->Clone(strTitXiSysEff);
      
      corrFacPtSysRes[sp][i] = (TH1F*) corrFacPt_100_100[sp][i]->Clone(strTitPtSysRes);
      corrFacZSysRes[sp][i]  = (TH1F*) corrFacZ_100_100[sp][i]->Clone(strTitZSysRes);
      corrFacXiSysRes[sp][i] = (TH1F*) corrFacXi_100_100[sp][i]->Clone(strTitXiSysRes);

 
      for(Int_t bin=1; bin<=corrFacPt_100_100[sp][i]->GetNbinsX(); bin++){
      
	Double_t errPtEff = 0.5*TMath::Abs(corrFacPt_095_100[sp][i]->GetBinContent(bin) - corrFacPt_105_100[sp][i]->GetBinContent(bin));
	Double_t errPtRes = 0.5*TMath::Abs(corrFacPt_100_080[sp][i]->GetBinContent(bin) - corrFacPt_100_120[sp][i]->GetBinContent(bin));
	
	corrFacPtSysEff[sp][i]->SetBinError(bin,errPtEff);
	corrFacPtSysRes[sp][i]->SetBinError(bin,errPtRes); 
      }

      for(Int_t bin=1; bin<=corrFacZ_100_100[sp][i]->GetNbinsX(); bin++){

	Double_t errZEff = 0.5*TMath::Abs(corrFacZ_095_100[sp][i]->GetBinContent(bin) - corrFacZ_105_100[sp][i]->GetBinContent(bin));
	Double_t errZRes = 0.5*TMath::Abs(corrFacZ_100_080[sp][i]->GetBinContent(bin) - corrFacZ_100_120[sp][i]->GetBinContent(bin)); 
	
	corrFacZSysEff[sp][i]->SetBinError(bin,errZEff);
	corrFacZSysRes[sp][i]->SetBinError(bin,errZRes);
      }
      
      
      for(Int_t bin=1; bin<=corrFacXi_100_100[sp][i]->GetNbinsX(); bin++){
	
	Double_t errXiEff = 0.5*TMath::Abs(corrFacXi_095_100[sp][i]->GetBinContent(bin) - corrFacXi_105_100[sp][i]->GetBinContent(bin));
	Double_t errXiRes = 0.5*TMath::Abs(corrFacXi_100_080[sp][i]->GetBinContent(bin) - corrFacXi_100_120[sp][i]->GetBinContent(bin)); 
	
	corrFacXiSysEff[sp][i]->SetBinError(bin,errXiEff);
	corrFacXiSysRes[sp][i]->SetBinError(bin,errXiRes);
      }
    }
  }


  // --------------------------------------------------------

  gStyle->SetOptStat(0);
  
  //gStyle->SetTitleStyle(1001);      
  //gStyle->SetTitleBorderSize(1);

  Int_t selectBin = 0;
  
  const Int_t nSpeciesPlot = 4;

  // ---

  TCanvas* c1 = createCanvas("c1");
  
  TLegend* leg1 = createLegend(0);

  for(Int_t sp=0; sp<nSpeciesPlot; sp++){
    for(Int_t i=0; i<nJetPtBins-1; i++){
      
      if(i != selectBin) continue;

      c1->cd(sp+1);

      TString strPlotTitleXi(Form("%s, #it{p}_{T}^{jet, ch} = %d-%d GeV/#it{c}",strTitSp[sp].Data(),(int)jetPtLim[i],(int)jetPtLim[i+1]));

      corrFacXi_100_100[sp][i]->SetTitle(strPlotTitleXi);
      corrFacXi_100_100[sp][i]->SetXTitle("#it{#xi}");
      corrFacXi_100_100[sp][i]->SetYTitle("Correction Factor");
      setHistoStyleColor(corrFacXi_100_100[sp][i],4);
      corrFacXi_100_100[sp][i]->GetXaxis()->SetRangeUser(0,6);
      corrFacXi_100_100[sp][i]->GetYaxis()->SetRangeUser(0.61,2.2);
      corrFacXi_100_100[sp][i]->SetLineColor(4); 
      corrFacXi_100_100[sp][i]->DrawCopy();

      setHistoStyleColor(corrFacXiSysEff[sp][i],4);
      corrFacXiSysEff[sp][i]->SetFillColor(7);
      corrFacXiSysEff[sp][i]->DrawCopy("same,E2");

      setHistoStyleColor(corrFacXi_10f6a[sp][i],2);
      corrFacXi_10f6a[sp][i]->SetMarkerStyle(24);
      corrFacXi_10f6a[sp][i]->DrawCopy("same");

      if(sp==0){
	leg1->AddEntry(corrFacXi_10f6a[sp][i],"Full simulation","P");
	leg1->AddEntry(corrFacXi_100_100[sp][i],"Fast simulation","P");
	leg1->AddEntry(corrFacXiSysEff[sp][i],"Efficiency +/- 5%","F");
	leg1->Draw();
      }
      
      gPad->RedrawAxis("");
      gPad->RedrawAxis("G");
    }
  }

  // ---

  TCanvas* c2 = createCanvas("c2");
  
  TLegend* leg2 = createLegend(0);
  
  for(Int_t sp=0; sp<nSpeciesPlot; sp++){
    for(Int_t i=0; i<nJetPtBins-1; i++){
      
      if(i != selectBin) continue;
      
      c2->cd(sp+1);

      TString strPlotTitleXi(Form("%s, #it{p}_{T}^{jet, ch} = %d-%d GeV/#it{c}",strTitSp[sp].Data(),(int)jetPtLim[i],(int)jetPtLim[i+1]));

      corrFacXi_100_100[sp][i]->SetTitle(strPlotTitleXi);
      corrFacXi_100_100[sp][i]->SetXTitle("#it{#xi}");
      corrFacXi_100_100[sp][i]->SetYTitle("Correction Factor");
      setHistoStyleColor(corrFacXi_100_100[sp][i],4);
      corrFacXi_100_100[sp][i]->GetXaxis()->SetRangeUser(0,6);
      corrFacXi_100_100[sp][i]->GetYaxis()->SetRangeUser(0.61,2.2);
      corrFacXi_100_100[sp][i]->SetLineColor(4); 
      corrFacXi_100_100[sp][i]->DrawCopy();
      
      setHistoStyleColor(corrFacXiSysRes[sp][i],4);
      corrFacXiSysRes[sp][i]->SetFillColor(7);
      corrFacXiSysRes[sp][i]->DrawCopy("same,E2");
      
      setHistoStyleColor(corrFacXi_10f6a[sp][i],2);
      corrFacXi_10f6a[sp][i]->SetMarkerStyle(24);
      corrFacXi_10f6a[sp][i]->DrawCopy("same");
      
      if(sp==0){
	leg2->AddEntry(corrFacXi_10f6a[sp][i],"Full simulation","P");
	leg2->AddEntry(corrFacXi_100_100[sp][i],"Fast simulation","P");
	leg2->AddEntry(corrFacXiSysRes[sp][i],"Resolution +/- 20%","F");
	leg2->Draw();
      }
      
      gPad->RedrawAxis("");
      gPad->RedrawAxis("G");
    }
  }
  
  // --------------

  TCanvas* c3 = createCanvas("c3");

  TLegend* leg3 = createLegend(1);
  
  for(Int_t sp=0; sp<nSpeciesPlot; sp++){
    for(Int_t i=0; i<nJetPtBins-1; i++){
      
      if(i != selectBin) continue;

      c3->cd(sp+1);
      
      TString strPlotTitlePt(Form("%s, #it{p}_{T}^{jet, ch} = %d-%d GeV/#it{c}",strTitSp[sp].Data(),(int)jetPtLim[i],(int)jetPtLim[i+1]));
    
      corrFacPt_100_100[sp][i]->SetTitle(strPlotTitlePt);
      corrFacPt_100_100[sp][i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      corrFacPt_100_100[sp][i]->SetYTitle("Correction Factor");
      setHistoStyleColor(corrFacPt_100_100[sp][i],4);
      corrFacPt_100_100[sp][i]->GetXaxis()->SetRangeUser(0,jetPtLim[i+1]);
      corrFacPt_100_100[sp][i]->GetYaxis()->SetRangeUser(0.41,2.2);
      corrFacPt_100_100[sp][i]->SetLineColor(4); 
      corrFacPt_100_100[sp][i]->DrawCopy();  

      setHistoStyleColor(corrFacPtSysEff[sp][i],4);
      corrFacPtSysEff[sp][i]->SetFillColor(7);
      corrFacPtSysEff[sp][i]->DrawCopy("same,E2");
      
      setHistoStyleColor(corrFacPt_10f6a[sp][i],2);
      corrFacPt_10f6a[sp][i]->SetMarkerStyle(24);
      corrFacPt_10f6a[sp][i]->DrawCopy("same");

      if(sp==0){
	leg3->AddEntry(corrFacPt_10f6a[sp][i],"Full simulation","P");
	leg3->AddEntry(corrFacPt_100_100[sp][i],"Fast simulation","P");
	leg3->AddEntry(corrFacPtSysEff[sp][i],"Efficiency +/- 5%","F");
	leg3->Draw();
      }
      
      gPad->RedrawAxis("");
      gPad->RedrawAxis("G");
    }
  }
  
  // --------------

  TCanvas* c4 = createCanvas("c4");

  TLegend* leg4 = createLegend(1);

  for(Int_t sp=0; sp<nSpeciesPlot; sp++){
    for(Int_t i=0; i<nJetPtBins-1; i++){
      
      if(i != selectBin) continue;
      
      c4->cd(sp+1);
      
      TString strPlotTitlePt(Form("%s, #it{p}_{T}^{jet, ch} = %d-%d GeV/#it{c}",strTitSp[sp].Data(),(int)jetPtLim[i],(int)jetPtLim[i+1]));
      
      corrFacPt_100_100[sp][i]->SetTitle(strPlotTitlePt);
      corrFacPt_100_100[sp][i]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      corrFacPt_100_100[sp][i]->SetYTitle("Correction Factor");
      setHistoStyleColor(corrFacPt_100_100[sp][i],4);
      corrFacPt_100_100[sp][i]->GetXaxis()->SetRangeUser(0,jetPtLim[i+1]);
      corrFacPt_100_100[sp][i]->GetYaxis()->SetRangeUser(0.41,2.2);
      corrFacPt_100_100[sp][i]->SetLineColor(4); 
      corrFacPt_100_100[sp][i]->DrawCopy();  
      
      setHistoStyleColor(corrFacPtSysRes[sp][i],4);
      corrFacPtSysRes[sp][i]->SetFillColor(7);
      corrFacPtSysRes[sp][i]->DrawCopy("same,E2");
      
      setHistoStyleColor(corrFacPt_10f6a[sp][i],2);
      corrFacPt_10f6a[sp][i]->SetMarkerStyle(24);
      corrFacPt_10f6a[sp][i]->DrawCopy("same");
      
      if(sp==0){
	leg4->AddEntry(corrFacPt_10f6a[sp][i],"Full simulation","P");
	leg4->AddEntry(corrFacPt_100_100[sp][i],"Fast simulation","P");
	leg4->AddEntry(corrFacPtSysRes[sp][i],"Resolution +/- 20%","F");
	leg4->Draw();
      }
      
      gPad->RedrawAxis("");
      gPad->RedrawAxis("G");
    }
  }


  // ---
  
  TCanvas* c5 = createCanvas("c5");
  
  TLegend* leg5 = createLegend(2);//new TLegend(0.54,0.65,0.83,0.83);

  for(Int_t sp=0; sp<nSpeciesPlot; sp++){
    for(Int_t i=0; i<nJetPtBins-1; i++){
      
      if(i != selectBin) continue;
      
      c5->cd(sp+1);
      
      TString strPlotTitleZ(Form("%s, #it{p}_{T}^{jet, ch} = %d-%d GeV/#it{c}",strTitSp[sp].Data(),(int)jetPtLim[i],(int)jetPtLim[i+1]));

      corrFacZ_100_100[sp][i]->SetTitle(strPlotTitleZ);
      corrFacZ_100_100[sp][i]->SetXTitle("#it{z}");
      corrFacZ_100_100[sp][i]->SetYTitle("Correction Factor");
      setHistoStyleColor(corrFacZ_100_100[sp][i],4);
      corrFacZ_100_100[sp][i]->GetXaxis()->SetRangeUser(0,1.0);
      corrFacZ_100_100[sp][i]->GetYaxis()->SetRangeUser(0.35,2.2);
      if(sp>=2) corrFacZ_100_100[sp][i]->GetYaxis()->SetRangeUser(0.35,4.2);
      corrFacZ_100_100[sp][i]->SetLineColor(4); 
      corrFacZ_100_100[sp][i]->DrawCopy();
      
      setHistoStyleColor(corrFacZSysEff[sp][i],4);
      corrFacZSysEff[sp][i]->SetFillColor(7);
      corrFacZSysEff[sp][i]->DrawCopy("same,E2");
      
      setHistoStyleColor(corrFacZ_10f6a[sp][i],2);
      corrFacZ_10f6a[sp][i]->SetMarkerStyle(24);
      corrFacZ_10f6a[sp][i]->DrawCopy("same");
      
      if(sp==0){
	leg5->AddEntry(corrFacZ_10f6a[sp][i],"Full simulation","P");
	leg5->AddEntry(corrFacZ_100_100[sp][i],"Fast simulation","P");
	leg5->AddEntry(corrFacZSysEff[sp][i],"Efficiency +/- 5%","F");
	leg5->Draw();
      }
      
      gPad->RedrawAxis("");
      gPad->RedrawAxis("G");
    }
  }
    
  // ---
    
  TCanvas* c6 = createCanvas("c6");
  
  TLegend* leg6 = createLegend(2);//new TLegend(0.14,0.14,0.45,0.35);
  
  for(Int_t sp=0; sp<nSpeciesPlot; sp++){
    for(Int_t i=0; i<nJetPtBins-1; i++){
      
      if(i != selectBin) continue;
      
      c6->cd(sp+1);
      
      TString strPlotTitleZ(Form("%s, #it{p}_{T}^{jet, ch} = %d-%d GeV/#it{c}",strTitSp[sp].Data(),(int)jetPtLim[i],(int)jetPtLim[i+1]));

      corrFacZ_100_100[sp][i]->SetTitle(strPlotTitleZ);
      corrFacZ_100_100[sp][i]->SetXTitle("#it{z}");
      corrFacZ_100_100[sp][i]->SetYTitle("Correction Factor");
      setHistoStyleColor(corrFacZ_100_100[sp][i],4);
      corrFacZ_100_100[sp][i]->GetXaxis()->SetRangeUser(0,1.0);
      corrFacZ_100_100[sp][i]->GetYaxis()->SetRangeUser(0.35,2.2);
      if(sp>=2) corrFacZ_100_100[sp][i]->GetYaxis()->SetRangeUser(0.35,4.2);
      corrFacZ_100_100[sp][i]->SetLineColor(4); 
      corrFacZ_100_100[sp][i]->DrawCopy();
      
      setHistoStyleColor(corrFacZSysRes[sp][i],4);
      corrFacZSysRes[sp][i]->SetFillColor(7);
      corrFacZSysRes[sp][i]->DrawCopy("same,E2");
      
      setHistoStyleColor(corrFacZ_10f6a[sp][i],2);
      corrFacZ_10f6a[sp][i]->SetMarkerStyle(24);
      corrFacZ_10f6a[sp][i]->DrawCopy("same");
      
      if(sp==0){
	leg6->AddEntry(corrFacZ_10f6a[sp][i],"Full simulation","P");
	leg6->AddEntry(corrFacZ_100_100[sp][i],"Fast simulation","P");
	leg6->AddEntry(corrFacZSysRes[sp][i],"Resolution +/- 20%","F");
	leg6->Draw();
      }
      
      gPad->RedrawAxis("");
      gPad->RedrawAxis("G");
    }
  }

  if(savePlots){
    c1->SaveAs(Form("%scorrFacModEff_dNdXi_%02d_%02d.pdf", saveDir.Data(), (int) jetPtLim[selectBin],(int) jetPtLim[selectBin+1]));
    c2->SaveAs(Form("%scorrFacModRes_dNdXi_%02d_%02d.pdf", saveDir.Data(), (int) jetPtLim[selectBin],(int) jetPtLim[selectBin+1]));
    
    c3->SaveAs(Form("%scorrFacModEff_dNdPt_%02d_%02d.pdf", saveDir.Data(), (int) jetPtLim[selectBin],(int) jetPtLim[selectBin+1]));
    c4->SaveAs(Form("%scorrFacModRes_dNdPt_%02d_%02d.pdf", saveDir.Data(), (int) jetPtLim[selectBin],(int) jetPtLim[selectBin+1]));
    
    c5->SaveAs(Form("%scorrFacModEff_dNdz_%02d_%02d.pdf", saveDir.Data(), (int) jetPtLim[selectBin],(int) jetPtLim[selectBin+1]));
    c6->SaveAs(Form("%scorrFacModRes_dNdz_%02d_%02d.pdf", saveDir.Data(), (int) jetPtLim[selectBin],(int) jetPtLim[selectBin+1]));
  }

  // ----
  // write syst. error to file

  if(writeSysErrors){

    for(Int_t sp=0; sp<nSpecies; sp++){
      for(Int_t i=0; i<nJetPtBins; i++){ // jet slices
	for(Int_t bin=1; bin<=corrFacPtSysEff[sp][i]->GetNbinsX(); bin++){

	  Int_t jetPtLoLim = (int) jetPtLim[i];
	  Int_t jetPtUpLim = (int) jetPtLim[i+1];
	  
	  TString strNameEffPt(Form("hSysErrEffPt_%02d_%02d%s",jetPtLoLim,jetPtUpLim,strSp[sp].Data()));
	  TString strNameEffZ(Form("hSysErrEffZ_%02d_%02d%s",jetPtLoLim,jetPtUpLim,strSp[sp].Data()));
	  TString strNameEffXi(Form("hSysErrEffXi_%02d_%02d%s",jetPtLoLim,jetPtUpLim,strSp[sp].Data()));
	  
	  TString strNameResPt(Form("hSysErrResPt_%02d_%02d%s",jetPtLoLim,jetPtUpLim,strSp[sp].Data()));
	  TString strNameResZ(Form("hSysErrResZ_%02d_%02d%s",jetPtLoLim,jetPtUpLim,strSp[sp].Data()));
	  TString strNameResXi(Form("hSysErrResXi_%02d_%02d%s",jetPtLoLim,jetPtUpLim,strSp[sp].Data()));
	  
    hSysErrEffPt[sp][i] = (TH1F*) corrFacPt_100_100[sp][i]->Clone(strNameEffPt); 
    hSysErrEffZ[sp][i]  = (TH1F*) corrFacZ_100_100[sp][i]->Clone(strNameEffZ);
    hSysErrEffXi[sp][i] = (TH1F*) corrFacXi_100_100[sp][i]->Clone(strNameEffXi);
    
    hSysErrResPt[sp][i] = (TH1F*) corrFacPt_100_100[sp][i]->Clone(strNameResPt); 
    hSysErrResZ[sp][i]  = (TH1F*) corrFacZ_100_100[sp][i]->Clone(strNameResZ);
    hSysErrResXi[sp][i] = (TH1F*) corrFacXi_100_100[sp][i]->Clone(strNameResXi);
    
    /*OLD: Wrong binning!
	  hSysErrEffPt[sp][i] = (TH1F*) corrFacPt_10f6a[sp][i]->Clone(strNameEffPt); 
	  hSysErrEffZ[sp][i]  = (TH1F*) corrFacZ_10f6a[sp][i]->Clone(strNameEffZ);
	  hSysErrEffXi[sp][i] = (TH1F*) corrFacXi_10f6a[sp][i]->Clone(strNameEffXi);
	  
	  hSysErrResPt[sp][i] = (TH1F*) corrFacPt_10f6a[sp][i]->Clone(strNameResPt); 
	  hSysErrResZ[sp][i]  = (TH1F*) corrFacZ_10f6a[sp][i]->Clone(strNameResZ);
	  hSysErrResXi[sp][i] = (TH1F*) corrFacXi_10f6a[sp][i]->Clone(strNameResXi);
	  */
	  
	  hSysErrEffPt[sp][i]->Reset();
	  hSysErrEffZ[sp][i]->Reset();
	  hSysErrEffXi[sp][i]->Reset();
	  
	  hSysErrResPt[sp][i]->Reset();
	  hSysErrResZ[sp][i]->Reset();
	  hSysErrResXi[sp][i]->Reset();
	}
      }
    }
    
    
    for(Int_t sp=0; sp<nSpecies; sp++){
      for(Int_t i=0; i<nJetPtBins; i++){ // jet slices
	for(Int_t bin=1; bin<=corrFacPtSysEff[sp][i]->GetNbinsX(); bin++){
	  
 	  if(fh1FFTrackPtGenPrim[sp][i]->GetBinContent(bin) == 0) continue; 

	  Double_t sysErrEffPt = 0;
	  Double_t err  = corrFacPtSysEff[sp][i]->GetBinError(bin);
	  Double_t cont = TMath::Abs(corrFacPtSysEff[sp][i]->GetBinContent(bin));
	  
	  if(cont) sysErrEffPt = err/cont;
	  
 	  hSysErrEffPt[sp][i]->SetBinContent(bin,sysErrEffPt); 
   	  hSysErrEffPt[sp][i]->SetBinError(bin,0); 	  
	}
      }
    }
  
    // ---
  
    for(Int_t sp=0; sp<nSpecies; sp++){
      for(Int_t i=0; i<nJetPtBins; i++){ // jet slices
	for(Int_t bin=1; bin<=corrFacXiSysEff[sp][i]->GetNbinsX(); bin++){
	  
	  if(fh1FFXiGenPrim[sp][i]->GetBinContent(bin) == 0) continue; 
	  
	  Double_t sysErrEffXi = 0;
	  Double_t err  = corrFacXiSysEff[sp][i]->GetBinError(bin);
	  Double_t cont = TMath::Abs(corrFacXiSysEff[sp][i]->GetBinContent(bin));
	  
	  if(cont) sysErrEffXi = err/cont;
	  
 	  hSysErrEffXi[sp][i]->SetBinContent(bin,sysErrEffXi); 
   	  hSysErrEffXi[sp][i]->SetBinError(bin,0); 	  
	}
      }
    }    
    
    // ---

    for(Int_t sp=0; sp<nSpecies; sp++){
      for(Int_t i=0; i<nJetPtBins; i++){ // jet slices
	for(Int_t bin=1; bin<=corrFacZSysEff[sp][i]->GetNbinsX(); bin++){
	  
	  if(fh1FFZGenPrim[sp][i]->GetBinContent(bin) == 0) continue; 
	  
	  Double_t sysErrEffZ = 0;
	  Double_t err  = corrFacZSysEff[sp][i]->GetBinError(bin);
	  Double_t cont = TMath::Abs(corrFacZSysEff[sp][i]->GetBinContent(bin));
	  
	  if(cont) sysErrEffZ = err/cont;

 	  hSysErrEffZ[sp][i]->SetBinContent(bin,sysErrEffZ); 
   	  hSysErrEffZ[sp][i]->SetBinError(bin,0); 	  
	}
      }      
    }
  
    // ---

    for(Int_t sp=0; sp<nSpecies; sp++){
      for(Int_t i=0; i<nJetPtBins; i++){ // jet slices
	for(Int_t bin=1; bin<=corrFacPtSysRes[sp][i]->GetNbinsX(); bin++){

	  if(fh1FFTrackPtGenPrim[sp][i]->GetBinContent(bin) == 0) continue; 
	  
	  Double_t sysErrResPt = 0;
	  Double_t err  = corrFacPtSysRes[sp][i]->GetBinError(bin);
	  Double_t cont = TMath::Abs(corrFacPtSysRes[sp][i]->GetBinContent(bin));
	  
	  if(cont) sysErrResPt = err/cont;
   	  
	  hSysErrResPt[sp][i]->SetBinContent(bin,sysErrResPt); 
   	  hSysErrResPt[sp][i]->SetBinError(bin,0); 	  
	}
      }
    }
     
    //  ---

    for(Int_t sp=0; sp<nSpecies; sp++){
      for(Int_t i=0; i<nJetPtBins; i++){ // jet slices
	for(Int_t bin=1; bin<=corrFacXiSysRes[sp][i]->GetNbinsX(); bin++){
	  
	  if(fh1FFXiGenPrim[sp][i]->GetBinContent(bin) == 0) continue; 
	  
	  Double_t sysErrResXi = 0;
	  Double_t err  = corrFacXiSysRes[sp][i]->GetBinError(bin);
	  Double_t cont = TMath::Abs(corrFacXiSysRes[sp][i]->GetBinContent(bin));
	  
	  if(cont) sysErrResXi = err/cont;
   	
	  hSysErrResXi[sp][i]->SetBinContent(bin,sysErrResXi); 
   	  hSysErrResXi[sp][i]->SetBinError(bin,0); 	  

	}
      }
    }
    
    // ---

    for(Int_t sp=0; sp<nSpecies; sp++){
      for(Int_t i=0; i<nJetPtBins; i++){ // jet slices
	for(Int_t bin=1; bin<=corrFacZSysRes[sp][i]->GetNbinsX(); bin++){
	  
	  if(fh1FFZGenPrim[sp][i]->GetBinContent(bin) == 0) continue; 
	  
	  Double_t sysErrResZ = 0;
	  Double_t err  = corrFacZSysRes[sp][i]->GetBinError(bin);
	  Double_t cont = TMath::Abs(corrFacZSysRes[sp][i]->GetBinContent(bin));
	  
	  if(cont) sysErrResZ = err/cont;

   	  hSysErrResZ[sp][i]->SetBinContent(bin,sysErrResZ); 	  
   	  hSysErrResZ[sp][i]->SetBinError(bin,0); 	  
	}
      }
    }    
    

    TString fnameEff = Form("%soutSysErr_eff.root", saveDir.Data());
    TFile fOutEff(fnameEff,"RECREATE");
    
    cout<<" write to file "<<fnameEff<<endl;
    
    for(Int_t sp=0; sp<nSpecies; sp++){
      for(Int_t i=0; i<nJetPtBins; i++){
	
	hSysErrEffPt[sp][i]->Write();
	hSysErrEffZ[sp][i]->Write();
	hSysErrEffXi[sp][i]->Write();
      }
    }
    
    fOutEff.Close();
    
    // ---
    
    TString fnameRes = Form("%soutSysErr_res.root", saveDir.Data());
    TFile fOutRes(fnameRes,"RECREATE");
    
    cout<<" write to file "<<fnameRes<<endl;
    
    for(Int_t sp=0; sp<nSpecies; sp++){
      for(Int_t i=0; i<nJetPtBins; i++){
	
	hSysErrResPt[sp][i]->Write();
	hSysErrResZ[sp][i]->Write();
      }
    }
    
    fOutRes.Close();
  } // writeSysErrors
  


}
