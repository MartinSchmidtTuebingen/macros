void DrawFromDifferentFiles() {
  const Int_t nOfFiles = 2;
  const TString fileNames[nOfFiles] = {"~/Documents/Promotion/Analyses/LHC13c.pass4_DedicatedMultSplines/Result/EfficiencyCorrection_JetMCData/output_EfficiencyCorrection_outputSystematicsTotal_SummedSystematicErrors_Jets_Inclusive_PureGauss__Pt___centrality60_100__2018_04_24__2018_04_24.root", "~/Documents/HEPData-ins1244523-v1-Table_8.root"};
  TString outputfile = "~/Documents/Promotion/Analyses/LHC13c.pass4_DedicatedMultSplines/ComparisonInclusive_YieldProtons_60_100_MTF_To_ins1244523.root";
  const Int_t nOfHistos = 4;
  const Int_t nOfOutHistos = 2;
  TString canvasName = "Protons_60_100";
  const TString histNames[nOfHistos] = {"hYieldProtons_corrected", "hYieldProtons_sysError_corrected", "Graph1D_y2", "Graph1D_y3"};
  const TString legendNames[nOfHistos+nOfOutHistos] = {"MTF 60 - 100 (V0A)", "", "Reference 60 - 80 (V0A)", "Reference 80 - 100 (V0A)", "MTF/Reference 60 - 80 (V0A)", "MTF/Reference 80 - 100 (V0A)"};
//   const TString legendNames[nOfOutHistos] = {"MTF Result", "", "Reference 0 - 5 (V0A)", "Reference 5 - 10 (V0A)"};
//   const TString drawOptions[nOfOutHistos] = {"", "E2", "", ""};
  const TString drawOptions[nOfHistos+nOfOutHistos] = {"", "E2", "", "", "", ""};
//   const Int_t drawColors[nOfOutHistos] = {kBlue, kBlue, kBlue+1, kBlue+2};
//   const Int_t drawColors[nOfOutHistos] = {kBlue, kBlue, kBlue+1, kBlue+2};
//   const Int_t drawColors[nOfOutHistos] = {kBlue, kBlue+3};
  const Int_t drawColors[nOfHistos+nOfOutHistos] = {kBlue, kBlue, kBlue+1, kBlue+2,kBlue+1, kBlue+2};
//   const Int_t drawColors[nOfOutHistos] = {kBlue, kBlue+3};
  const Bool_t axeLogs[2] = {kTRUE, kTRUE};
  const TString axeTitles[2] = {"",""};
  const TString title = "";
  TH1* histos[nOfHistos] = {0x0,};
  TH1* histosOut[nOfOutHistos] = {0x0,};
  const Int_t histCuts[nOfFiles] = {0,2};
  gROOT->LoadMacro("GetObjectOutOfDirectory.C");
  Double_t min = 0.0;
  Double_t max = 0.0;
  Double_t xmin = 0.0;
  Double_t xmax = 1.0;
  
  gROOT->LoadMacro("ConvertTGraphToHistogram.C");
  for (Int_t i=0;i<nOfFiles;i++) {
    TFile* f = new TFile(fileNames[i].Data()); 
    if (!f)
      continue;
    for (Int_t j=histCuts[i];j<(i == nOfFiles - 1 ? nOfHistos : histCuts[i+1]);j++) {
      TH1* h =  dynamic_cast<TH1*>GetObjectOutOfDirectory(f,histNames[j]);
      if (h && TString(h->ClassName()).Contains("TH")) {
        histos[j] = dynamic_cast<TH1*>(h->Clone(histNames[j]+("_cloned")));
//         if (histos[j]) {
  //         min = TMath::Min(min, histos[j]->GetBinContent(histos[j]->GetMinimumBin()));
  //         max = TMath::Max(max, histos[j]->GetBinContent(histos[j]->GetMaximumBin()));
  //         xmin = TMath::Min(xmin, histos[j]->GetXaxis()->GetXmax());
  //         xmax = TMath::Max(xmax, histos[j]->GetXaxis()->GetXmin());
//         }
      }
      if (!histos[j]) {
        TGraph* gr = (TGraph*)GetObjectOutOfDirectory(f,histNames[j]);
        if (gr) {
          histos[j] = ConvertTGraphToHistogram(gr);
//           histos[j] = SetNameTitle(histNames[j]+("_cloned"),histNames[j]+("_cloned"));
//           histos[j] = SetNameTitle("test");
        }
      }
      if (!histos[j])
        printf("No Histogram or Graph %s found in File %\n",fileNames[i].Data(), histNames[j].Data());
      else
        histos[j]->SetDirectory(0x0);
    }
    f->Close();
  }
  
  gROOT->LoadMacro("DivideHistograms.C");
  histosOut[0] = DivideHistograms(histos[0],histos[2]);
  histosOut[0]->SetNameTitle("h1Ratio_0_5", "h1Ratio_0_5");
    
  histosOut[1] = DivideHistograms(histos[0],histos[3]);
  histosOut[1]->SetNameTitle("h1Ratio_5_10", "h1Ratio_5_10");
  
//   for (Int_t i=0;i<4;++i)
//     histosOut[i]=histos[i];

  TCanvas* c = new TCanvas(canvasName.Data(), canvasName.Data());
  c->Divide(1,2);
  c->cd(1);
  c->cd(1)->SetLogx(axeLogs[0]);
  c->cd(1)->SetLogy(axeLogs[1]);  
  Bool_t namesSetted = kFALSE;
  for (Int_t j=0;j<nOfHistos;j++) {
    if (!namesSetted && histos[j]) {
      histos[j]->SetTitle(title.Data());
      if (!axeTitles[0].IsNull())
        histos[j]->GetXaxis()->SetTitle(axeTitles[0].Data());
      if (!axeTitles[1].IsNull())
        histos[j]->GetYaxis()->SetTitle(axeTitles[1].Data());
      histos[j]->GetXaxis()->SetRangeUser(0.1,50);
      namesSetted = kTRUE;
    }
    if (histos[j]) {
      if (drawColors[j] != kWhite) {
        histos[j]->SetLineColor(drawColors[j]);
        histos[j]->SetMarkerColor(drawColors[j]);
        histos[j]->SetFillColor(drawColors[j]);
      }
      TString drawOption = drawOptions[j] + TString(" ") + TString(j==0 ? "" : "same");
      histos[j]->Draw(drawOption.Data());
    }
  }
  TLegend* l = new TLegend(0.1,0.7,0.48,0.9);
  for (Int_t j=0;j<nOfHistos;j++) {
    if (histos[j] && !legendNames[j].IsNull())
      l->AddEntry(histos[j], legendNames[j].Data());
  }
  l->Draw("same");
  
  c->cd(2);
  c->cd(2)->SetLogx(axeLogs[0]);
  c->cd(2)->SetLogy(kFALSE);  
  namesSetted = kFALSE;
  for (Int_t j=0;j<nOfOutHistos;j++) {
    if (!namesSetted && histosOut[j]) {
      histosOut[j]->SetTitle(title.Data());
      histosOut[j]->GetYaxis()->SetTitle("Ratio");
      histosOut[j]->GetXaxis()->SetRangeUser(0.1,50);
      histosOut[j]->GetYaxis()->SetRangeUser(0.01,2.2);
      histosOut[j]->GetYaxis()->SetLabelSize(0.05);
      histosOut[j]->GetYaxis()->SetTitleSize(0.05);
      histosOut[j]->GetYaxis()->SetTitleOffset(0.75);
      namesSetted = kTRUE;
    }
    if (histosOut[j]) {
      if (drawColors[j] != kWhite) {
        histosOut[j]->SetLineColor(drawColors[j+4]);
        histosOut[j]->SetMarkerColor(drawColors[j+4]);
        histosOut[j]->SetFillColor(drawColors[j+4]);
      }
      TString drawOption = drawOptions[j+4] + TString(" ") + TString(j==0 ? "" : "same");
      histosOut[j]->Draw(drawOption.Data());
    }
  }
  l1 = new TLegend(0.1,0.7,0.48,0.9);
  for (Int_t j=0;j<nOfOutHistos;j++) {
    if (histosOut[j] && !legendNames[j+4].IsNull())
      l1->AddEntry(histosOut[j], legendNames[j+4].Data());
  }
  l1->Draw("same");
  
  c->SaveAs(outputfile.Data());
}