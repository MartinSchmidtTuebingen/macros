void DrawFromDifferentFiles() {
  const Int_t nOfFiles = 2;
//   const TString fileNames[nOfFiles] = {"~/Documents/Projects/LHC13c.pass4_DedicatedMultSplines/Result/EfficiencyCorrection_JetMCData/output_EfficiencyCorrection_outputSystematicsTotal_SummedSystematicErrors_Jets_Inclusive_PureGauss__Pt___centrality0_100__2018_04_24__2018_04_24.root"};
  const TString fileNames[nOfFiles] = {"~/Documents/Projects/LHC13c.pass4_fullJets/AnalysisResults_charged_10GeVjetCut.root", "~/Documents/Projects/LHC13c.pass4_fullJets/AnalysisResults_full_10GeVjetCut.root"};
  const Int_t nOfHistos = 2;
  const TString histNames[nOfHistos] = {"fh1nRecJetsCuts", "fh1nRecJetsCuts"};
  const TString legendNames[nOfHistos] = {"Charged Jets", "Full Jets"};
  const TString drawOptions[nOfHistos] = {"", ""};
  const Int_t drawColors[nOfHistos] = {kRed, kGreen};
  const Bool_t axeLogs[2] = {kTRUE, kTRUE};
  const TString axeTitles[2] = {"","Jet yield/Total number of Jets"};
  const TString title = "";
  TH1* histos[nOfHistos];
  const Int_t histCuts[nOfFiles] = {0,1};
  gROOT->LoadMacro("GetObjectOutOfDirectory.C");
  Double_t min = 0.0;
  Double_t max = 0.0;
  Double_t xmin = 0.0;
  Double_t xmax = 1.0;
  for (Int_t i=0;i<nOfFiles;i++) {
    TFile* f = new TFile(fileNames[i].Data()); 
    if (!f)
      continue;
    for (Int_t j=histCuts[i];j<(i == nOfFiles - 1 ? nOfHistos : histCuts[i+1]);j++) {
      TH1* h = (TH1*)GetObjectOutOfDirectory(f,histNames[j]);
      if (h) {
        histos[j] = dynamic_cast<TH1*>(h->Clone(histNames[j]+("_cloned")));
        if (histos[j]) {
          histos[j]->SetDirectory(0x0);
  //         min = TMath::Min(min, histos[j]->GetBinContent(histos[j]->GetMinimumBin()));
  //         max = TMath::Max(max, histos[j]->GetBinContent(histos[j]->GetMaximumBin()));
  //         xmin = TMath::Min(xmin, histos[j]->GetXaxis()->GetXmax());
  //         xmax = TMath::Max(xmax, histos[j]->GetXaxis()->GetXmin());
        }
      }
//       if (!histos[j]) {
//         TGraph* gr = (TGraph*)GetObjectOutOfDirectory(f,histNames[j]);
//         if (gr)
//           histos[j]->dynamic_cast<TH1*>(gr->GetHistogram()->Clone(histNames[j]+("_cloned")));          
//       }
      if (!histos[j])
        printf("No Histogram or Graph %s in File %s not found\n",fileNames[i].Data(), histNames[j].Data());
    }
    f->Close();
  }
  TCanvas* c = new TCanvas();
  c->SetLogx(axeLogs[0]);
  c->SetLogy(axeLogs[1]);  
  Bool_t namesSetted = kFALSE;
  for (Int_t j=0;j<nOfHistos;j++) {
    if (!namesSetted && histos[j]) {
      histos[j]->SetTitle(title.Data());
      if (!axeTitles[0].IsNull())
        histos[j]->GetXaxis()->SetTitle(axeTitles[0].Data());
      if (!axeTitles[1].IsNull())
        histos[j]->GetYaxis()->SetTitle(axeTitles[1].Data());
//       histos[j]->GetXaxis()->SetRangeUser(xmin,xmax);
//       histos[j]->GetYaxis()->SetRangeUser(min,max);
      namesSetted = kTRUE;
    }
    if (histos[j]) {
      if (drawColors[j] != kWhite) {
        histos[j]->SetLineColor(drawColors[j]);
        histos[j]->SetMarkerColor(drawColors[j]);
      }
      histos[j]->Scale(1.0/histos[j]->GetEntries());
      TString drawOption = drawOptions[j] + TString(" ") + TString(j==0 ? "" : "same");
      histos[j]->Draw(drawOption.Data());
    }
  }
  TLegend* l = new TLegend(0.1,0.7,0.48,0.9);
  for (Int_t j=0;j<nOfHistos;j++) {
    if (histos[j])
      l->AddEntry(histos[j], legendNames[j].Data());
  }
  l->Draw("same");
  c->SaveAs("~/Documents/Projects/LHC13c.pass4_fullJets/JetComparison_10GeVjetCut.root");
}