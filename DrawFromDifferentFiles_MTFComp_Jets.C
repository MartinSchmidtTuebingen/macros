void DrawFromDifferentFiles() {
  const Int_t nOfFiles = 5;
  const TString fileNames[nOfFiles] = { "~/Documents/Promotion/Analyses/LHC13c.pass4_DedicatedMultSplines/Result/EfficiencyCorrection_JetMCData/output_EfficiencyCorrection_outputSystematicsTotal_SummedSystematicErrors_Jets_Inclusive_PureGauss__Pt___centrality0_10__2018_04_24__2018_04_24.root", 
    "~/Documents/Promotion/Analyses/LHC13c.pass4_DedicatedMultSplines/Result/EfficiencyCorrection_JetMCData/output_EfficiencyCorrection_outputSystematicsTotal_SummedSystematicErrors_Jets_PureGauss__Pt___centrality0_10_jetPt5.0_10.0__2018_04_24__2018_04_24.root", "~/Documents/Promotion/Analyses/LHC13c.pass4_DedicatedMultSplines/Result/EfficiencyCorrection_JetMCData/output_EfficiencyCorrection_outputSystematicsTotal_SummedSystematicErrors_Jets_PureGauss__Pt___centrality0_10_jetPt10.0_15.0__2018_04_24__2018_04_24.root", "~/Documents/Promotion/Analyses/LHC13c.pass4_DedicatedMultSplines/Result/EfficiencyCorrection_JetMCData/output_EfficiencyCorrection_outputSystematicsTotal_SummedSystematicErrors_Jets_PureGauss__Pt___centrality0_10_jetPt15.0_20.0__2018_04_24__2018_04_24.root", "~/Documents/Promotion/Analyses/LHC13c.pass4_DedicatedMultSplines/Result/EfficiencyCorrection_JetMCData/output_EfficiencyCorrection_outputSystematicsTotal_SummedSystematicErrors_Jets_PureGauss__Pt___centrality0_10_jetPt20.0_30.0__2018_04_24__2018_04_24.root"};
  TString outputfile = "~/Documents/Promotion/Analyses/LHC13c.pass4_DedicatedMultSplines/Comp_Jet_Inclusive_Proton_Cent0_10.root";
  const Int_t nOfHistos = 15;
  const Int_t nDiffHistNames = 3;
  const TString histNames[nDiffHistNames] = {"hRatioToPiProtons_corrected", "hRatioToPiProtons_sysError_corrected", "hMCtoPi_p"};
  const TString legendNamesInclusive[nDiffHistNames] = {"p/#pi, Inclusive", "", "p/#pi, Inclusive (MC)"};
  const TString legendNames[nDiffHistNames] = {"p/#pi, p_{T}^{jet} %d - %d GeV/#it{c}", "", "p/#pi, p_{T}^{jet} %d - %d GeV/#it{c} (MC)"};
  const TString drawOptions[nDiffHistNames] = {"", "E2", ""};
  Int_t jetPt[5] = {5,10,15,20,30};
  const Int_t drawColors[5] = {kGreen, kRed, kBlue, kMagenta, kBlack};
  const Bool_t axeLogs[2] = {kTRUE, kFALSE};
  const TString axeTitles[2] = {"",""};
  const TString title = "";
  TH1* histos[nOfHistos];
  for (Int_t j=0;j<nOfHistos;j++) {
    histos[j] = 0x0;
  }
  const Int_t histCuts[nOfFiles] = {0,3,6,9,12};
  gROOT->LoadMacro("GetObjectOutOfDirectory.C");
  Double_t min = 0.0;
  Double_t max = 0.0;
  Double_t xmin = 0.0;
  Double_t xmax = 1.0;
  for (Int_t i=0;i<nOfFiles;i++) {
    TFile* f = new TFile(fileNames[i].Data());
    if (f->IsZombie()) {
      continue;
    }
    for (Int_t j=histCuts[i];j<(i < nOfFiles - 1 ? histCuts[i+1] : nOfHistos);j++) {
      Int_t nameNumber = j - (TMath::FloorNint(j/nDiffHistNames) * nDiffHistNames);
//       Int_t nameNumber = j;
      TH1* h = (TH1*)GetObjectOutOfDirectory(f,histNames[nameNumber]);
      if (h) {
        histos[j] = dynamic_cast<TH1*>(h->Clone(histNames[nameNumber]+TString::Format("_cloned_%d",jetPt[i])));
        if (histos[j]) {
          histos[j]->SetDirectory(0x0);
        }
      }
      if (!histos[j]) {
        TGraph* gr = (TGraph*)GetObjectOutOfDirectory(f,histNames[nameNumber]);
        if (gr)
          histos[j] = dynamic_cast<TH1*>(gr->GetHistogram()->Clone(histNames[nameNumber]+("_cloned")));          
      }
      if (!histos[j])
        printf("No Histogram or Graph %s found in File %s /n", histNames[nameNumber].Data(),fileNames[i].Data());
    }
    f->Close();
  }
  TCanvas* c = new TCanvas("c2","c2");
  c->SetLogx(axeLogs[0]);
  c->SetLogy(axeLogs[1]);  
  Bool_t namesSetted = kFALSE;
  TLegend* l = new TLegend(0.1,0.7,0.48,0.9);
  for (Int_t j=0;j<nOfHistos;j++) {
    if (!histos[j])
      continue;
    Int_t nameNumber = j - (TMath::FloorNint(j/nDiffHistNames) * nDiffHistNames);
    Int_t fileNumber = TMath::FloorNint(j/nDiffHistNames);
    if (!namesSetted) {
      histos[j]->SetTitle(title.Data());
      if (!axeTitles[0].IsNull())
        histos[j]->GetXaxis()->SetTitle(axeTitles[0].Data());
      if (!axeTitles[1].IsNull())
        histos[j]->GetYaxis()->SetTitle(axeTitles[1].Data());
      namesSetted = kTRUE;
    }
    histos[j]->SetFillStyle(0);
    histos[j]->SetLineColor(drawColors[TMath::FloorNint(j/nDiffHistNames)]);
    histos[j]->SetMarkerColor(drawColors[TMath::FloorNint(j/nDiffHistNames)]);
    histos[j]->SetFillColor(drawColors[TMath::FloorNint(j/nDiffHistNames)]);
    TString drawOption = drawOptions[nameNumber] + TString(" ") + TString(j==0 ? "" : "same");
    histos[j]->Draw(drawOption.Data());
    if (!legendNames[nameNumber].IsNull()) {
      if (fileNumber == 0)
        l->AddEntry(histos[j], legendNamesInclusive[nameNumber].Data(),"lp");
      else {
        l->AddEntry(histos[j], TString::Format(legendNames[nameNumber].Data(),jetPt[fileNumber - 1],jetPt[fileNumber]).Data(),"lp");
      }
    }
  }
  l->Draw("same");
  c->SaveAs(outputfile.Data());
}