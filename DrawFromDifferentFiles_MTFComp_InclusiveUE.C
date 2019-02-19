void DrawFromDifferentFiles() {
  const Int_t nOfFiles = 5;
  const TString fileNames[nOfFiles] = {"~/Documents/Projects/LHC13c.pass4_DedicatedMultSplines/Result/EfficiencyCorrection_JetMCData/output_EfficiencyCorrection_outputSystematicsTotal_SummedSystematicErrors_Jets_Inclusive_PureGauss__Pt___centrality60_100__2018_04_24__2018_04_24.root",
    "~/Documents/Projects/LHC13c.pass4_DedicatedMultSplines/Result/EfficiencyCorrection_JetMCData/output_EfficiencyCorrection_outputSystematics_Jets_UE_PC_PureGauss_Pt_Method_centrality60_100_jetPt5.0_10.0_nSigma0.0__2018_04_24.root", "~/Documents/Projects/LHC13c.pass4_DedicatedMultSplines/Result/EfficiencyCorrection_JetMCData/output_EfficiencyCorrection_outputSystematics_Jets_UE_PC_PureGauss_Pt_Method_centrality60_100_jetPt10.0_15.0_nSigma0.0__2018_04_24.root", "~/Documents/Projects/LHC13c.pass4_DedicatedMultSplines/Result/EfficiencyCorrection_JetMCData/output_EfficiencyCorrection_outputSystematics_Jets_UE_PC_PureGauss_Pt_Method_centrality60_100_jetPt15.0_20.0_nSigma0.0__2018_04_24.root", "~/Documents/Projects/LHC13c.pass4_DedicatedMultSplines/Result/EfficiencyCorrection_JetMCData/output_EfficiencyCorrection_outputSystematics_Jets_UE_PC_PureGauss_Pt_Method_centrality60_100_jetPt20.0_30.0_nSigma0.0__2018_04_24.root"};
  TString outputfile = "~/Documents/Projects/LHC13c.pass4_DedicatedMultSplines/Comp_UE_Inclusive_Proton_Cent60_100.root";
  const Int_t nOfHistos = 15;
  const TString histNames[3] = {"hRatioToPiKaons_corrected", "hRatioToPiKaons_sysError_corrected", "hMCtoPi_K"};
  const TString legendNamesInclusive[3] = {"K/#pi, Inclusive", "", "K/#pi (MC), Inclusive"};
  const TString legendNames[3] = {"K/#pi, p_{T}^{jet} %d - %d GeV/#it{c}", "", "K/#pi (MC), p_{T}^{jet} %d - %d GeV/#it{c}"};
  const TString drawOptions[3] = {"", "E2",""};
  Int_t jetPt[5] = {5,10,15,20,30};
  const Int_t drawColors[5] = {kGreen, kRed, kBlue, kMagenta, kBlack};
  const Bool_t axeLogs[2] = {kTRUE, kFALSE};
  const TString axeTitles[2] = {"",""};
  const TString title = "";
  TH1* histos[nOfHistos];
  for (Int_t j=0;j<nOfFiles;j++) {
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
    if (!f)
      continue;
    for (Int_t j=histCuts[i];j<(TMath::Min(histCuts[i+1],nOfHistos));j++) {
      Int_t nameNumber = j - (TMath::FloorNint(j/3) * 3);
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
        printf("No Histogram or Graph %s found in File %s /n",fileNames[i].Data(), histNames[nameNumber].Data());
    }
    f->Close();
  }
  TCanvas* c = new TCanvas();
  c->SetLogx(axeLogs[0]);
  c->SetLogy(axeLogs[1]);  
  Bool_t namesSetted = kFALSE;
  TLegend* l = new TLegend(0.1,0.7,0.48,0.9);
  for (Int_t j=0;j<9;j++) {
    Int_t nameNumber = j - (TMath::FloorNint(j/3) * 3);
    Int_t fileNumber = TMath::FloorNint(j/3);
    if (nameNumber == 1)
      continue;
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
//       if (drawColors[j] != kWhite) {
//         histos[j]->SetLineColor(drawColors[j]);
//         histos[j]->SetMarkerColor(drawColors[j]);
        histos[j]->SetFillStyle(0);
//         histos[j]->SetLineColor(histos[j]->GetLineColor()+TMath::FloorNint(j/6)*2);
//         histos[j]->SetMarkerColor(histos[j]->GetLineColor());
//         histos[j]->SetFillColor(histos[j]->GetLineColor());
        histos[j]->SetLineColor(drawColors[TMath::FloorNint(j/3)]);
        histos[j]->SetMarkerColor(drawColors[TMath::FloorNint(j/3)]);
        histos[j]->SetFillColor(drawColors[TMath::FloorNint(j/3)]);
//       }
//       histos[j]->Scale(1.0/histos[j]->GetEntries());
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
  }
  l->Draw("same");
  c->SaveAs(outputfile.Data());
}