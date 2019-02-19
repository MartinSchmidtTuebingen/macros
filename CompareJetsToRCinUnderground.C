void CompareJetsToRCinUnderground(TString filepath, TString filename) {
  
  TFile* f = new TFile((filepath+filename).Data(),"UPDATE");
  gROOT->LoadMacro("/home/nqm/mschmidt/Documents/macros/GetObjectOutOfDirectory.C");
  TH1D* jets = (TH1D*)GetObjectOutOfDirectory(f, "fh1nRecJetsCuts");
  if (!jets) {
    cout << "fh1nRecJetsCuts not found"
    return;
  }
  Int_t n = jets->FindFirstBinAbove(0.0);
  TH1D* jetsNormalized = (TH1D*)(jets->Clone());
  jetsNormalized->SetNameTitle("fh1nRecJetsCuts_normalized", "fh1nRecJetsCuts_normalized");
  jetsNormalized->Scale(1.0/jetsNormalized->Integral(n, 240));
  
  TH1D* RCinUE = (TH1D*)GetObjectOutOfDirectory(f, "fh1nRCinUnderground");
  if (!RCinUE) {
    cout << "fh1nRCinUnderground not found"
    return;
  }  
  TH1D* RCinUE_normalized = (TH1D*)(RCinUE->Clone());
  RCinUE_normalized->SetNameTitle("fh1nRCinUnderground_normalized", "fh1nRCinUnderground_normalized");
  RCinUE_normalized->Scale(1.0/RCinUE_normalized->Integral(n, 240));  
  
  TH1D* relation = (TH1D*)(jetsNormalized->Clone());
  relation->Divide(RCinUE_normalized);
  relation->SetNameTitle("fh1Relation", "fh1Relation");

  jetsNormalized->Write(jetsNormalized->GetName, TObject::kOverwrite);
  RCinUE_normalized->Write(RCinUE_normalized->GetName, TObject::kOverwrite);
  relation->Write(relation->GetName, TObject::kOverwrite);
  f->Close();
}