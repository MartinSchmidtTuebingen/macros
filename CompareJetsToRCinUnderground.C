void CompareJetsToRCinUnderground(TString filename) {
  
  TObject* saveDir;
  
  TFile* f = new TFile(filename.Data(),"UPDATE");
  gROOT->LoadMacro("/home/martin/Documents/macros/GetObjectOutOfDirectory.C");
  TH1D* jets = (TH1D*)GetObjectOutOfDirectory(f, "fh1nRecJetsCuts", saveDir);
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
  
  cout << saveDir << endl;
  /*
  if (saveDir->InheritsFrom("TCollection")) {
    TCollection* saveCollection = (TCollection*)saveDir;
    saveCollection->Add(jetsNormalized);
    saveCollection->Add(RCinUE_normalized);
    saveCollection->Add(relation);
  }
  else {*/
    jetsNormalized->Write(jetsNormalized->GetName(), TObject::kOverwrite);
    RCinUE_normalized->Write(RCinUE_normalized->GetName(), TObject::kOverwrite);
    relation->Write(relation->GetName(), TObject::kOverwrite);
//   }
  f->Close();
}