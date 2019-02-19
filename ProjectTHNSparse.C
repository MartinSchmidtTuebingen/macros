void ProjectTHNSparse(TString filename) {
  TFile* f = new TFile(filename.Data(),"READ");
  THnSparse* h = (THnSparse*)f->FindObjectAny("TPCPIDEtaQA")->FindObject("hPIDdataAll");
  //Identification
  h->GetAxis(0)->SetRange(4,4);
  //Species 
  h->GetAxis(1)->SetRange(2,2);
  //pT
  h->GetAxis(2)->SetRangeUser(1.0,2.0);
  //Multiplicity
//   h->GetAxis(3)->SetRangeUser(0,500);
  //Eta
  h->GetAxis(5)->SetRangeUser(-0.2,0.2);
  TCanvas* c = new TCanvas();
  
  const Int_t numberOfProjections = 4;
  const Int_t multiplicities[2*numberOfProjections] = {0,500,500,1000,1000,1500, 1500, 2000};
  TH1* hProj[numberOfProjections] = {0x0,0x0,0x0,0x0};
  
  TString name;
  
  for (Int_t i=0;i<numberOfProjections;++i) {
    h->GetAxis(3)->SetRangeUser(multiplicities[2*i],multiplicities[2*i+1]);
    hProj[i] = h->Projection(4);
    name = TString::Format("Multiplicity%i_%i",multiplicities[2*i],multiplicities[2*i+1]);
    hProj[i]->SetNameTitle(name.Data(),name.Data());
    hProj[i]->SetLineColor(i+1);
    hProj[i]->Scale(1.0/hProj[i]->GetEntries());
    cout << hProj[i]->GetSkewness() << " " << hProj[i]->GetKurtosis() << endl;
    hProj[i]->Fit("gaus");
    hProj[i]->GetFunction("gaus")->SetLineColor(i+1);
    hProj[i]->GetFunction("gaus")->SetLineWidth(1);
  }
  for (Int_t i=0;i<numberOfProjections;++i) {
    hProj[i]->Draw(i==0 ? "" : "same");
  }
  c->BuildLegend();
}