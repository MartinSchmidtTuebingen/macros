Double_t smear (Double_t truevalue);

void Unfolding_Examples() {
  gSystem->Load("libRooUnfold.so");
  
  RooUnfoldResponse response (40, -10.0, 10.0);
  TH2F* resp = new TH2F("response", "response", 40, -10, 10, 40, -10, 10);
  TH1D* hTrue = new TH1D("true","TRUTH",40,-10.0,10.0);
  hTrue->SetLineColor(kGreen);
  TH1D* hMeas = new TH1D("meas","Measured",40, -10.0, 10.0);
  hMeas->SetLineColor(kBlue);
  
  for (Int_t i=0;i<100000;i++) {
    Double_t truevalue = gRandom->BreitWigner(0.3,2.5);
    Double_t meas = smear(truevalue);
    if (meas!=-99999.0) {
      response.Fill(meas,truevalue);
      resp->Fill(meas,truevalue);
    }
    else
      response.Miss(truevalue);
  }
    
  for (Int_t i=0;i<10000;++i) {
    Double_t truevalue = gRandom->BreitWigner(0.3,2.5);
    hTrue->Fill(truevalue);
    Double_t meas = smear(truevalue);
    if (meas != -99999.0)
      hMeas->Fill(meas);
  }
  
//   TCanvas* c1 = new TCanvas("c1","Response Matrix");
//   resp->Draw("colz");
//   
//   TCanvas* c2 = new TCanvas("c2","Distributions");
//   c2->SetLogy(kTRUE);
//   hTrue->Draw();
//   hMeas->Draw("same");

  RooUnfoldBayes unfoldBayes(&response,hMeas,4);
  RooUnfoldSvd unfoldSVD(&response,hMeas,35);
  
  TH1D* hRecoBayes = (TH1D*)unfoldBayes.Hreco();
  hRecoBayes->SetLineColor(kRed);
  hRecoBayes->SetLineStyle(7);
  TH1D* hRecoSVD = (TH1D*)unfoldSVD.Hreco();
  hRecoSVD->SetLineColor(kRed);
  hRecoSVD->SetLineStyle(7);
  
  TCanvas* c1 = new TCanvas("c1","SVD");
  c1->SetLogy(kTRUE);
  hMeas->Draw();
  hTrue->Draw("same");
  hRecoSVD->Draw("same");
  
  TCanvas* c2 = new TCanvas("c2","Bayes");
  c2->SetLogy(kTRUE);
  hMeas->Draw();
  hTrue->Draw("same");
  hRecoBayes->Draw("same");
}


Double_t smear (Double_t truevalue)
{
  // f u n c t i o n to smear the true value
  // to a m e a s u r e d value with limited e f f i c i e n c y
  // first we cook up some e f f i c i e n c y ...

  Double_t xeff = 0.3 + (1.0 -0.3) /20*( truevalue +10.0) ;
  Double_t x = gRandom->Rndm();
  // random number
  if (x > xeff ) return -99999.0; // s i m u l a t e e f f i c i e n c y

  // then s i m u l a t e a g a u s s i a n r e s o l u t i o n effect
  Double_t resolution = gRandom->Gaus( -2.5 ,0.2) ;
  return truevalue + resolution;
}