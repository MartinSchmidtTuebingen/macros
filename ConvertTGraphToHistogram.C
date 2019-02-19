TH1F* ConvertTGraphToHistogram(TGraph* gr) {
  
  Int_t nbins = gr->GetN();
  
  Double_t* xbins = new Double_t[nbins+1];
  
  for (Int_t i=0;i<nbins;++i)
    xbins[i] = gr->GetX()[i] - gr->GetErrorXlow(i);
  
  xbins[nbins] = gr->GetX()[nbins-1] + gr->GetErrorXhigh(nbins-1);
  
  TH1F* h = gr->GetHistogram()->Clone("cloned");
  h->SetBins(nbins,xbins);
  
  h->SetContent(gr->GetY());
  
  for (Int_t i=1;i<=nbins;++i) {
//     h->SetBinContent(i,gr->GetY()[i-1]);
    h->SetBinError(i,gr->GetErrorY(i-1));
  }
    
  
  return h;
  
}