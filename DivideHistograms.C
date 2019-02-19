TH1* DivideHistograms(TH1* h1, TH1* h2, TAxis* axe = 0x0) {
  
  
  TAxis* ax1 = h1->GetXaxis();
  TAxis* ax2 = h2->GetXaxis();
  
  axe = ax1;
  
  Double_t* xbins = new Double_t[axe->GetNbins()+1];
  
  for (Int_t i=0;i<axe->GetNbins();++i)
    xbins[i]=axe->GetBinLowEdge(i+1);
  
  xbins[axe->GetNbins()]=axe->GetBinUpEdge(axe->GetNbins());
  
  TH1F* h1divided = (TH1F*)h1->Clone("h1divided");
  h1divided->SetBins(axe->GetNbins(),xbins);
  
  Double_t min1 = ax1->GetXmin();
  Double_t min2 = ax2->GetXmin();
  Double_t max1 = ax1->GetXmax();
  Double_t max2 = ax2->GetXmax();
  
  for (Int_t n=1;n<=h1divided->GetNbinsX();++n) {
    Double_t binlowEdge = axe->GetBinLowEdge(n);
    Double_t binUpEdge = axe->GetBinUpEdge(n);
    if (binlowEdge < min1 || binlowEdge < min2)
      continue;
    if (binUpEdge > max1 || binUpEdge > max2)
      break;
    Int_t bin1low = h1->FindFixBin(binlowEdge);
    Int_t bin1Up = h1->FindFixBin(binUpEdge);
    Int_t bin2low = h2->FindFixBin(binlowEdge);
    Int_t bin2Up = h2->FindFixBin(binUpEdge);
    Double_t content1 = 0.0;
    if (bin1low == bin1Up)
      content1 += h1->GetBinContent(bin1low) * (binUpEdge - binlowEdge) * (binUpEdge - binlowEdge)/h1->GetBinWidth(bin1low);
    else {
      content1 += h1->GetBinContent(bin1low) * (h1->GetXaxis()->GetBinUpEdge(bin1low) - binlowEdge) * (h1->GetXaxis()->GetBinUpEdge(bin1low) - binlowEdge)/h1->GetBinWidth(bin1low);
      for (Int_t i=bin1low+1;i<bin1Up;++i) {
        content1 += h1->GetBinContent(i) * h1->GetBinWidth(i);
      }
      content1 += h1->GetBinContent(bin1Up) * (binUpEdge - h1->GetBinLowEdge(bin1Up)) * (binUpEdge - h1->GetBinLowEdge(bin1Up))/h1->GetBinWidth(bin1Up);
    }
    Double_t content2 = 0.0;
    if (bin2low == bin2Up)
      content2 += h2->GetBinContent(bin2low) * (binUpEdge - binlowEdge) * (binUpEdge - binlowEdge)/h2->GetBinWidth(bin2low);
    else {
      content2 += h2->GetBinContent(bin2low) * (h2->GetXaxis()->GetBinUpEdge(bin2low) - binlowEdge) * (h2->GetXaxis()->GetBinUpEdge(bin2low) - binlowEdge)/h2->GetBinWidth(bin2low);
      for (Int_t i=bin2low+1;i<bin2Up;++i) {
        content2 += h2->GetBinContent(i) * h2->GetBinWidth(i);
      }
      content2 += h2->GetBinContent(bin2Up) * (binUpEdge - h2->GetBinLowEdge(bin2Up)) * (binUpEdge - h2->GetBinLowEdge(bin2Up))/h2->GetBinWidth(bin2Up);
    }
    if (TMath::Abs(content2) < 1e-40)
      h1divided->SetBinContent(n,0.0);
    else
      h1divided->SetBinContent(n,content1/content2);
  }
  
  return h1divided;
}