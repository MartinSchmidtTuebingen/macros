void changeCanvas(TCanvas* canv) {
  if (!canv)
    return;
  
  Double_t xmin = 0.0;
  Double_t xmax = 2.0;
  Double_t ymin = 0.0;
  Double_t ymax = 1.5;
  
  Float_t leftmargin = 0.09;
  Float_t rightmargin = 0.02;
  Float_t bottommargin = 0.09;
  Float_t topmargin = 0.02;
  
  UInt_t windowwidth = 800;
  UInt_t windowheight = 660;
  
  //Legend Position in percentage of the whole Canvas for jT Fractions
//   Double_t legendx1 = 0.7;    //Left 
//   Double_t legendx2 = 0.95;   //Right
//   Double_t legendy1 = 0.75;    //Bottom
//   Double_t legendy2 = 0.95;    //Top
  
  //Ratios
  Double_t legendx1 = 0.1;    //Left 
  Double_t legendx2 = 0.4;   //Right
  Double_t legendy1 = 0.65;    //Bottom
  Double_t legendy2 = 0.95;    //Top 
  
    //R fractions
//   Double_t legendx1 = 0.2;    //Left 
//   Double_t legendx2 = 0.5;   //Right
//   Double_t legendy1 = 0.35;    //Bottom
//   Double_t legendy2 = 0.65;    //Top 
  
  TList* prim = canv->GetListOfPrimitives();
  TH1* h = (TH1*)(prim->At(0));
  h->GetXaxis()->SetRangeUser(xmin,xmax);
  h->GetYaxis()->SetRangeUser(ymin,ymax);
  h->GetYaxis()->SetTitleOffset(1.1);
  h->GetYaxis()->SetTitle("Corrected Ratio");
  h->SetTitle("");
  canv->Draw();
  h->Draw("same axis");
  TLegend* legend = (TLegend*)prim->FindObject("TPave");
  
  if (legend) {
    canv->Draw();
    legend->SetX1NDC(legendx1);
    legend->SetX2NDC(legendx2);
    legend->SetY1NDC(legendy1);
    legend->SetY2NDC(legendy2);
//     legend->SetNColumns(1);
    legend->AddEntry((TH1*)(prim->At(4)),"K/#pi (MC)");
    legend->AddEntry((TH1*)(prim->At(5)),"p/K (MC)");    
    // For fractions
//     legend->AddEntry((TH1*)(prim->At(6)),"#pi (MC)");
//     legend->AddEntry((TH1*)(prim->At(7)),"K (MC)");
//     legend->AddEntry((TH1*)(prim->At(8)),"p (MC)");
  }
  
  canv->Draw();

  canv->SetLeftMargin(leftmargin);
  canv->SetRightMargin(rightmargin);
  canv->SetBottomMargin(bottommargin);
  canv->SetTopMargin(topmargin);  
  
  canv->Draw();
  
  canv->SetWindowSize(windowwidth, windowheight); 
}