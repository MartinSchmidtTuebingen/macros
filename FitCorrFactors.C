#include "TF1.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TSpline.h"
#include <iostream>

class PieceWisePoly {
public:
  PieceWisePoly(Int_t parts, Double_t* cutxvalues, Int_t* polys, Double_t* params = 0x0, Int_t smooth = 2);
  ~PieceWisePoly();
  double Eval(double *x, double* p = 0x0;);
  void SetParam(Double_t* params);
  
private:
  Int_t doSmoothing;
  Int_t nParts;
  Double_t* cuts;
  Int_t* polygrades;
  Double_t* freeparams;
};

PieceWisePoly::PieceWisePoly(Int_t parts, Double_t* cutxvalues, Int_t* polys, Double_t* params, Int_t smooth) {
  delete cuts;
  cuts = new Double_t[parts-1];
  delete polygrades;
  polygrades = new Int_t[parts];
  
  for (Int_t i=0;i<parts-1;i++) 
    cuts[i] = cutxvalues[i];
  
  for (Int_t i=0;i<parts;i++)
    polygrades[i] = polys[i];
  
  doSmoothing = smooth;
  
  SetParam(params);
}

PieceWisePoly::~PieceWisePoly() {
  delete cuts;
  cuts = 0x0;
  delete polygrades;
  polygrades = 0x0;
}

void PieceWisePoly::SetParam(Double_t* params) {
  if (!params)
    return;
  
  Int_t nOfFreeParams = 0;
  for (Int_t i=0;i<parts;i++)
    nOfFreeParams += polygrades[i];
  
  if (smooth==1)
    nOfFreeParams -= parts-1;
  else if (smooth==1)
    nOfFreeParams -= 2 * (parts-1);
  
  if (nOfFreeParams < 1) {
    std::cout << "Warning: Function does not have free parameters" << std::endl;
    return;
  }
  
  if (!freeparams) 
    freeparams = new Double_t[nOfFreeParams];
  
}

double PieceWisePoly::Eval(double* x, double* p) {
  return (2.0*x[0]);
}

void FitCorrFactors(TH1* h, Int_t particle = 0) {
  
  PieceWisePoly* poly = new PieceWisePoly();
  TF1* func = new TF1("func",poly,&PieceWisePoly::Eval,0,1,0);
  std::cout << func->Eval(2) << std::endl;
//   Int_t nOfEmptyBins = 0;
//   for (Int_t i=1;i<=h->GetNbinsX();++i) {
//     if (h->GetBinContent(i) == 1 && h->GetBinError(i) == 0) {
//       h->SetBinContent(i,0);
//       nOfEmptyBins++;
//     }
//   }
//   
//   Double_t* x = new Double_t[h->GetNbinsX() - nOfEmptyBins];
//   Double_t* y = new Double_t[h->GetNbinsX() - nOfEmptyBins];
//   
//   Int_t xnumber = 0;
//   for (Int_t i=1;i<=h->GetNbinsX();++i) {
//     if (h->GetBinContent(i) != 0 || h->GetBinError(i) != 0) {
//       x[xnumber]=h->GetBinCenter(i);
//       y[xnumber]=h->GetBinContent(i);
//       xnumber++;
//     }
//   }
//   TSpline3* spl = new TSpline3("eleff",x,y,xnumber);
//   
//   
//   //For electrons
// //   TF1* funcElectrons = new TF1("funcElectrons", "x < 0.6 ? pol6(0) : x < 3.2 ? pol3(6) : x < 20 ? pol5(9) : 0" , 0.15, 25);
// // //   TF1* funcElectrons = new TF1("funcElectrons", "x < 0.6 ? pol6(0) : x < 3.2 ? pol3(6) : x < 9 ? pol4(9) : pol1(13)" , 0.15, 50);
// //   h->Fit(funcElectrons, "R");
// //   h->GetXaxis()->SetRangeUser(0.15,50);
//   TCanvas* c = new TCanvas();
//   c->SetLogx(kTRUE);
//   spl->Draw();
//   h->Draw("same");
}

