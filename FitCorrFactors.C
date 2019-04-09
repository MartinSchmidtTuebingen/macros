#include "TF1.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TSpline.h"
#include "TString.h"
#include <iostream>

class PieceWisePoly {
public:
  PieceWisePoly(Int_t parts, Double_t* cutxvalues, Int_t* polys, Double_t* params = 0x0, Int_t smooth = 2);
  ~PieceWisePoly();
  double Eval(double *x, double* p = 0x0;);
  void SetParam(Double_t* params);
  
  Int_t GetNOfParam() const {return nOfFreeParams;};
  
private:
  Int_t doSmoothing;
  Int_t nParts;
  Double_t* cuts;
  Int_t nOfFreeParams;
  TF1* piecewisepolynom;
};

PieceWisePoly::PieceWisePoly(Int_t parts, Double_t* cutxvalues, Int_t* polys, Double_t* params, Int_t smooth) {
  if (parts > 1) {
    cuts = new Double_t[parts-1];
    
    for (Int_t i=0;i<parts-1;i++) 
      cuts[i] = cutxvalues[i];
  }
  
  doSmoothing = smooth;
  
  TString functionString = "";
  
  Int_t remParts = parts;
  Int_t 
  while (remParts > 1) {
    functionString += "x < " + cuts[parts-remParts] + " ? pol" + polys[parts-remParts] + "() :";
    remParts--;
  }
  functionString += " pol" + polys[parts] + "()";
  
  
  
  
  Int_t nOfFreeParams = 0;
  for (Int_t i=0;i<parts;i++)
    nOfFreeParams += polygrades[i];
  
  if (doSmoothing==1)
    nOfFreeParams -= parts-1;
  else if (doSmoothing==2)
    nOfFreeParams -= 2 * (parts-1);
  
  if (nOfFreeParams < 1) {
    std::cout << "Warning: Function does not have free parameters" << std::endl;
    return;
  }  
  
  SetParam(params);
  
  
  
}

PieceWisePoly::~PieceWisePoly() {
  delete cuts;
  cuts = 0x0;
  delete piecewisepolynom;
}

void PieceWisePoly::SetParam(Double_t* params) {
  if (!params)
    return;
  
  if (!freeparams) 
    freeparams = new Double_t[nOfFreeParams];
  
  for (Int_t i=0;i<nOfFreeParams;++i)
    freeparams[i] = params[i];
  
  return;
}

double PieceWisePoly::Eval(double* x, double* p) {
  SetParam(p);
  Double_t xx = x[0];
  Int_t polyNumber = 0;
  while (cuts[polyNumber] < x)
    polyNumber++;
  
  return polynom[polyNumber]->Eval(xx);
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

