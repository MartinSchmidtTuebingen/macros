#include "TF1.h"
#include "TString.h"
#include <cmath>

#include "PieceWisePoly.h"


PieceWisePoly::PieceWisePoly(Int_t parts, Double_t* cutxvalues, Int_t* polys, Double_t xmin, Double_t xmax,  Double_t* params, Int_t smooth)
  : doSmoothing(2)
  , nParts(parts)
  , cuts(0x0)
  , nOfFreeParams(0)
  , polyParameters(0x0)
  , piecewisepolynom(0x0)
  {
  nParts = parts;
  if (parts > 1) {
    cuts = new Double_t[nParts-1];
    
    for (Int_t i=0;i<nParts-1;i++) 
      cuts[i] = cutxvalues[i];
  }
  
  polyParameters = new Int_t[nParts];
  for (Int_t i=0;i<nParts;++i)
    polyParameters[i]=polys[i];
  
  doSmoothing = smooth;
  
  TString functionString = "";
  
  Int_t remParts = nParts;
  Int_t previousNOfParameters = 0;
  while (remParts > 1) {
    functionString += TString::Format("x < %f ? pol%d(%d) : ",cuts[nParts-remParts],polyParameters[nParts-remParts]-1,previousNOfParameters);
    previousNOfParameters += polyParameters[nParts-remParts];
    remParts--;
  }
  functionString += TString::Format("pol%d(%d)",polyParameters[nParts-remParts]-1,previousNOfParameters);
  piecewisepolynom = new TF1("piecewisepolynom",functionString.Data(),xmin,xmax);
  nOfFreeParams = previousNOfParameters += polyParameters[nParts-remParts];
  
  if (doSmoothing==1)
    nOfFreeParams -= nParts-1;
  else if (doSmoothing==2)
    nOfFreeParams -= 2 * (nParts-1);
  
  SetParam(params);
}

PieceWisePoly::~PieceWisePoly() {
  delete cuts;
  cuts = 0x0;
  delete piecewisepolynom;
  piecewisepolynom = 0x0;
}

void PieceWisePoly::SetParam(Double_t* params) {
  if (!params)
    return;
  
  Bool_t changed = kFALSE;
  
  for (Int_t j=0;j<polyParameters[0];++j) {
    if (changed) {
      piecewisepolynom->SetParameter(j,params[j]);
    }
    else if (piecewisepolynom->GetParameter(j) != params[j]) {
        changed = kTRUE;
        piecewisepolynom->SetParameter(j,params[j]);
    }
  }

  Int_t parNumber = polyParameters[0];
  Int_t freeparNumber = polyParameters[0];
  
  for (Int_t i=1;i<nParts;++i) {
    for (Int_t j=doSmoothing;j<polyParameters[i];++j) {
      Int_t internalparNumber = parNumber + j;
      if (changed) {
        piecewisepolynom->SetParameter(internalparNumber,params[freeparNumber]);
      }
      else if (piecewisepolynom->GetParameter(internalparNumber) != params[freeparNumber]) {
          changed = kTRUE;
          piecewisepolynom->SetParameter(internalparNumber,params[freeparNumber]);
      }
      freeparNumber++;
    }
    if (changed) {
      for (Int_t j=doSmoothing-1;j>=0;--j) 
        piecewisepolynom->SetParameter(parNumber+j,SumUp(cuts[i-1],parNumber-polyParameters[i-1],parNumber-1,j,j)-SumUp(cuts[i-1],parNumber,parNumber+polyParameters[i]-1,j,j+1));
    }
    
    parNumber += polyParameters[i];
  }
  return;
}

Double_t PieceWisePoly::SumUp(Double_t constant, Int_t start, Int_t end, Int_t derivative, Int_t startn) {
  Double_t sum = 0.0;
  if (derivative == 1) {
    for (Int_t i=startn;i<=end-start;++i) {
      Double_t potence = 1.0;
      for (Int_t j=1;j<i;++j)
        potence *= constant;
      
      sum += i * piecewisepolynom->GetParameter(start+i) * potence;
    }
  }
  else {
    for (Int_t i=startn;i<=end-start;++i) {
      Double_t potence = 1.0;
      for (Int_t j=1;j<=i;++j)
        potence *= constant;
      
      sum += piecewisepolynom->GetParameter(start+i) * potence;
    }
  }
  return sum;
}

TF1* PieceWisePoly::GetPartFunction(Int_t i) {
  if (i>nParts)
    return 0x0;
  
  TString functionString = TString::Format("pol%d(0)",polyParameters[i]-1);
  Int_t previousNOfParameters = 0;
  for (Int_t j=0;j<i;++j) 
    previousNOfParameters += polyParameters[j];
  
  TF1* func = new TF1(TString::Format("func_%d",i).Data(),functionString.Data(),0,1);
  
  for (Int_t j=0;j<polyParameters[i];++j)
    func->SetParameter(j,piecewisepolynom->GetParameter(previousNOfParameters+j));
  
  return func;
}


double PieceWisePoly::Eval (double x, double* p) {
  SetParam(p);
  return piecewisepolynom->Eval(x);
}
  