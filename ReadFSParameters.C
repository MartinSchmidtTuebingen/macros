#include "TFile.h"
#include "TString.h"
#include "TNamed.h"
#include "TObjArray.h"
#include "TF1.h"
#include "TObjString.h"

#include "AliPID.h"

#include "./PieceWisePoly.cxx"

#include <iostream>
using namespace std;

void ReadFSParameters(TString parameterFile) {
  TFile* f = new TFile(parameterFile.Data(),"UPDATE");
  
  for (Int_t species=0;species<AliPID::kSPECIES;++species) {
    for (Int_t charge=0;charge<=1;++charge) {
      TString name = TString::Format("fastSimulationParameters_%s_%d",AliPID::ParticleShortName(species),charge);
      TNamed* cont = (TNamed*)f->FindObjectAny(name.Data());
      if (!name)
        cout << "Could not find " << name << "." << endl;
      
      TString* parString = new TString(cont->GetTitle());
      TObjArray* arrPar = parString->Tokenize(";");
      Int_t nOfParts = (((TObjString*)(arrPar->At(0)))->GetString()).Atoi();
      Double_t* cuts = new Double_t[nOfParts-1];
      Int_t* nparameters = new Int_t[nOfParts]; 
      for (Int_t part = 0;part<nOfParts - 1;++part) {
        nparameters[part] = (((TObjString*)(arrPar->At(nOfParts + part)))->GetString()).Atoi();
        cuts[part] = (((TObjString*)(arrPar->At(1 + part)))->GetString()).Atof();
      }
      nparameters[nOfParts - 1] = (((TObjString*)(arrPar->At(2*nOfParts -1)))->GetString()).Atoi();
      
      PieceWisePoly* pwp = new PieceWisePoly(nOfParts,cuts,nparameters,0,50,0x0,2);
      TString nameFunction = TString::Format("fastSimulationFunction_%s_%s",AliPID::ParticleShortName(species),charge ? "neg" : "pos");
      TF1* func = new TF1(nameFunction.Data(),pwp,0,50,pwp->GetNOfParam());
      for (Int_t param=0;param<pwp->GetNOfParam();++param) {
        func->SetParameter(param,(((TObjString*)(arrPar->At(2*nOfParts + param)))->GetString()).Atof());
        cout << func->GetParameter(param) << endl;
      }
      func->Write();
      delete nparameters;
      delete cuts;
      delete arrPar;
      delete parString;
    }
  }
}