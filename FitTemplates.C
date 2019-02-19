#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TCanvas.h"
#include "TObjArray.h"
#include "TFractionFitter.h"
#include "TFile.h"
#include "TGraph.h"
#include "TKey.h"
#include "TString.h"
#include <iostream>
using namespace std;

TObject* GetObjectOutOfTrainFile(const char* fileName, const char* objName);

void FitTemplates() {
  TH1D *data = 0x0;                              //data histogram
  TH1D *mc0 = 0x0;                               // first MC histogram
  TH1D *mc1 = 0x0;                               // second MC histogram
  TH1D *mc2 = 0x0;                               // third MC histogram
  
  TString dataFileName = "/home/nqm/mschmidt/Documents/Projects/pp13TeV/Data/DCASigmaVariation/maschmid_PID_PureGauss_PositiveNegative.root";
  const char* mcFileName = "/home/nqm/mschmidt/Documents/Projects/pp13TeV/MCData/DCASigmaVariation/maschmid_PID_PureGauss_PositiveNegative.root";
  
  TString dataHistName ="fh3DCA_XY_";
  TString mc0HistName ="fh3DCA_XY_Primaries_";
  TString mc1HistName ="fh3DCA_XY_WeakDecays_";
  TString mc2HistName ="fh3DCA_XY_Material_";
  
  TString chargeString[2] = { "Positive", "Negative" };
  
  TH2D* dataAll = 0x0;
  TFile* f = 0x0;

  TFile* output = 0x0;
  
  for (Int_t chargeNumber=0;chargeNumber<2;++chargeNumber) {
    dataHistName = "fh3DCA_XY_" + chargeString[chargeNumber];
    TString mc0HistName = TString("fh3DCA_XY_Primaries_") + chargeString[chargeNumber];
    TString mc1HistName = TString("fh3DCA_XY_WeakDecays_") + chargeString[chargeNumber];
    TString mc2HistName = TString("fh3DCA_XY_Material_") + chargeString[chargeNumber];
      
    if (dataFileName.Contains("output")) {
      f = new TFile(dataFileName.Data(), "READ");
      f->GetObject(dataHistName, dataAll);
    }
    else {
      dataAll = (TH2D*)((TH3D*)GetObjectOutOfTrainFile(dataFileName.Data(), dataHistName))->Project3D("yx");
    }
    
    if (!dataAll) {
      printf("Could not retrieve Histogram %s from file %s\n",dataHistName, dataFileName.Data());
      return;
    }
    TH2D* primariesAll = (TH2D*)((TH3D*)GetObjectOutOfTrainFile(mcFileName, mc0HistName))->Project3D("yx");
    TH2D* weakDecaysAll = (TH2D*)((TH3D*)GetObjectOutOfTrainFile(mcFileName, mc1HistName))->Project3D("yx");
    TH2D* materialAll = (TH2D*)((TH3D*)GetObjectOutOfTrainFile(mcFileName, mc2HistName))->Project3D("yx");
    const TAxis* ptaxe = dataAll->GetXaxis();
    Int_t npTBins = ptaxe->GetNbins();
    const TAxis* dcaaxe = dataAll->GetYaxis();
    Int_t nDCABins = dcaaxe->GetNbins();
    
    TH1D* hPrimariesDCA_p = new TH1D(TString::Format("hPrimariesDCA_p_%s",chargeString[chargeNumber].Data()).Data(), TString::Format("hPrimariesDCA_p_%s",chargeString[chargeNumber].Data()).Data(), npTBins, ptaxe->GetXbins()->GetArray()); 
    TH1D* hSecondariesDCA_p = new TH1D(TString::Format("hSecondariesDCA_p_%s",chargeString[chargeNumber].Data()).Data(), TString::Format("hSecondariesDCA_p_%s",chargeString[chargeNumber].Data()).Data(), npTBins, ptaxe->GetXbins()->GetArray()); 
    TH1D* hMaterialDCA_p = new TH1D(TString::Format("hMaterialDCA_p_%s",chargeString[chargeNumber].Data()).Data(), TString::Format("hMaterialDCA_p_%s",chargeString[chargeNumber].Data()).Data(), npTBins, ptaxe->GetXbins()->GetArray());
    
    TH2D* fitResult = (TH2D*)(dataAll->Clone((TString("FitResult_") + chargeString[chargeNumber]).Data()));
    fitResult->SetTitle((TString("Sum of the fitted Fractions - ")+chargeString[chargeNumber]).Data());
    for (Int_t pTBin = 0;pTBin < npTBins+1;pTBin++) {
      for (Int_t DCABin =0;DCABin < nDCABins;DCABin++) {
        fitResult->SetBinContent(pTBin,DCABin,0.0);
      }
    }
    
    Double_t value1,error1,value2,error2,value3,error3;
    
    TObjArray* outputDFiles = new TObjArray();
    
    for (Int_t pTBin = 8;pTBin < 16/*npTBins*/;pTBin++) {
      TDirectoryFile* df = new TDirectoryFile(TString::Format("Fractions_pT_%d_%d",ptaxe->GetBinLowEdge(pTBin),ptaxe->GetBinUpEdge(pTBin)).Data(),TString::Format("Fractions_pT_%d_%d",ptaxe->GetBinLowEdge(pTBin),ptaxe->GetBinUpEdge(pTBin)).Data());
      Double_t pT = (ptaxe->GetBinLowEdge(pTBin) + ptaxe->GetBinUpEdge(pTBin))/2.0;
      TObjArray *mc = new TObjArray(3);        // MC histograms are put in this array
      data = dataAll->ProjectionY("fh2DCA_XY_Proj",pTBin,pTBin);
      if (data->GetEntries() < 500)
        continue;
      
      mc0 = primariesAll->ProjectionY((TString("fh2DCA_XY_Prim_Proj_") + chargeString[chargeNumber]).Data(),pTBin,pTBin);
      mc1 = weakDecaysAll->ProjectionY((TString("fh2DCA_XY_WeakDecays_Proj_") + chargeString[chargeNumber]).Data(),pTBin,pTBin);
      mc2 = materialAll->ProjectionY((TString("fh2DCA_XY_Proj_Material_") + chargeString[chargeNumber]).Data(),pTBin,pTBin);
      mc->Add(mc0);
      mc->Add(mc1);
      mc->Add(mc2);
      TFractionFitter* fit = new TFractionFitter(data, mc); // initialise
      fit->Constrain(1,0.0,1.0);               // constrain fraction 1 to be between 0 and 1
      fit->Constrain(2,0.0,1.0);
      fit->Constrain(3,0.0,1.0);
      Int_t status = fit->Fit();               // perform the fit
      std::cout << "fit status: " << status << std::endl;
      if (status == 0) {                       // check on fit status
        TH1F* result = (TH1F*)(fit->GetPlot());
        
        fit->GetResult(0,value1,error1);
        fit->GetResult(1,value2,error2);
        fit->GetResult(2,value3,error3);
        
        hPrimariesDCA_p->SetBinContent(pTBin, value1);
        hPrimariesDCA_p->SetBinError(pTBin, error1);
        hSecondariesDCA_p->SetBinContent(pTBin, value2);
        hSecondariesDCA_p->SetBinError(pTBin, error2);
        hMaterialDCA_p->SetBinContent(pTBin, value3);
        hMaterialDCA_p->SetBinError(pTBin, error3);

        for (Int_t DCABin = 1;DCABin < nDCABins+1;DCABin++) {
          fitResult->SetBinContent(pTBin,DCABin,result->GetBinContent(DCABin));
        }
        cout << "Fractions for pT: " << pT << " in Bin " << pTBin << ": "  << value1 << " " << value2 << " " << value3 << endl;
      } 
    }
    
    output = new TFile("~/Documents/Projects/pp13TeV/DCAAnalysis/DCA_Analysis.root","UPDATE");
    
    dataAll->Write();
    fitResult->Write();
    hPrimariesDCA_p->Write();
    hSecondariesDCA_p->Write();
    hMaterialDCA_p->Write();
  }
  output->Close();
  if (f)
    f->Close();  
  return;
}

TObject* GetObjectOutOfTrainFile(const char* fileName, const char* objName) {
  TFile* file = new TFile(fileName,"READ");
  TList* keys = file->GetListOfKeys();
  TKey* key = (TKey*)keys->At(0);
  TObjArray* arr = new TObjArray();
  key->Read(arr);
  TObject* obj = arr->FindObject(objName);
  file->Close();
  return obj;
}