//Write out the results from fast MC simulation, not calculating correction factors, scaled with nOfJets

#include "TFile.h"
#include "AliCFContainer.h"
#include "AliPID.h"
#include "TH2.h"

#include "THnSparseDefinitions.h"

#include <iostream>
using namespace std;

Int_t iPt     = 0;
Int_t iMCid   = 0;
Int_t iEta    = 0;
Int_t iCharge = 0;
Int_t iMult   = 0;
Int_t iJetPt  = 0;
Int_t iZ = 0;
Int_t iXi = 0;
Int_t iDistance = 0;
Int_t ijT = 0;

Int_t iObsAxis = 0;

Int_t writeOutCorrectionFiles(TString effFile, TString outfilepath, TString addoutFileName) {
  TFile* fileEff = new TFile(effFile.Data());
  if (!fileEff) {
    printf("Failed to open efficiency file \"%s\"\n", effFile.Data());
    return -1;
  }
  
  AliCFContainer* data = (AliCFContainer*)(fileEff->Get("containerEff"));
  if (!data) {
    printf("Failed to load efficiency container!\n");
    return -1;
  }  
  
  const Int_t nOfEffSteps = 2;
  EffSteps usedEffSteps[nOfEffSteps] = {kStepGenWithGenCuts, kStepRecWithRecCutsPrimaries};
  TString dirNameEffSteps[nOfEffSteps] = {"Gen", "RecCuts"};
  
  // For backward compatibility:
  // Check whether "P_{T}" or "p_{T}" is used
  TString momentumString = "p";
  for (Int_t i = 0; i < data->GetNVar(); i++) {
    TString temp = data->GetVarTitle(i);
    if (temp.Contains("P_{")) {
      momentumString = "P";
      break;
    }
    else if (temp.Contains("p_{")) {
      momentumString = "p";
      break;
    }
  }
  
//   AliPID::ParticleShortName(species)
  
  iPt     = data->GetVar(Form("%s_{T} (GeV/c)", momentumString.Data()));
  iMCid   = data->GetVar("MC ID");
  iEta    = data->GetVar("#eta");
  iCharge = data->GetVar("Charge (e_{0})");
  iMult   = data->GetVar("Centrality Percentile");
  iJetPt = data->GetVar(Form("%s_{T}^{jet} (GeV/c)", momentumString.Data()));
  iZ     = data->GetVar(Form("z = %s_{T}^{track} / %s_{T}^{jet}", momentumString.Data(), momentumString.Data()));
  iXi    = data->GetVar(Form("#xi = ln(%s_{T}^{jet} / %s_{T}^{track})", momentumString.Data(), momentumString.Data()));
  iDistance = data->GetVar("R");
  ijT = data->GetVar("j_{T} (GeV/c)");  
  
  const Int_t nOfTrackObservables = 5;
  Int_t trackObservableBins[nOfTrackObservables] = {iPt, iZ, iXi, iDistance, ijT};
  TString observableNames[nOfTrackObservables] = {"TrackPt", "Z", "Xi", "R", "jT"};
  
  //Setting jet limits and getting the number of rec/gen jets from the file associated with the efficiency file
  const Int_t nOfJetBins = 5;
  Double_t jetPtLimits[2*nOfJetBins] = {5.0,10.0,10.0,15.0,15.0,20.0,20.0,30.0,30.0,80.0};
  Int_t nOfJets[2][nOfJetBins] = {0};
  Int_t jetBinLimits[2*nOfJetBins] = {0};
  
  TString pathNameDataMC = effFile;
  pathNameDataMC.ReplaceAll("_efficiency", "");
  
  TFile* fDataMC = TFile::Open(pathNameDataMC.Data());
  if (!fDataMC)  {
    std::cout << std::endl;
    std::cout << "Failed to open file \"" << pathNameDataMC.Data() << "\" to obtain num of rec/gen jets!" << std::endl;   
    return -1;
  }
  
  TString listName = pathNameDataMC;
  listName.Replace(0, listName.Last('/') + 1, "");
  listName.ReplaceAll(".root", "");
    
  TObjArray* histList = (TObjArray*)(fDataMC->Get(listName.Data()));
  if (!histList) {
    std::cout << std::endl;
    std::cout << "Failed to load list \"" << listName.Data() << "\" to obtain num of rec/gen jets!" << std::endl;
    return -1;
  }    
  
  TH2* hNjetsGen = (TH2D*)histList->FindObject("fh2FFJetPtGen");
  TH2* hNjetsRec = (TH2D*)histList->FindObject("fh2FFJetPtRec");
  
  Int_t lowerCentralityBinLimit = -1;
  Int_t upperCentralityBinLimit = -1;

  for (Int_t i=0;i<nOfJetBins;++i) {
    jetBinLimits[2*i] = data->GetAxis(iJetPt, 0)->FindFixBin(jetPtLimits[2*i] + 0.001);
    jetBinLimits[2*i+1] = data->GetAxis(iJetPt, 0)->FindFixBin(jetPtLimits[2*i+1] - 0.001);   
    
    nOfJets[0][i] = hNjetsGen ? hNjetsGen->Integral(lowerCentralityBinLimit, upperCentralityBinLimit, jetBinLimits[2*i], jetBinLimits[2*i+1]) : 1.;
    nOfJets[1][i] = hNjetsRec ? hNjetsRec->Integral(lowerCentralityBinLimit, upperCentralityBinLimit, jetBinLimits[2*i], jetBinLimits[2*i+1]) : 1.;
  }
  
  TFile* outFile = new TFile((outfilepath + TString("/outCorrections_PythiaFastJet_") + addoutFileName + TString(".root")).Data(),"RECREATE");
  
  for (Int_t species=-1;species<AliPID::kSPECIES;++species) {
    TString speciesString = species >= 0 ? (TString("_") + TString(AliPID::ParticleShortName(species))) : TString("");
    data->SetRangeUser(iMCid,species+1,species+1,kTRUE);
    for (Int_t effStep=0;effStep<nOfEffSteps;++effStep) {
      TString dirName = dirNameEffSteps[effStep] + speciesString;
      outFile->mkdir(dirName.Data());
      outFile->cd(dirName.Data());
      for (Int_t jetPtStep = 0;jetPtStep<nOfJetBins;++jetPtStep) {
        data->SetRangeUser(iJetPt,jetBinLimits[2*jetPtStep],jetBinLimits[2*jetPtStep+1],kTRUE);
        for (Int_t observable = 0;observable<nOfTrackObservables;++observable) {
          TH1* h = data->Project(usedEffSteps[effStep],trackObservableBins[observable]);
          h->SetNameTitle(TString::Format("fh1FF%s%s%s_%02d_%02d",observableNames[observable].Data(),dirNameEffSteps[effStep].Data(),speciesString.Data(),(Int_t)jetPtLimits[jetPtStep*2],(Int_t)jetPtLimits[jetPtStep*2+1]),"");
    
          if (nOfJets[TMath::Min(effStep,1)][jetPtStep])
            h->Scale(1.0/nOfJets[TMath::Min(effStep,1)][jetPtStep]);
          
          for (Int_t binNumber = 0;binNumber<=h->GetNbinsX();binNumber++) 
            h->SetBinContent(binNumber,h->GetBinContent(binNumber)/h->GetBinWidth(binNumber));
            
          h->Write();
        }
      }
    }
  }
  fileEff->Close();
  outFile->Close();
  return 0;
}