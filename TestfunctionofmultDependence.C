void FitAndSetPointsInGraph(TH1D* h, TGraphErrors* gSlvsInvOffset, TGraphErrors* gCurvvsInvOffset, Int_t pTBin, TGraphErrors* gSlvsTanTheta = 0x0, TGraphErrors* gCurvvsTanTheta = 0x0, Int_t ThetaBin = -1);

const Int_t numTanThetaBins = 10;
  
const  Double_t tanThetaBinEdges[numTanThetaBins * 2] = {
    -1.0, -0.8,
    -0.8, -0.6,
    -0.6, -0.4,
    -0.4, -0.2,
    -0.2, 0.0,
    0.0, 0.2,
    0.2, 0.4,
    0.4, 0.6,
    0.6, 0.8,
    0.8, 1.0 
};

void TestfunctionofmultDependence() {
  gROOT->LoadMacro("./GetObjectOutOfDirectory.C");
  gROOT->LoadMacro("./GetObjectOutOfCanvas.C+");
  
  TFile* f = new TFile("~/Documents/Promotion/SplineCreation/LHC15o_pass2_lowIR_MultDifferentiated/bhess_PIDetaTree_multiplicityDependence_2018_09_30__19_10.root","READ");
  TFile* saveFile = new TFile("~/Documents/Promotion/SplineCreation/LHC15o_pass2_lowIR_MultDifferentiated/testpol2multDependence_Simple.root","RECREATE");  
  
//   TFile* f = new TFile("~/Documents/SplineCreation/LHC15o_pass2_MultDifferentiated/bhess_PIDetaTree_multiplicityDependence_2018_09_30__19_10.root","READ");
//   TFile* saveFile = new TFile("~/Documents/SplineCreation/LHC15o_pass2_MultDifferentiated/testpol2multDependence.root","RECREATE");
  
  saveFile->mkdir("allTanTheta");
  saveFile->cd("allTanTheta");
  
  TGraphErrors* gSlvsInvOffset = new TGraphErrors(24);
  gSlvsInvOffset->SetName("slopes2_allEta");
  TGraphErrors* gCurvvsInvOffset = new TGraphErrors(24);
  gCurvvsInvOffset->SetName("curves2_allEta");
  
  TH1D* hSlVarvsTanTheta = new TH1D("hSlVarvsTanTheta", "Slope of slopes vs TanTheta",10,-1.0,1.0};
  TH1D* hCurvVarvsTanTheta = new TH1D("hCurvVarvsTanTheta", "Slope of curvature vs TanTheta",10,-1.0,1.0};
  
  TH1D* hSlopeAllEta = new TH1D("hSlopeAllEta", "Slope of slopes for all Eta",1,-1.0,1.0);
  TH1D* hCurvAllEta = new TH1D("hCurvAllEta", "Slope of curvature for all Eta",1,-1.0,1.0);
  
  for (Int_t i=0;i<=23;++i) {
    gSlvsInvOffset->SetPoint(i,-10,0);
    gCurvvsInvOffset->SetPoint(i,-10,0);
    
    TCanvas* c = (TCanvas*)GetObjectOutOfDirectory(f,TString::Format("canv_pBin%d_allEta",i));
    TH1D* h = (TH1D*)GetObjectOutOfCanvas(c,TString::Format("hMean_allEta_pBin%d",i));
    h->SetDirectory(0x0); 
    FitAndSetPointsInGraph(h,gSlvsInvOffset, gCurvvsInvOffset, i);
    h->Write();
    delete h;
    h = 0x0;
  }
  TFitResultPtr fitRes = gSlvsInvOffset->Fit("pol1","s","same",0,10000);
  if (((Int_t)fitRes) == 0) {
    Double_t p1 = fitRes->GetParams()[1];
    Double_t err1 = fitRes->GetErrors()[1];   
    hSlopeAllEta->SetBinContent(1,p1);
    hSlopeAllEta->SetBinError(1,err1);   
  } 
  TF1* func = new TF1("pol1withparameters","pol1");
  Double_t firstx, firsty, lastx,lasty;
  gCurvvsInvOffset->GetPoint(0,firstx,firsty);
  gCurvvsInvOffset->GetPoint(23,lastx,lasty);
//   func->SetParameter(0,firsty);
//   func->SetParameter(1,(lasty-firsty)/(lastx-firstx));
  TFitResultPtr fitRes = gCurvvsInvOffset->Fit(func,"s","same",0,10000);
  Double_t CurvAllEtaaxis = 0.0;
  if (((Int_t)fitRes) == 0) {
    CurvAllEtaaxis = fitRes->GetParams()[0];
    Double_t p1 = fitRes->GetParams()[1];
    Double_t err1 = fitRes->GetErrors()[1];   
    hCurvAllEta->SetBinContent(1,p1);
    hCurvAllEta->SetBinError(1,err1);
  }
  hCurvAllEta->Write();
  hSlopeAllEta->Write();
  gCurvvsInvOffset->Write();
  gSlvsInvOffset->Write();
  delete gSlvsInvOffset;
  delete gCurvvsInvOffset;
  gSlvsInvOffset = 0x0;
  gCurvvsInvOffset = 0x0;
  
//   TGraphErrors* gSlvsTanTheta[24];
//   TGraphErrors* gCurvvsTanTheta[24];
//   
//   for (Int_t i=0;i<24;i++) {
//     gSlvsTanTheta[i] =  new TGraphErrors(numTanThetaBins);
//     gSlvsTanTheta[i]->SetName(Form("tanThetaSlopes_pBin%d", i));
//     gCurvvsTanTheta[i] =  new TGraphErrors(numTanThetaBins);
//     gCurvvsTanTheta[i]->SetName(Form("tanThetaCurves_pBin%d", i));
//   }
//   
//   for (Int_t j=0;j<=9;j++) {
//     saveFile->mkdir(TString::Format("tanThetaBin%d",j).Data());
//     saveFile->cd(TString::Format("tanThetaBin%d",j).Data());    
//     gSlvsInvOffset = new TGraphErrors(24);
//     gSlvsInvOffset->SetName(Form("slopes2_tanTheta%d", j));
//     gCurvvsInvOffset = new TGraphErrors(24);
//     gCurvvsInvOffset->SetName(Form("curves2_tanTheta%d", j));
//     for (Int_t i=0;i<=23;++i) {
//       gSlvsInvOffset->SetPoint(i,-10,0);
//       gCurvvsInvOffset->SetPoint(i,-10,0);
//       
//       TCanvas* c = (TCanvas*)GetObjectOutOfDirectory(f,TString::Format("canv_pBin%d_tanThetaBin%d",i,j));
//       TH1D* h = (TH1D*)GetObjectOutOfCanvas(c,TString::Format("hMean_pBin%d_tanThetaBin_%d",i,j));
//       h->SetDirectory(0x0);
//       FitAndSetPointsInGraph(h,gSlvsInvOffset, gCurvvsInvOffset, i, gSlvsTanTheta[i], gCurvvsTanTheta[i], j);
//       h->Write();
//       delete h;
//       h = 0x0;      
//     }
//     TFitResultPtr fitRes = gSlvsInvOffset->Fit("pol1","s","same",0,10000);
//     if (((Int_t)fitRes) == 0) {
//       Double_t p1 = fitRes->GetParams()[1];
//       Double_t err1 = fitRes->GetErrors()[1];   
//       hSlVarvsTanTheta->SetBinContent(j+1,p1);
//       hSlVarvsTanTheta->SetBinError(j+1,err1);   
//     } 
//     TF1* fcurvvsInvOffset = new TF1("fcurvvsInvOffset", "pol1");
//     fcurvvsInvOffset->SetParameter(0,CurvAllEtaaxis);
//     fcurvvsInvOffset->SetParameter(1,hCurvAllEta->GetBinContent(1));
//     TFitResultPtr fitRes = gCurvvsInvOffset->Fit(fcurvvsInvOffset,"s","same",0,10000);
//     if (((Int_t)fitRes) == 0) {
//       Double_t p1 = fitRes->GetParams()[1];
//       Double_t err1 = fitRes->GetErrors()[1];   
//       hCurvVarvsTanTheta->SetBinContent(j+1,p1);
//       hCurvVarvsTanTheta->SetBinError(j+1,err1);
//     }    
//     gSlvsInvOffset->Write();
//     gCurvvsInvOffset->Write();
//     delete gSlvsInvOffset;
//     delete gCurvvsInvOffset;
//     gSlvsInvOffset = 0x0;
//     gCurvvsInvOffset = 0x0;    
//   }
//   
//   saveFile->mkdir("tanThetaFits");
//   saveFile->cd("tanThetaFits");
//   for (Int_t i=0;i<24;i++) {
//     gCurvvsTanTheta[i]->Write();
//     gSlvsTanTheta[i]->Write();    
//   }
//   hSlVarvsTanTheta->Write();
//   hCurvVarvsTanTheta->Write();
  f->Close();
  saveFile->Close();
}

void FitAndSetPointsInGraph(TH1D* h, TGraphErrors* gSlvsInvOffset, TGraphErrors* gCurvvsInvOffset, Int_t pTBin, TGraphErrors* gSlvsTanTheta, TGraphErrors* gCurvvsTanTheta, Int_t ThetaBin) {
  TFitResultPtr fitRes = h->Fit("pol2","s","",0,10000);
  TF1* stepfunc = new TF1("stepfunc","(x < (-[1]/(2*[2])) ? pol2(0) : [0] - [1]*[1]/(4*[2]))");
  stepfunc->SetParameter(0,fitRes->GetParams()[0]);
  stepfunc->SetParameter(1,fitRes->GetParams()[1]);
  stepfunc->SetParameter(2,fitRes->GetParams()[2]);
  TFitResultPtr fitRes = h->Fit(stepfunc,"s","same");
  if (((Int_t)fitRes) == 0) {
    Double_t p0 = fitRes->GetParams()[0];
    Double_t p1 = fitRes->GetParams()[1];
    Double_t p2 = fitRes->GetParams()[2];
    Double_t err0 = fitRes->GetErrors()[0];
    Double_t err1 = fitRes->GetErrors()[1];   
    Double_t err2 = fitRes->GetErrors()[2];
    
    Double_t totErrorSlope = 99999;
    if (p0 != 0 && p1 != 0) {
      totErrorSlope = TMath::Abs(p1 / p0 * TMath::Sqrt(TMath::Power(err0 / p0, 2) + TMath::Power(err1 / p1, 2)));
    }
    Double_t totErrorCurv = 99999;
    if (p0 != 0 && p2 != 0) {
      totErrorCurv = TMath::Abs(p2 / p0 * TMath::Sqrt(TMath::Power(err0 / p0, 2) + TMath::Power(err2 / p2, 2)));
    }
//     if (totErrorSlope/(TMath::Abs(p1/p0)) <= 1) {
      gSlvsInvOffset->SetPoint(pTBin,1/p0,p1/p0);
      gSlvsInvOffset->SetPointError(pTBin,TMath::Abs(err0/p0/p0),totErrorSlope);
//     }
//     if (totErrorCurv/(TMath::Abs(p2/p0)) <= 1) {
      gCurvvsInvOffset->SetPoint(pTBin,1/p0,p2/p0);
      gCurvvsInvOffset->SetPointError(pTBin,TMath::Abs(err0/p0/p0),totErrorCurv);
//     }
    
    if (gSlvsTanTheta) {    
      gSlvsTanTheta->SetPoint(ThetaBin, (tanThetaBinEdges[2*ThetaBin] + tanThetaBinEdges[2*ThetaBin+1]) / 2., p1/p0);
      gSlvsTanTheta->SetPointError(ThetaBin, (tanThetaBinEdges[2*ThetaBin+1] - tanThetaBinEdges[2*ThetaBin]) / 2., totErrorSlope); 
      gCurvvsTanTheta->SetPoint(ThetaBin, (tanThetaBinEdges[2*ThetaBin] + tanThetaBinEdges[2*ThetaBin+1]) / 2., p2/p0);
      gCurvvsTanTheta->SetPointError(ThetaBin, (tanThetaBinEdges[2*ThetaBin+1] - tanThetaBinEdges[2*ThetaBin]) / 2., totErrorCurv); 
    }
  } 
}