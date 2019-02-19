void SplineLoad() {
  AliTPCPIDResponse* resp = new AliTPCPIDResponse();
  resp->InitFromOADB(195517,4,"4","~/Documents/SplineCreation/LHC13c_pass4_MultDifferentiated/Iteration6/TPCPIDResponseOADB_2018_02_14_13C_pass4_It6.root");
  AliESDtrack* track = new AliESDtrack();
  Double_t dEdx = 0.0;
  Int_t nPoints = 0;
  AliTPCPIDResponse::ETPCgainScenario gain = 0;
  TSpline3* func = 0x0;
  resp->GetEnableMultSplines();
  TSpline3* func = 0x0
  TSpline3* func1 = 0x0
  resp->SetCurrentEventMultiplicity(100);
  resp->ChooseSplineForMultiplicity();  
  resp->ResponseFunctiondEdxN(track, AliPID::kProton, 0, dEdx, nPoints, gain, &func)

}