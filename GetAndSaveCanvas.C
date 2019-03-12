void GetAndSaveCanvas(TString filename, TString savename) {
  TString canvasName = "cCorrDataToPiRatioMCTruth";
//   TString canvasName = "cCorrFractionsWithMC";
  gROOT->LoadMacro("GetObjectOutOfDirectory.C");
  gROOT->LoadMacro("changeCanvas.C");
  TCanvas* c = (TCanvas*)GetObjectOutOfDirectory(filename, canvasName)->Clone();
  changeCanvas(c);
  c->SaveAs(savename.Data());
  delete c;
}