// For normal commits: specialVersion = kFALSE, commitVersion = kTRUE
Int_t createMapFile(Bool_t specialVersion = kFALSE, Bool_t commitVersion = kFALSE,
                    TString normalisation = "NoNormalisation"/*, TString pathMapPackage = "etaMaps",
                    TString fileNameMapPackage = "TPCetaMaps.root"*/)
{
  AliLog::SetGlobalLogLevel(AliLog::kWarning);
  
  gROOT->LoadMacro("addMapToFile.C+");
  
  // 2015 pp, 13 TeV, pass2
  addMapToFile("../pp13TeV/Parameters/outputCheckShapeEtaTree_2016_05_03__12_07__NewSplines___2016_04_29__08_39.root", normalisation, "LHC15f", 2, kFALSE, kFALSE);
  addMapToFile("../pp13TeV/Parameters/outputCheckShapeEtaTree_2016_05_03__12_09__CorrectedWithMap_outputCheckShapeEtaTree_2016_05_03__12_07__NewSplines___2016_04_29__08_39___2016_05_03__12_07.root", normalisation, "LHC15f", 2, kFALSE, kTRUE);

  return 0;
}
