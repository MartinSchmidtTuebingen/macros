enum modes { kPt = 0, kZ = 1, kXi = 2};

//_____________________________________________________________________________________________________________________________
Bool_t CentralityHasDecimalsPlaces(Double_t cent)
{
  // Check, whether the centrality has decimals places or is an integer
  const Double_t temp1 = ((Int_t)cent)*1e6;
  const Double_t temp2 = cent*1e6;
  
  return TMath::Abs(temp1 - temp2) > 1;
}


//_____________________________________________________________________________________________________________________________
void setMeanThreshold(Double_t& setMeanLowerThreshold, Double_t& setMeanUpperThreshold, Int_t mode, Double_t lowerJetPt, Double_t upperJetPt)
{
  const Double_t setMeanLowerThresholdPt = 4;
  const Double_t setMeanUpperThresholdPt = 999.;
  
  // Take effective jetPt for z and xi
  Double_t effectiveJetPt = 0.5 * (lowerJetPt + upperJetPt);
  
  if (mode == kPt) {
    setMeanLowerThreshold = setMeanLowerThresholdPt;
    setMeanUpperThreshold = setMeanUpperThresholdPt;
  }
  else if (mode == kZ) {
    setMeanLowerThreshold = setMeanLowerThresholdPt / effectiveJetPt;
    setMeanUpperThreshold = setMeanUpperThresholdPt / effectiveJetPt;
  }
  else if (mode == kXi) {
    // Thresholds are swapped!
    setMeanLowerThreshold = TMath::Log(effectiveJetPt / setMeanUpperThresholdPt);
    setMeanUpperThreshold = TMath::Log(effectiveJetPt / setMeanLowerThresholdPt);
    if (setMeanUpperThresholdPt > 900)
      setMeanLowerThreshold = 0.;
  }
  else {
    printf("ERROR: Unknown mode!\n");
    exit(-1);
  }
  
}
/* "", "_negCharge", "_posCharge"*/
void runSystematicErrorEstimation(TString jetString = "", TString chargeString = "", TString referencepath = "", TString sysuemethodpath = "", TString sysUEMethodString = "", TString syswoprePIDpath ="", TString sysshapepath = "" , TString sysshapeString = "", TString syssplinespath = "", TString syssigmapath = "", TString sysetapath = "", TString sysmultpath = "", TString outfilepath = "")
{
  // NOTE:
  // If mean instead of default is to be used, just set setMean to kTRUE in the following.
  // ALSO: In addUpSystematicErrors the reference must be adapted to the corresponding mean (not the default).
  // Currently, ONLY the graph and the plot with the single outputs for the sys errors are changed.
  // TODO for average instead of default for spline systetematics at high pT > 4 GeV/c:
  // - Seems to make no difference more or less (at least for the old "symmetric" systetematics)
  // - Questionable to set the threshold as function of momentum! Could mean that K and p are move away from the best estimated
  //   because only the different models only affect pions in the rel. rise but for quite some pT not the other species
  // - Also a problem:
  // ATTENTION: Normally, weighted average should be done either for the yields or should be WEIGHTED with the total yield
  // (must be done anyway for the to-pion-ratios!). BUT: The whole approach relies on the assumption that the overall statistics
  // is the same for all compared data points. This must be checked anyway and is true within percent typically!
  //
  // => Looks like much to do and many problems with most likely little gain. Better don't use the average instead of the default!
	printf("Start\n");
  gROOT->LoadMacro("~/Documents/macros/SystematicErrorEstimation.C+");
  
  const Double_t nSigma = 0;
  
  const Bool_t ignoreSigmaErrors = kTRUE;
  
  const Int_t numCentralities = 3;
  const Int_t centralities[2*numCentralities] = { 0, 100, 0, 10, 60, 100};
  
  const Int_t numCentralitiesRefMultOld = 9;
  const Int_t centralitiesRefMultOldLowEdge[numCentralitiesRefMultOld] = { 0,  7, 13, 20, 29, 40, 50, 60,   0 };
  const Int_t centralitiesRefMultOldUpEdge[numCentralitiesRefMultOld]  = { 7, 13, 20, 29, 40, 50, 60, 95, 125 };
  
  const Int_t numCentralitiesRefMultNew = 18;
  const Int_t centralitiesRefMultNewLowEdge[numCentralitiesRefMultNew] =
    { 1, 4,  7, 10, 15, 20, 25, 30, 40, 50, 60,  70,  100, 15,   25,   40,   60,    0 };
  const Int_t centralitiesRefMultNewUpEdge[numCentralitiesRefMultNew]  =
    { 4, 7, 10, 15, 20, 25, 30, 40, 50, 60, 70, 100, 9999, 25, 9999, 9999, 9999, 9999 };
  
  
  const Int_t numCentralitiesV0Mpp = 15;
  const Double_t centralitiesV0MppLowEdge[numCentralitiesV0Mpp] = {    0, 0.01, 0.1, 1,  5, 10, 15, 20, 30, 40, 50,  70,   0, 0,   0 };
  const Double_t centralitiesV0MppUpEdge[numCentralitiesV0Mpp]  = { 0.01,  0.1,   1, 5, 10, 15, 20, 30, 40, 50, 70, 100, 0.1, 1, 100 };
  
  TString centralityString = "";
  
  const Int_t numJetPtBins = 4;
  const Int_t jetPt[numJetPtBins + 1] = { 5, 10, 15, 20, 30};
  
//   const Int_t numJetPtBins = 1;
//   const Int_t jetPt[numJetPtBins + 1] = { 10, 30};
  
  TString jetPtString = "";
  
  Double_t setMeanLowerThreshold = 999.;
  Double_t setMeanUpperThreshold = 999.;

  Bool_t useCentralities = kTRUE;
  Bool_t useJetPt = !(jetString.Contains("Inclusive")) && jetString.Contains("Jets");

  const TString modeString[3] = {"Pt", "Z", "Xi"};
  const Int_t maxJetPtBins = useJetPt ? numJetPtBins : 1;
  const Int_t maxMode = useJetPt ? 2 : 1;
  // Inclusive means: Only pT // In case of taking the average of the systematics: Only do this within certain range, HERE: pT only
  if (useJetPt) 
    setMeanThreshold(setMeanLowerThreshold, setMeanUpperThreshold, kPt, -1, -1);
  
  for (Int_t iCent = 0; iCent < numCentralities; iCent++) {
    centralityString = useCentralities ? Form("_centrality%d_%d", centralities[2*iCent], centralities[2*iCent + 1]) : ""; 

    for (Int_t iJetPt = 0; iJetPt <= maxJetPtBins-1; iJetPt++)  {
      if (iJetPt == maxJetPtBins) {
        // 10-40 bin
        iJetPt++;
        }
      jetPtString = useJetPt ? Form("_jetPt%d.0_%d.0", jetPt[iJetPt], jetPt[iJetPt + 1]) : "";
      jetPtString.Append(chargeString);        
      for (Int_t mode = 0; mode <= maxMode-1; mode++) {  
        if (!sysuemethodpath.IsNull())
        {
          Bool_t setMean = kFALSE;
          const Int_t numFiles = 2;
          
          TString outFileTitle = Form("Method%s%s", centralityString.Data(), jetPtString.Data());
          
          TString fileNames[numFiles];
          fileNames[0] = Form("%s/maschmid_PID_%sresults_LLFit__%s_2_reg1_regFac1.00_noMuons_idSpectra%s%s.root", referencepath.Data(), jetString.Data(), modeString[mode].Data(), centralityString.Data(), jetPtString.Data());
          fileNames[1] = Form("%s/maschmid_PID_%sresults_LLFit__%s_2_reg1_regFac1.00_noMuons_idSpectra%s%s.root", sysuemethodpath.Data(), sysUEMethodString.Data(), modeString[mode].Data(), centralityString.Data(), jetPtString.Data());
          
          TString histTitles[numFiles];
          histTitles[0] = "Perpendicular Cone (Default)";
          histTitles[1] = "Random Cone";
                  
          outFileTitle = Form("%s%s_%s", jetString.Data(), modeString[mode].Data(), outFileTitle.Data());
          SystematicErrorEstimation(outfilepath, outFileTitle, fileNames, histTitles, numFiles, nSigma, ignoreSigmaErrors, setMean, setMeanLowerThreshold, setMeanUpperThreshold);  
        }
        if (!syswoprePIDpath.IsNull())
        {
          Bool_t setMean = kFALSE;
          const Int_t numFiles = 2;
          
          TString outFileTitle = Form("PrePID%s%s", centralityString.Data(), jetPtString.Data());
          
          TString fileNames[numFiles];
          fileNames[0] = Form("%s/maschmid_PID_%sresults_LLFit__%s_2_reg1_regFac1.00_noMuons_idSpectra%s%s.root", referencepath.Data(), jetString.Data(), modeString[mode].Data(), centralityString.Data(), jetPtString.Data());
          fileNames[1] = Form("%s/maschmid_PID_%sresults_LLFit__%s_2_reg1_regFac1.00_noMuons%s%s.root", syswoprePIDpath.Data(), jetString.Data(), modeString[mode].Data(), centralityString.Data(), jetPtString.Data());

          TString histTitles[numFiles];
          histTitles[0] = "PID combined (default)";
          histTitles[1] = "No PID";
  
                  
          outFileTitle = Form("%s%s_%s", jetString.Data(), modeString[mode].Data(), outFileTitle.Data());
          SystematicErrorEstimation(outfilepath, outFileTitle, fileNames, histTitles, numFiles, nSigma, ignoreSigmaErrors, setMean, setMeanLowerThreshold, setMeanUpperThreshold);  
        }
        if (!sysshapepath.IsNull())
        {
          Bool_t setMean = kFALSE;
          const Int_t numFiles = 2;
          
          TString outFileTitle = Form("Shape%s%s", centralityString.Data(), jetPtString.Data());
          
          TString fileNames[numFiles];
          fileNames[0] = Form("%s/maschmid_PID_%sresults_LLFit__%s_2_reg1_regFac1.00_noMuons_idSpectra%s%s.root", referencepath.Data(), jetString.Data(), modeString[mode].Data(), centralityString.Data(), jetPtString.Data());
          fileNames[1] = Form("%s/maschmid_PID_%sresults_LLFit__%s_2_reg1_regFac1.00_noMuons_idSpectra%s%s.root", sysshapepath.Data(), sysshapeString.Data(), modeString[mode].Data(), centralityString.Data(), jetPtString.Data());

          TString histTitles[numFiles];
          histTitles[0] = "Pure Gauss (default)";
          histTitles[1] = "Asymmetric shape ";
          
          outFileTitle = Form("%s%s_%s", jetString.Data(), modeString[mode].Data(), outFileTitle.Data());
          SystematicErrorEstimation(outfilepath, outFileTitle, fileNames, histTitles, numFiles, nSigma, ignoreSigmaErrors, setMean, setMeanLowerThreshold,setMeanUpperThreshold);  
        }
        if (!syssigmapath.IsNull())
        {
          Bool_t setMean = kFALSE;
          const Int_t numFiles = 3;
          
          TString outFileTitle = Form("Sigma%s%s", centralityString.Data(), jetPtString.Data());
          
          TString fileNames[numFiles];
          fileNames[0] = Form("%s/maschmid_PID_%sresults_LLFit__%s_2_reg1_regFac1.00_noMuons_idSpectra%s%s.root", referencepath.Data(), jetString.Data(), modeString[mode].Data(), centralityString.Data(), jetPtString.Data());
          fileNames[1] = Form("%s/maschmid_PID_%sSystematicsSigmaDown_results_LLFit__%s_2_reg1_regFac1.00_noMuons_idSpectra%s%s.root", syssigmapath.Data(), jetString.Data(), modeString[mode].Data(), centralityString.Data(), jetPtString.Data());
          fileNames[2] = Form("%s/maschmid_PID_%sSystematicsSigmaUp_results_LLFit__%s_2_reg1_regFac1.00_noMuons_idSpectra%s%s.root", syssigmapath.Data(), jetString.Data(), modeString[mode].Data(), centralityString.Data(), jetPtString.Data());

          TString histTitles[numFiles];
          histTitles[0] = "Standard #sigma map";
          histTitles[1] = "#sigma correction - 4%";
          histTitles[2] = "#sigma correction + 4%";
          
          outFileTitle = Form("%s%s_%s", jetString.Data(), modeString[mode].Data(), outFileTitle.Data());
          SystematicErrorEstimation(outfilepath, outFileTitle, fileNames, histTitles, numFiles, nSigma, ignoreSigmaErrors, setMean, setMeanLowerThreshold,setMeanUpperThreshold);  
        }
        if (!sysetapath.IsNull())
        {
          Bool_t setMean = kFALSE;
          const Int_t numFiles = 3;
          
          TString outFileTitle = Form("Eta%s%s", centralityString.Data(), jetPtString.Data());
          
          TString fileNames[numFiles];
          fileNames[0] = Form("%s/maschmid_PID_%sresults_LLFit__%s_2_reg1_regFac1.00_noMuons_idSpectra%s%s.root", referencepath.Data(), jetString.Data(), modeString[mode].Data(), centralityString.Data(), jetPtString.Data());
          fileNames[1] = Form("%s/maschmid_PID_%sSystematicsEtaDown_results_LLFit__%s_2_reg1_regFac1.00_noMuons_idSpectra%s%s.root", sysetapath.Data(), jetString.Data(), modeString[mode].Data(), centralityString.Data(), jetPtString.Data());
          fileNames[2] = Form("%s/maschmid_PID_%sSystematicsEtaUp_results_LLFit__%s_2_reg1_regFac1.00_noMuons_idSpectra%s%s.root", sysetapath.Data(), jetString.Data(), modeString[mode].Data(), centralityString.Data(), jetPtString.Data());

          TString histTitles[numFiles];
          histTitles[0] = "Standard #eta map";
          histTitles[1] = "#eta correction - 1.5%, p > 0.45 GeV/c; - %3, p < 0.45 GeV/c";
          histTitles[2] = "#eta correction + 1.5%, p > 0.45 GeV/c; + %3, p < 0.45 GeV/c";
          
          outFileTitle = Form("%s%s_%s", jetString.Data(), modeString[mode].Data(), outFileTitle.Data());
          SystematicErrorEstimation(outfilepath, outFileTitle, fileNames, histTitles, numFiles, nSigma, ignoreSigmaErrors, setMean, setMeanLowerThreshold,setMeanUpperThreshold);  
        }
        if (!syssplinespath.IsNull())
        {
          Bool_t setMean = kFALSE;
          const Int_t numFiles = 3;
                  
          TString outFileTitle = Form("Splines%s%s", centralityString.Data(), jetPtString.Data());
          TString fileNames[numFiles];
          fileNames[0] = Form("%s/maschmid_PID_%sresults_LLFit__%s_2_reg1_regFac1.00_noMuons_idSpectra%s%s.root", referencepath.Data(), jetString.Data(), modeString[mode].Data(), centralityString.Data(), jetPtString.Data());
          fileNames[1] = Form("%s/maschmid_PID_%sSystematicsSplinesDown_results_LLFit__%s_2_reg1_regFac1.00_noMuons_idSpectra%s%s.root", syssplinespath.Data(), jetString.Data(), modeString[mode].Data(), centralityString.Data(), jetPtString.Data());
          fileNames[2] = Form("%s/maschmid_PID_%sSystematicsSplinesUp_results_LLFit__%s_2_reg1_regFac1.00_noMuons_idSpectra%s%s.root", syssplinespath.Data(), jetString.Data(), modeString[mode].Data(), centralityString.Data(), jetPtString.Data());

          TString histTitles[numFiles];
          histTitles[0] = "Standard splines";
          histTitles[1] = "Splines -0.6% - -1.4% for #beta#gamma > 50; -0.3% for #beta#gamma < 50%";
          histTitles[2] = "Splines +0.8% - +2% for #beta#gamma > 50; +0.3% for #beta#gamma < 50%";
          
          outFileTitle = Form("%s%s_%s", jetString.Data(), modeString[mode].Data(), outFileTitle.Data());
          SystematicErrorEstimation(outfilepath, outFileTitle, fileNames, histTitles, numFiles, nSigma, ignoreSigmaErrors, setMean, setMeanLowerThreshold,setMeanUpperThreshold);
        }
        if (!sysmultpath.IsNull())
        {
          Bool_t setMean = kFALSE;
          const Int_t numFiles = 3;
                  
          TString outFileTitle = Form("Mult%s%s", centralityString.Data(), jetPtString.Data());
          TString fileNames[numFiles];
          fileNames[0] = Form("%s/maschmid_PID_%sresults_LLFit__%s_2_reg1_regFac1.00_noMuons_idSpectra%s%s.root", referencepath.Data(), jetString.Data(), modeString[mode].Data(), centralityString.Data(), jetPtString.Data());
          fileNames[1] = Form("%s/maschmid_PID_%sSystematicsMultDown_results_LLFit__%s_2_reg1_regFac1.00_noMuons_idSpectra%s%s.root", sysmultpath.Data(), jetString.Data(), modeString[mode].Data(), centralityString.Data(), jetPtString.Data());
          fileNames[2] = Form("%s/maschmid_PID_%sSystematicsMultUp_results_LLFit__%s_2_reg1_regFac1.00_noMuons_idSpectra%s%s.root", sysmultpath.Data(), jetString.Data(), modeString[mode].Data(), centralityString.Data(), jetPtString.Data());

          TString histTitles[numFiles];
          histTitles[0] = "Standard Multiplicity correction";
          histTitles[1] = "Splines + 0.2%";
          histTitles[2] = "Splines - 0.2%";
          
          outFileTitle = Form("%s%s_%s", jetString.Data(), modeString[mode].Data(), outFileTitle.Data());
          SystematicErrorEstimation(outfilepath, outFileTitle, fileNames, histTitles, numFiles, nSigma, ignoreSigmaErrors, setMean, setMeanLowerThreshold, setMeanUpperThreshold);
        }
/*          {
          Bool_t setMean = kFALSE;
          const Int_t numFiles = 3;
          
          TString outFileTitle = Form("MethodAverageRef%s", jetPtString.Data());
          
          TString fileNames[numFiles];
          fileNames[0] = 
                  Form("%sResults/maschmid_PID_%s%sresults_LLFit__%s_2_reg1_regFac1.00_noMuons_idSpectra%s.root"
                                  , path.Data(), jetString.Data(), shaperef.Data(), modeString[mode].Data(), jetPtString.Data());
          fileNames[1] = 
                  Form("%sResults/maschmid_PID_Jets_UE_PC_%sresults_LLFit__%s_2_reg1_regFac1.00_noMuons_idSpectra%s.root"
                                  , path.Data(), shaperef.Data(), modeString[mode].Data(), jetPtString.Data());
          fileNames[2] = 
                  Form("%sResults/maschmid_PID_Jets_UE_RC_%sresults_LLFit__%s_2_reg1_regFac1.00_noMuons_idSpectra%s.root"
                                  , path.Data(), shaperef.Data(), modeString[mode].Data(), jetPtString.Data());
                  
          TString histTitles[numFiles];
          histTitles[0] = "Averaged Result of Methods";
          histTitles[1] = "Perpendicular Cone";
          histTitles[2] = "Random Cone";
  
                  
          outFileTitle = Form("%s%s_%s", jetString.Data(), modeString[mode].Data(), outFileTitle.Data());
          SystematicErrorEstimation(path, outFileTitle, fileNames, histTitles, numFiles, nSigma, ignoreSigmaErrors, setMean, setMeanLowerThreshold, setMeanUpperThreshold);  
        }    */          
      }
    }
  }
}