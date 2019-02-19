//_____________________________________________________________________________________________________________________________
Bool_t CentralityHasDecimalsPlaces(Double_t cent)
{
  // Check, whether the centrality has decimals places or is an integer
  const Double_t temp1 = ((Int_t)cent)*1e6;
  const Double_t temp2 = cent*1e6;
  
  return TMath::Abs(temp1 - temp2) > 1;
}


//_____________________________________________________________________________________________________________________________
void runAddUpSystematicErrors(Bool_t useCentralities = kFALSE,
                              TString referencepath = "",
                              TString date = "2013_12_11",
                              TString jetString = "Jets_",
                              TString chargeString = "",
                              TString path = "/hera/alice/maschmid/analysis/10d_e.pass2_merged/nclCut",
                              Bool_t multRef = kTRUE,
                              Bool_t newRefBinning = kTRUE)
{
  gROOT->LoadMacro("~/Documents/macros/AddUpSystematicErrors.C+");
  
  const Double_t nSigma = 0;
  
  const TString sigmaString = Form("_nSigma%.1f", nSigma);// nSigma > 0 ? Form("_nSigma%.1f", nSigma) : "";
  
  //const TString sigmaString = Form("%s_nSigma%.1f", chargeString.Data(),  nSigma);// nSigma > 0 ? Form("_nSigma%.1f", nSigma) : "";

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
  const Double_t centralitiesV0MppLowEdge[numCentralitiesV0Mpp] = {    0, 0.01, 0.1, 1,  5, 10, 15, 20, 30, 40, 50,  70,   0, 0,   0};
  const Double_t centralitiesV0MppUpEdge[numCentralitiesV0Mpp]  = { 0.01,  0.1,   1, 5, 10, 15, 20, 30, 40, 50, 70, 100, 0.1, 1, 100 };
  
  TString centralityString = "";
  
  const Int_t numJetPtBins = 4;
  const Int_t jetPt[numJetPtBins + 1] = { 5, 10, 15, 20, 30};
  
//   const Int_t numJetPtBins = 1;
//   const Int_t jetPt[numJetPtBins + 1] = { 10, 30};
  
  TString jetPtString = "";
  
  Bool_t useJetPt = !(jetString.Contains("Inclusive")) && jetString.Contains("Jets");

  const TString modeString[3] = {"Pt", "Z", "Xi"};
  const Int_t maxJetPtBins = useJetPt ? numJetPtBins : 1;
  const Int_t maxMode = useJetPt ? 2 : 1;  
  
  //TString path = "/hera/alice/maschmid/analysis/13b.pass3";  
  //TString date = "2013_09_13";

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
    
        const Int_t numFiles = 5;
        
        TString outFileTitle = Form("SummedSystematicErrors_%s_%s__%s%s__%s", jetString.Data(), modeString[mode].Data(), centralityString.Data(), jetPtString.Data(), date.Data());
        
        TString fileNames[numFiles];
        fileNames[0] = Form("%s/outputSystematics_%s%s_Splines%s%s%s__%s.root", path.Data(), jetString.Data(), modeString[mode].Data(), centralityString.Data(), jetPtString.Data(), sigmaString.Data(), date.Data());
        fileNames[1] = Form("%s/outputSystematics_%s%s_Sigma%s%s%s__%s.root", path.Data(), jetString.Data(), modeString[mode].Data(), centralityString.Data(), jetPtString.Data(), sigmaString.Data(), date.Data());
//         fileNames[2] = Form("%s/outputSystematics_%s%s_Mult%s%s%s__%s.root", path.Data(), jetString.Data(), modeString[mode].Data(), centralityString.Data(), jetPtString.Data(), sigmaString.Data(), date.Data());
        fileNames[2] = Form("%s/outputSystematics_%s%s_Eta%s%s%s__%s.root", path.Data(), jetString.Data(), modeString[mode].Data(), centralityString.Data(), jetPtString.Data(), sigmaString.Data(), date.Data());
        fileNames[3] = Form("%s/outputSystematics_%s%s_Shape%s%s%s__%s.root", path.Data(), jetString.Data(), modeString[mode].Data(), centralityString.Data(), jetPtString.Data(), sigmaString.Data(), date.Data());
        fileNames[4] = Form("%s/outputSystematics_%s%s_PrePID%s%s%s__%s.root", path.Data(), jetString.Data(), modeString[mode].Data(), centralityString.Data(), jetPtString.Data(), sigmaString.Data(), date.Data());
//         fileNames[6] = Form("%s/outputSystematics_%s%s_Method%s%s%s__%s.root", path.Data(), jetString.Data(), modeString[mode].Data(), centralityString.Data(), jetPtString.Data(), sigmaString.Data(), date.Data()); 

        TString fileNameReference = Form("%s/maschmid_PID_%sresults_LLFit__%s_2_reg1_regFac1.00_noMuons_idSpectra%s%s.root", referencepath.Data(), jetString.Data(), modeString[mode].Data(), centralityString.Data(), jetPtString.Data());
        
        if (jetString.Contains("UE"))
          AddUpSystematicErrors(path, outFileTitle, fileNames, numFiles, fileNameReference);  
        else
          AddUpSystematicErrors(path, outFileTitle, fileNames, numFiles-1, fileNameReference);
      }
    }
    
    if (!useCentralities)
      break;
  }
  
  
//   // pp inclusive/Jets
//   
//   
//   TString reference = "";
//   
//   if (noMuons) {
//     reference =  useJetPt ? Form("/Results/maschmid_PID_%sresults_LLFit__Pt_2_reg1_regFac1.00_noMuons_idSpectra", jetString.Data())
//                           : Form("maschmid_PID_PureGauss_results_LLFit__Pt_2_reg1_regFac1.00_noMuons_idSpectra%s.root",
//                                  chargeString.Data());
//   }
//   else {
//     reference =  useJetPt ? "maschmid_PID_Jets_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCHybridTrackCutsJets_idSpectra"
//               : Form("maschmid_PID_Jets_Inclusive_results_LLFit__Pt_2_reg1_regFac1.00_muonToElTunedOnMCHybridTrackCuts_idSpectra%s.root",
//                      chargeString.Data());
//   }
//   
//   TString referenceCurrent = reference;
//   
//   TString modeString[3] = { "_pT", "_z", "_xi" };
//   TString modeString2[3] = { "Pt", "Z", "Xi" };
//   
//   for (Int_t mode = -1; mode <= 1; mode++) {
//     if (mode == -1 && useJetPt)
//       continue;
//     
//     if (!useJetPt && mode >= 0)
//       break;
//     
//     for (Int_t iJetPt = 0; iJetPt <= numJetPtBins-1; iJetPt++) {
//       if (iJetPt == numJetPtBins) {
//         // 10-40 bin
//         iJetPt++;
//       }
//       
//       TString modeTotString = "";
//       if (mode >= 0 && !useDirectFit) 
//         modeTotString = Form("__centrality_all_jetPt%d.0_%d.0%s", jetPt[iJetPt], jetPt[iJetPt + 1], modeString[mode].Data());
//       
//       jetPtString = "";
//       if (mode >= 0) // => useJetPt
//         jetPtString = Form("_jetPt%d.0_%d.0%s", jetPt[iJetPt], jetPt[iJetPt + 1], modeTotString.Data());
//       
//       jetPtString.Append(chargeString);
//       
//       if (mode >= 0) {// => useJetPt
//         if (useDirectFit) {
//           // TODO no muons for direct fit vs. z, xi yet!
//           referenceCurrent = Form("Results/maschmid_PID_%sresults_LLFit__%s_2_reg1_regFac1.00_%s_idSpectra%s.root", jetString.Data(), 
//                                   modeString2[mode].Data(),
//                                   (mode == 0 && !noMuons) ? "muonToElTunedOnMCHybridTrackCutsJets" : "noMuons",
//                                   jetPtString.Data());
//         }
//         else
//           referenceCurrent =  Form("output_extractedFFs_%s_jetPt%d.0_%d.0__centrality_all_jetPt%d.0_%d.0.root", reference.Data(),
//                                   jetPt[iJetPt], jetPt[iJetPt + 1], jetPt[iJetPt], jetPt[iJetPt + 1]);
//       }
//       
//       TString prefix = "outputSystematics";
//       if (mode >= 0 && useDirectFit)
//         prefix = Form("%s_%s%s", prefix.Data(), jetString.Data(), modeString2[mode].Data());
//       
//       const Int_t numFilesMax = 8;
//       
//       Int_t numFiles = 1;
//       // TODO no muons for direct fit vs. z, xi yet!
//       if ((!useDirectFit || mode < 1) && !noMuons)
//         numFiles++;
//       
//       TString modeSuffix = "";
//       if (mode >= 0 && useDirectFit)
//         modeSuffix = Form("_%s_", modeString2[mode].Data());
//       
//       TString outFileTitle = Form("SummedSystematicErrors_%s%s%s", jetString.Data(), modeSuffix.Data(), jetPtString.Data());
//       
//       TString fileNames[numFilesMax];
//       fileNames[0] = Form("%s/%s_MethodAverageRef%s%s__%s.root", path.Data(), prefix.Data(), jetPtString.Data(),
//                           sigmaString.Data(), date.Data());
// //       fileNames[0] = Form("%s/%s_Splines%s%s__%s.root", path.Data(), prefix.Data(), jetPtString.Data(),
// //                           sigmaString.Data(), date.Data());
// //       fileNames[1] = Form("%s/%s_Sigma%s%s__%s.root", path.Data(), prefix.Data(), jetPtString.Data(),
// //                           sigmaString.Data(), date.Data());
// //       fileNames[2] = Form("%s/%s_Eta%s%s__%s.root", path.Data(), prefix.Data(), jetPtString.Data(),
// //                           sigmaString.Data(), date.Data());
// //       fileNames[3] = Form("%s/%s_Shape%s%s__%s.root", path.Data(), prefix.Data(), jetPtString.Data(),
// //                           sigmaString.Data(), date.Data());
// //       fileNames[4] = Form("%s/%s_PrePID%s%s__%s.root", path.Data(), prefix.Data(), jetPtString.Data(),
// //                           sigmaString.Data(), date.Data());
// //       fileNames[5] = Form("%s/%s_Mult%s%s__%s.root", path.Data(), prefix.Data(), jetPtString.Data(),
// //                           sigmaString.Data(), date.Data());
// //       fileNames[6] = Form("%s/%s_Method%s%s__%s.root", path.Data(), prefix.Data(), jetPtString.Data(),
// //                           sigmaString.Data(), date.Data());      
// //       if ((!useDirectFit || mode < 1) && !noMuons)
// //         fileNames[7] = Form("%s/%s_Muons%s%s__%s.root", path.Data(), prefix.Data(), jetPtString.Data(),
// //                             sigmaString.Data(), date.Data());
//         
//       TString fileNameReference = Form("%s/%s", path.Data(), referenceCurrent.Data());
//     
//       AddUpSystematicErrors(path, outFileTitle, fileNames, numFiles, fileNameReference, mode, useDirectFit); 
//       
//       if (!useJetPt)
//         break;
//     }
//   }
//   
  
  
  
  /*
  // multRef/V0M: pp inclusive
  // e.g. a 'runAddUpSystematicErrors.C(kTRUE, kFALSE, kTRUE, kTRUE, "2014_07_18", "", "/hera/alice/maschmid/analysis/10d_e.pass2_merged/nclCut/multiplicityDependencePP/trackletsMult/pileUpRejectionSPD/")' -b -q
  useCentralities = kTRUE;
  useJetPt = kFALSE;
  useDirectFit = kTRUE;
  noMuons = kTRUE;
  
  Int_t numCentralitiesRefMult = newRefBinning ? numCentralitiesRefMultNew : numCentralitiesRefMultOld;
  Int_t* centralitiesRefMultLowEdge = 0x0;
  Int_t* centralitiesRefMultUpEdge = 0x0;
  
  if (multRef) {
    centralitiesRefMultLowEdge = new Int_t[numCentralitiesRefMult];
    centralitiesRefMultUpEdge = new Int_t[numCentralitiesRefMult];
    
    for (Int_t i = 0; i < numCentralitiesRefMult; i++) {
      centralitiesRefMultLowEdge[i] = newRefBinning ? centralitiesRefMultNewLowEdge[i] : centralitiesRefMultOldLowEdge[i];
      centralitiesRefMultUpEdge[i]  = newRefBinning ? centralitiesRefMultNewUpEdge[i]  : centralitiesRefMultOldUpEdge[i];
    }
  }
  
  Int_t mode = 0;
  
  const Int_t numMults = multRef ? numCentralitiesRefMult : numCentralitiesV0Mpp;
  for (Int_t iMult = -1; iMult < numMults; iMult++) {
    if (multRef) {
      if (iMult == -1)
        centralityString = "";
      else
        centralityString = Form("_centrality%d_%d", centralitiesRefMultLowEdge[iMult], centralitiesRefMultUpEdge[iMult]);
    }
    else {
      centralityString = "";
      
      if (iMult >= 0) {
        Bool_t centralityHasDecimalsPlaces = CentralityHasDecimalsPlaces(centralitiesV0MppLowEdge[iMult]) ||
                                             CentralityHasDecimalsPlaces(centralitiesV0MppUpEdge[iMult]);
        
        centralityString = Form(centralityHasDecimalsPlaces ? "_centrality%.0fem2_%.0fem2" : "_centrality%.0f_%.0f", 
                                centralityHasDecimalsPlaces ? 100.*centralitiesV0MppLowEdge[iMult] : centralitiesV0MppLowEdge[iMult],
                                centralityHasDecimalsPlaces ? 100.*centralitiesV0MppUpEdge[iMult] : centralitiesV0MppUpEdge[iMult]);
      }
    }
    centralityString.Append(chargeString);
    
    TString reference = Form("maschmid_PID_TPCdefaultPriors_ITS_results_LLFit__Pt_2_reg1_regFac1.00_noMuons_idSpectra%s.root",
                             centralityString.Data());
    
    
    TString prefix = "outputSystematics";
    
    const Int_t numFilesMax = 5;
    
    Int_t numFiles = 5;
    
    TString outFileTitle = Form("SummedSystematicErrors%s", centralityString.Data());
    
    TString fileNames[numFilesMax];
    fileNames[0] = Form("%s/%s_Splines%s%s__%s.root", path.Data(), prefix.Data(), centralityString.Data(),
                        sigmaString.Data(), date.Data());
    fileNames[1] = Form("%s/%s_Sigma%s%s__%s.root", path.Data(), prefix.Data(), centralityString.Data(),
                        sigmaString.Data(), date.Data());
    fileNames[2] = Form("%s/%s_Eta%s%s__%s.root", path.Data(), prefix.Data(), centralityString.Data(),
                        sigmaString.Data(), date.Data());
    fileNames[3] = Form("%s/%s_Shape%s%s__%s.root", path.Data(), prefix.Data(), centralityString.Data(),
                        sigmaString.Data(), date.Data());
    fileNames[4] = Form("%s/%s_PrePID%s%s__%s.root", path.Data(), prefix.Data(), centralityString.Data(),
                        sigmaString.Data(), date.Data());
      
    TString fileNameReference = Form("%s/%s", path.Data(), reference.Data());
  
    AddUpSystematicErrors(path, outFileTitle, fileNames, numFiles, fileNameReference, mode, useDirectFit); 
  }
  */
}