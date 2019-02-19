#ifndef __CINT__
#include "TFileMerger.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TString.h"
#include "TSystem.h"
#endif

void recursiveMerging(TString pathNameWithWildCards = "/lustre/nyx/alice/users/maschmid/train/V012.pPb/2016-11-14_2315.23840/mergedRuns/pPb/5.023ATeV/LHC13c.pass2/1*.ana/maschmid_PID_Jets_UE_PC_SystematicsEtaDown.root", TString outputFile = "/lustre/nyx/alice/users/maschmid/train/V012.pPb/2016-11-14_2315.23840/mergedPeriods/pPb/5.023ATeV/LHC13c.pass2/maschmid_PID_Jets_UE_PC_SystematicsEtaDown", Int_t nmerge = 5, Int_t startwithindex = 0/*, 
                      Bool_t cleanUpIntermediateFiles = kTRUE*/) {
  
  //e.g. pathNameWithWildCards = "/lustre/alice/wiechula/train/PIDtests/2012-04-26_1312.13238.newMaps/mergedRuns/pp/7TeV/LHC10d.pass2/125*.ana/TPCresidualPID.root";
  //e.g. pathNameWithWildCards = "/hera/alice/wiechula/train/V009.pp/2012-04-26_1312.13238/mergedRuns/pp/7TeV/LHC10d.pass2/125*.ana/TPCresidualPID.root";
  TString s=gSystem->GetFromPipe(Form("ls %s", pathNameWithWildCards.Data()));

  TObjArray *arr=s.Tokenize("\n");

  Int_t depth = 0;
  
  printf("List of files:\n");
  for (Int_t ifile = 0; ifile < arr->GetEntriesFast(); ifile++) {
    if (arr->At(ifile)) {
      Bool_t skip = ifile < startwithindex;
      printf("%d (%d)%s: %s\n", (Int_t)(ifile/nmerge), ifile, skip  ? " skipped" : "", arr->At(ifile)->GetName());
      
      if (skip)
        arr->RemoveAt(ifile);
    }
  }
  
  printf("\n\n");
  
  arr->Compress();
  
  
  
  
  printf("List of files after removal:\n");
  for (Int_t ifile = 0; ifile < arr->GetEntriesFast(); ifile++) {
    if (arr->At(ifile)) {
      printf("%d (%d): %s\n", (Int_t)(ifile/nmerge), ifile, arr->At(ifile)->GetName());
    }
  }
  
  printf("\n\n");
  
  
  
  while (arr->GetEntries() > 1){
    printf("depth: %d\n",depth);
    for (Int_t iIter=0; iIter<TMath::Ceil((Double_t)arr->GetEntries()/((Double_t)nmerge)); ++iIter){
      printf("Iter: %d\n",iIter);
      TFileMerger m(0);
      m.OutputFile(Form("%s.%d.%d.root",outputFile.Data(),depth,iIter));
      printf("writing output file: %s\n", Form("%s.%d.%d.root",outputFile.Data(),depth,iIter));
      for (Int_t ifile=iIter*nmerge; ifile<(iIter+1)*nmerge; ++ifile){
        if (!arr->At(ifile)) continue;
        printf("Adding file: %s\n",arr->At(ifile)->GetName());
        m.AddFile(arr->At(ifile)->GetName());
      }
      m.Merge();
    }
    delete arr;
    arr=0x0;
    s=gSystem->GetFromPipe(Form("ls %s.%d.[0-9]*.root",outputFile.Data(),depth));
    arr=s.Tokenize("\n");
    ++depth;
    printf("\n-----------\n");
  }
  
  gSystem->Exec(Form("mv %s.%d.0.root %s.root",outputFile.Data(),depth-1,outputFile.Data()));
  /*
  if (cleanUpIntermediateFiles) {
    delete arr;
    arr=0x0;
    s=gSystem->GetFromPipe(Form("ls %s.[0-9]*.[0-9]*.root",outputPathName.Data()));
    arr=s.Tokenize("\n");
    
    for (Int_t i=0; i<arr->GetEntries(); ++i){
      printf("Removing file: %s\n", arr->At(i)->GetName());
      gSystem->Exec(Form("rm %s", arr->At(i)->GetName()));
    }
  }*/
}

