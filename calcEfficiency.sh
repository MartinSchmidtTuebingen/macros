#!/bin/bash

#Par1 is etaAbs cut (currently 0.8)
#Par2 date - if set, sys err output will be processed. Otherwise: Usual output
#Par3 is the path to data
#Par4 is the path to MCData
#Par5 JetString processed
#Par6: lower Centrality
#par7: upper Centrality
#par8: Jet Loop?
#par9: Single Systematics or noJet?
#par10: MC truth file for UE
folder_for_macros=~/Documents/macros

source ali master

date=$2
charge=0
#charge=$3

nSigmaString=nSigma0.0

chargeString=""
if [ "x$charge" = "x-1" ]; then
  chargeString="_negCharge"
fi
if [ "x$charge" = "x1" ]; then
  chargeString="_posCharge"
fi

#Calculate errors, exclude multiplicity-dependent errors
sysErrTypeMC=6
# sysErrTypeMC=6  #For Merging without MultDep
# sysErrTypeMC=4  #For Merging

pathData=$3

fileEfficiency=$4

if [ ! -f $fileEfficiency ];then
  echo "Efficiency file not found"
  exit 0
fi

if [ "$5" != "" ];then
  jetString=$5
else
  jetString="Jets_PureGauss_"
fi

pathMCUETruth=${10}

pathSysErr="ONLY_USED_FOR_JETS"

#Bool_t correctGeantFluka, Bool_t newGeantFluka, Bool_t scaleStrangeness, Bool_t applyMuonCorrection, Int_t chargeMode,
stringArgs1="kTRUE,kFALSE,kTRUE,kTRUE,$charge"

#Int_t rebinEfficiencyObs,
stringArgs2="1"


# CentralityString: lowerCentralityData,upperCentralityData,lowerCentrality, upperCentrality

stringCent=$6,$7,$6,$7
if [ "$6" != "-2" ];then
  centString=_centrality$6"_"$7
else
  centString=""
fi

jetLoop=$8



#Double_t etaAbsCut,
#Double_t eps_trigger,typeMCSysErrors sysErrorTypeMC,Bool_t normaliseToNInel
#Assume eps_trigger = 1 for mult analysis (careful: May be wrong for peripheral bins)
stringArgs3="1.,0,kFALSE"
stringArgs3MB="1.,$sysErrTypeMC,kTRUE"
# stringArgs3MB="1.,0,kTRUE"



dateString=""
if [ "$date" = "" ]; then
  if [ $jetLoop = "true" ]; then
    filename1="maschmid_PID_Jets_PureGauss_results_LLFit__"
    filename2="2_reg1_regFac1.00_noMuons_idSpectra"
  else
    filename="maschmid_PID_Jets_PureGauss_results_LLFit__Pt_2_reg1_regFac1.00_noMuons_idSpectra_centrality60_100";
  fi  
  echo "No sys errors";
else
  if [ "$9" != "" ];then
    if [ "$9" = "noJet" ];then
      filename1="outputSystematicsTotal_SummedSystematicErrors"
    else
      filename1="outputSystematics_$jetString";
    fi
    dateString="__$date"
  else  
    filename1="outputSystematicsTotal_SummedSystematicErrors_$jetString";
    dateString="__$date""__$date";
  fi  
  echo "With sys errors";
fi


# With fixed centrality range for corr factors
corrMultLow=-2
corrMultUp=-2

# Expand if desired
#Bool_t correctMCID, Bool_t individualMultCorr, TString multCorrPathName, Bool_t correctEff10d10e, Bool_t isLowestMultBin
# stringArgsMB=$stirngArgs3MB,..."
#stringArgs3=$stringArgs3,kFALSE,kFALSE,'"'$fileMultDepSec'"',kTRUE
if [ "$jetLoop" = "true" ]; then
  echo "Processing Jet Loop"
  # Modes: 0=pt, 1=z, 2=xi, 3=R, 4=jT
  for mode in R Jt #pT Z xi R Jt
  do
    if [ "$mode" = "pT" ];then
      iOBS=0
    elif [ "$mode" = "z" ];then
      iOBS=1
    elif [ "$mode" = "xi" ];then
      iOBS=2
    elif [ "$mode" = "R" ];then
      iOBS=3
    elif [ "$mode" = "Jt" ];then
      iOBS=4
    fi
    for jetPtString in 10.0_15.0 15.0_20.0 20.0_30.0
    do
      if [ "$jetPtString" = "5.0_10.0" ];then
        lowerJetPt=5.0
        upperJetPt=10.0
      elif [ "$jetPtString" = "10.0_15.0" ];then
        lowerJetPt=10.0
        upperJetPt=15.0
      elif [ "$jetPtString" = "15.0_20.0" ];then
        lowerJetPt=15.0
        upperJetPt=20.0
      elif [ "$jetPtString" = "20.0_30.0" ];then
        lowerJetPt=20.0
        upperJetPt=30.0
      elif [ "$jetPtString" = "10.0_30.0" ];then
        lowerJetPt=10.0
        upperJetPt=30.0        
      else
        echo "Missing definition for lower/upper JetPt"
      fi 
      if [ "$9" != "" ];then
        filename=$filename1""$mode"_"$9"$centString""_jetPt"$jetPtString"_"$nSigmaString$chargeString$dateString.root
      else  
        filename=$filename1""$mode""_$filename2$centString""_jetPt$jetPtString$chargeString$dateString.root
      fi
      if [ ! -f $pathData/$filename ];then
        echo "File $pathData/$filename not found"
      else
        aliroot $folder_for_macros/calcEfficiency.C+"(\"$fileEfficiency\",\"$pathData/$filename\", \"$pathSysErr\",$stringArgs1,$stringCent,$lowerJetPt,$upperJetPt,$iOBS,0,$stringArgs2,$1,$stringArgs3,\"$pathMCUETruth\")" -l -b -q
      fi
    done
  done
else
  if [ "$9" != "" ];then
    if [ "$9" = "noJet" ];then
      filename=$filename1$centString$chargeString$dateString.root
    else
      filename=$filename1"Pt_"$9"$centString"_$nSigmaString$chargeString$dateString.root
    fi  
  else
    filename=$filename1""_Pt__$centString""$chargeString$dateString.root
  fi  
  if [ ! -f $pathData/$filename ];then
    echo "File $pathData/$filename not found"
  else  
    aliroot $folder_for_macros/calcEfficiency.C+"(\"$fileEfficiency\",\"$pathData/$filename\", \"$pathSysErr\",$stringArgs1,$stringCent,-1,-1,0,0,$stringArgs2,$1,$stringArgs3MB)" -l -b -q
  fi
fi  







#   filename="$filename1"Pt_2"$filename2"
#   aliroot 'calcEfficiency.C+("'$fileEfficiency'","'$pathData'/'$filename'_jetPt5.0_10.0'$chargeString$dateString'.root","'$pathSysErr'",'"$stringArgs1"','"$stringCent"','"$stringArgs2"','$1','"$stringArgs3MB"',kFALSE)' -l -b -q
#   aliroot 'calcEfficiency.C+("'$fileEfficiency'","'$pathData'/'$filename'_jetPt10.0_15.0'$chargeString$dateString'.root","'$pathSysErr'",'"$stringArgs1"','"$stringCent"','"$stringArgs2"','$1','"$stringArgs3MB"',kFALSE)' -l -b -q
#   aliroot 'calcEfficiency.C+("'$fileEfficiency'","'$pathData'/'$filename'_jetPt15.0_20.0'$chargeString$dateString'.root","'$pathSysErr'",'"$stringArgs1"','"$stringCent"','"$stringArgs2"','$1','"$stringArgs3MB"',kFALSE)' -l -b -q
#   aliroot 'calcEfficiency.C+("'$fileEfficiency'","'$pathData'/'$filename'_jetPt20.0_30.0'$chargeString$dateString'.root","'$pathSysErr'",'"$stringArgs1"','"$stringCent"','"$stringArgs2"','$1','"$stringArgs3MB"',kFALSE)' -l -b -q 
#   filename="$filename1"Z_2"$filename2"
#   aliroot 'calcEfficiency.C+("'$fileEfficiency'","'$pathData'/'$filename'_jetPt5.0_10.0'$chargeString$dateString'.root","'$pathSysErr'",'"$stringArgs1"','"$stringCent"','"$stringArgs2"','$1','"$stringArgs3MB"',kFALSE)' -l -b -q
#   aliroot 'calcEfficiency.C+("'$fileEfficiency'","'$pathData'/'$filename'_jetPt10.0_15.0'$chargeString$dateString'.root","'$pathSysErr'",'"$stringArgs1"','"$stringCent"','"$stringArgs2"','$1','"$stringArgs3MB"',kFALSE)' -l -b -q
#   aliroot 'calcEfficiency.C+("'$fileEfficiency'","'$pathData'/'$filename'_jetPt15.0_20.0'$chargeString$dateString'.root","'$pathSysErr'",'"$stringArgs1"','"$stringCent"','"$stringArgs2"','$1','"$stringArgs3MB"',kFALSE)' -l -b -q
#   aliroot 'calcEfficiency.C+("'$fileEfficiency'","'$pathData'/'$filename'_jetPt20.0_30.0'$chargeString$dateString'.root","'$pathSysErr'",'"$stringArgs1"','"$stringCent"','"$stringArgs2"','$1','"$stringArgs3MB"',kFALSE)' -l -b -q   
# else
#   aliroot 'calcEfficiency.C+("'$fileEfficiency'","'$pathData'/'$filename$chargeString$dateString'.root","'$pathSysErr'",'"$stringArgs1"','"$stringCent"','"$stringArgs2"','$1','"$stringArgs3MB"',kFALSE)' -l -b -q
# fi