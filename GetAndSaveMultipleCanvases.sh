source ali master

for FILE in /home/martin/Documents/Promotion/Analyses/LHC13c.pass4_JetShape/Data/output_EfficiencyCorrection*__Jt*
do
  saveName=/home/martin/Documents/Promotion/Talks/2019_03_JEMeeting_JetShape/img/LHC13c_pass4_ToPiRatios
  if [[ "$FILE" = *"__Jt_2"* ]];then
    saveName=$saveName"_jT"
  elif [[ "$FILE" = *"__R_2"* ]];then
    saveName=$saveName"_R"
  fi
  if [[ "$FILE" = *"_jetPt10.0_15.0"* ]];then
    saveName=$saveName"_jetpT10_15"
  elif [[ "$FILE" = *"_jetPt15.0_20.0"* ]];then
    saveName=$saveName"_jetpT15_20"
  elif [[ "$FILE" = *"_jetPt20.0_30.0"* ]];then
    saveName=$saveName"_jetpT20_30"
  fi
  if [[ "$FILE" = *"_centrality0_10_"* ]];then
    saveName=$saveName"_Cent0_10"
  elif [[ "$FILE" = *"_centrality0_100_"* ]];then
    saveName=$saveName"_Cent0_100"
  elif [[ "$FILE" = *"_centrality60_100_"* ]];then
    saveName=$saveName"_Cent60_100"
  fi   
  saveName=$saveName.pdf
  
  aliroot GetAndSaveCanvas.C"(\"$FILE\", \"$saveName\")" -l -b -q
done
