# BPHParkingAnalysis
Analysis package for analysis of BPH parked data

### Instructions for 10_2_5
```
cmsrel CMSSW_10_2_5
cd CMSSW_10_2_5/src/
cmsenv

git clone https://github.com/ICBPHCMS/BPHParkingAnalysis.git
scram b
```

### Producing ntuples
```
cd BPHParkingAnalysis/NtupleProducer/bin
BToKstllNtupleProducer --isMC (0,1,2) --isResonant (0, 1, -1) --isEleFS (0, 1) --isKstFS (0, 1) --isLT (0, 1) --output ("outfile") --input ("inputFile")

isMC == 2 is for (old) Thomas'MC; isMC == 1 is new; isMC == 0 is for DATA
isEleFS == 1 is for electron final state; isEleFS == 0 is for muon final state
isKstFS == 1 is for Kst channel
isLT == 1 is for lepton-track-track configuration; isLT == 0 is for lepton-lepton-track configuration
```

### Analyzer
```
cd BPHParkingAnalysis/NtupleProducer/macro
g++ -Wall -o analyzeCharged_fastDATA_Kstll `root-config --cflags --glibs` -lRooFitCore analyzeCharged_fastDATA_Kstll.cpp


Commands to process your ntuples in parallel with the analyzer:

cd ../scripts

python cmsSplit.py --anType scriptAndJOBID --cfg config_runAnalyis.sh --tag YOUR-TAG --listquery -i PATH-TO-YOUR-NTUPLE-LIST.txt --filesperjob 10 --storeArea PATH-TO-YOUR-STORE-AREA

source launch_YOUR-TAG.sh
```

### Fit and plots
```
cd BPHParkingAnalysis/NtupleProducer/macro

Ele final state: root fitBmass_fromHistos.C'(1, "PATH-TO-YOUR-ANALYZER-OUTPUT-FILE")'
Mu final state:  root fitBmass_fromHistos.C'(0, "PATH-TO-YOUR-ANALYZER-OUTPUT-FILE")'
```

