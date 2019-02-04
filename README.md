# BPHParkingAnalysis
Analysis package for analysis of BPH parked data

### Instructions for 10_2_5
```
cmsrel CMSSW_10_2_5
cd CMSSW_10_1_5/src/
cmsenv

git clone https://github.com/ICBPHCMS/BPHParkingAnalysis.git
scram b
```

### Producing ntuples
```
cd BPHParkingAnalysis/NtupleProducer/bin

BToKstllNtupleProducer --isMC (0,1,2) --isResonant (0, 1, -1) --isEleFS (0, 1) --isKstFS (0, 1) --isLT (0, 1) --output ("outfile") --input ("inputFile")
```
