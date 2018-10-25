##with input a DASpath in case of cmsRun 
python ~/cmsSplit.py --anType cmsRun --cfg test94X_NANO_template.py  --tag BToKJPsiee  --dasquery  --das /Bu_KJPsi_ee_Pythia/tstreble-BuToKJPsiee_Pythia_MINIAODSIM_18_06_05-393c5b9eb90ffa47a3c9a0f5562ec979/USER --filesperjob 2  --storeArea /vols/cms/amartell/BParking/nanoPROD


##with input txt in case of simple script
python ~/cmsSplit.py --anType script --cfg config_data_runNtuples.sh --tag BToKee_2018B_BPH1 --listquery -i list_data_BToKee_2018B_BPH1_PR22.txt --filesperjob 1 --storeArea /vols/cms/amartell/BParking/ntuPROD/BPHrun2018B/


##with input txt in case of simple script and JOBID needed as input
python ~/cmsSplit.py --anType scriptAndJOBID --cfg config_data_runB.sh --tag BToKee_2018B_BPH1 --listquery -i list_data_BToKee_2018B_BPH1_PR22.txt --filesperjob 1 --storeArea /vols/cms/amartell/BParking/ntuPROD/BPHrun2018B/
