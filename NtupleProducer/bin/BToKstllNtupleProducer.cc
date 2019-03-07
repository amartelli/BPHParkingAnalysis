// BToKstllNtupleProducer --isMC (0,1,2) --isResonant (0, 1, -1) --isEleFS (0, 1) --isKstFS (0, 1) --isLT (0, 1) --output ("outfile") --input ("inputFile")

// new features => gen_index refers to the closest in dR - provided gen B and decays are within acceptance -, 
//                 no minimum dR required
//                 saved branch with dR of gen-reco lep1, lep2, kaon
//              => for gen_tag_muon do not require minimum pT
//                 saved branch with highest pT of extra muon (tag candidate)
//              => matching for reco tag muon is dR < 0.01
//              => isMC == 2 is for old - Thomas' - MC; isMC == 1 is new; isMC == 0 is DATA 
// 

#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TLorentzVector.h"

#include "NanoAODTree.h"


float JPsiMass_ = 3.0969;
float BuMass_ = 5.279;
float KaonMass_ = 0.493677;
float MuonMass_ = 0.10565837;
float ElectronMass_ = 0.5109989e-3;


float maxEtacceptance_ = 2.45; //2.4
float minPtacceptance_ = 0.7;  //1.

bool comparePairs(const std::pair<int, float>& i, const std::pair<int, float>& j){
  return i.second > j.second;
}


int main(int argc, char **argv){

  if(argc < 2) {
    std::cout << " Missing arguments " << std::endl;
    return -1;
  }
  int isMC = -1;
  int isResonant = -1;
  int isEleFinalState = -1;
  int isKstFinalState = -1;
  int isLeptonTrack = -1;
  string output = "";
  string input = "";
  string listFilesTXT = "";
  bool inputTXT = false;
  float genMuPtCut = 5.;
  bool overwrite = false;
  for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "--isMC") {
      if (i + 1 < argc) {
	isMC = atoi(argv[i+1]);
	break;
      } else {
	std::cerr << " --isMC option requires one argument " << std::endl;
	return 1;
      }
    }
  }for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "--isResonant") {
      if (i + 1 < argc) {
	isResonant = atoi(argv[i+1]);
	break;
      } else {
	std::cerr << " --isResonant option requires one argument " << std::endl;
	return 1;
      }
    }
  }for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "--isEleFS") {
      if (i + 1 < argc) {
	isEleFinalState = atoi(argv[i+1]);
	break;
      } else {
	std::cerr << " --isElFinalState option requires one argument " << std::endl;
	return 1;
      }
    }
  }for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "--isKstFS") {
      if (i + 1 < argc) {
	isKstFinalState = atoi(argv[i+1]);
	break;
      } else {
	std::cerr << " --isKstFinalState option requires one argument " << std::endl;
	return 1;
      }
    }
  }for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "--isLT") {
      if (i + 1 < argc) {
	isLeptonTrack = atoi(argv[i+1]);
	break;
      } else {
	std::cerr << " --isLeptonTrack option requires one argument " << std::endl;
	return 1;
      }
    }
  }for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "--output") {
      if (i + 1 < argc) {
	output = argv[i+1];
	break;
      } else {
	std::cerr << "--output option requires one argument." << std::endl;
	return 1;
      }      
    }  
  }for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "--input") {
      if (i + 1 < argc) {
        input = argv[i+1];
        break;
      } else {
	std::cerr << "--intput option requires one argument." << std::endl;
        return 1;
      }
    }
    else if(std::string(argv[i]) == "--inputTXT") {
      if (i + 1 < argc) {
        inputTXT = true;
        listFilesTXT = argv[i+1];
        break;
      } else {
	std::cerr << "--intputTXT option requires one argument." << std::endl;
        return 1;
      }
    }
  }

  if(output == ""){
    std::cerr << "--output argument required" << std::endl;
    return 1;
  }
  if(listFilesTXT == "" && input == ""){
    std::cerr << "--input argument required" << std::endl;
    return 1;
  }

  /*
  bool overwrite = false;
  for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "--overwrite") {
      overwrite = true;
      break;
    }
  }
  */



  //Always saveFullNanoAOD info because we are skimming
  bool saveFullNanoAOD = true;
  /*for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "--saveFullNanoAOD") {
      saveFullNanoAOD = true;
      break;
    }
    }*/
 


  if(isKstFinalState){
    std::cout << " implementation missing " << std::endl;
    return -10;
  }

  TFile* f_new = TFile::Open(output.c_str());
  if(f_new!=0 && !overwrite){
    cout<<output<<" already exists, please delete it before converting again"<<endl;
    return 0;
  }
  f_new = TFile::Open(output.c_str(),"RECREATE");


  TChain* oldtree = new TChain("Events");
  if(inputTXT){
    string reader;
    std::ifstream inFileLong;
    inFileLong.open(listFilesTXT.c_str(), std::ios::in);

    while(!inFileLong.eof()){
      inFileLong >> reader;
      if( inFileLong.eof() ) break;
      std::cout << " Adding " << reader << std::endl;
      oldtree->Add(reader.c_str());
    }
  }
  else   oldtree->Add(input.c_str());
  NanoAODTree* tree = new NanoAODTree(oldtree);

  TTree* tree_new=new TTree("BToKstllTree","BToKstllTree");
  if(saveFullNanoAOD)
    tree_new=tree->GetTree()->CloneTree(0);


  //New branches
  int _BToKstll_sel_index = -1;
  std::vector<int> _BToKstll_order_index;
  int _Muon_sel_index = -1; //Probe muon with selection algo.
  std::vector<int> _Muon_tag_index; //Probe muon with selection algo.
  int _Muon_probe_index = -1; //Probe muon for Acc.xEff. = _Muon_sel_index in data

  tree_new->Branch("BToKstll_sel_index",&_BToKstll_sel_index,"BToKstll_sel_index/I");
  tree_new->Branch("BToKstll_order_index", &_BToKstll_order_index);
  tree_new->Branch("Muon_sel_index",&_Muon_sel_index,"Muon_sel_index/I");
  tree_new->Branch("Muon_tag_index",&_Muon_tag_index);
  tree_new->Branch("Muon_probe_index",&_Muon_probe_index,"Muon_probe_index/I");


  int _GenPart_BToKstll_index = -1;
  int _GenPart_JPsiFromB_index = -1;
  int _GenPart_KFromB_index = -1;
  int _GenPart_lep1FromB_index = -1;
  int _GenPart_lep2FromB_index = -1;
  int _BToKstll_gen_index = -1;  //index of reco closest to gen (each daughter within 0.1)
  float _BToKstll_gendR_lep1FromB;
  float _BToKstll_gendR_lep2FromB;
  float _BToKstll_gendR_KFromB;
  float _BToKstll_gen_llMass = -1;
  float _BToKstll_gen_mass = -1;
  float _BToKstll_gen_muonTag_hpT = -1;

  //the 3 following are the index of the reco closest to the gen (dR < 0.1)
  //only for no LT
  int _Lep1_gen_index = -1;
  int _Lep2_gen_index = -1;
  int _KPFCand_gen_index = -1;  

  if(isMC != 0){
    tree_new->Branch("GenPart_BToKstll_index",&_GenPart_BToKstll_index,"GenPart_BToKstll_index/I");
    tree_new->Branch("GenPart_JPsiFromB_index",&_GenPart_JPsiFromB_index,"GenPart_JPsiFromB_index/I");
    tree_new->Branch("GenPart_KFromB_index",&_GenPart_KFromB_index,"GenPart_KFromB_index/I");
    tree_new->Branch("GenPart_lep1FromB_index",&_GenPart_lep1FromB_index,"GenPart_lep1FromB_index/I");
    tree_new->Branch("GenPart_lep2FromB_index",&_GenPart_lep2FromB_index,"GenPart_lep2FromB_index/I");
    tree_new->Branch("BToKstll_gen_index",&_BToKstll_gen_index,"BToKstll_gen_index/I");
    tree_new->Branch("BToKstll_gendR_lep1FromB",&_BToKstll_gendR_lep1FromB,"BToKstll_gendR_lep1FromB/F");
    tree_new->Branch("BToKstll_gendR_lep2FromB",&_BToKstll_gendR_lep2FromB,"BToKstll_gendR_lep2FromB/F");
    tree_new->Branch("BToKstll_gendR_KFromB",&_BToKstll_gendR_KFromB,"BToKstll_gendR_KFromB/F");
    tree_new->Branch("BToKstll_gen_llMass",&_BToKstll_gen_llMass,"BToKstll_gen_llMass/F");
    tree_new->Branch("BToKstll_gen_mass",&_BToKstll_gen_mass,"BToKstll_gen_mass/F");
    tree_new->Branch("BToKstll_gen_muonTag_hpT",&_BToKstll_gen_muonTag_hpT,"BToKstll_gen_muonTag_hpT/F");
    tree_new->Branch("Lep1_gen_index",&_Lep1_gen_index,"Lep1_gen_index/I");
    tree_new->Branch("Lep2_gen_index",&_Lep2_gen_index,"Lep2_gen_index/I");
    tree_new->Branch("KPFCand_gen_index",&_KPFCand_gen_index,"KPFCand_gen_index/I");
  }


  bool _HLT_Mu8p5_IP3p5 = false;
  bool _HLT_Mu10p5_IP3p5 = false;
  bool _HLT_Mu9_IP6 = false;
  bool _HLT_Mu8_IP3 = false;
  bool _HLT_BPHParking = false;

  //for MC
  bool _HLT_Mu8_IP6 = false;
  bool _HLT_Mu8_IP5 = false;
  bool _HLT_Mu9_IP4 = false;
  bool _HLT_Mu7_IP4 = false;
  bool _HLT_Mu9_IP5 = false;
  bool _HLT_Mu12_IP6 = false;

  bool _Muon_isHLT_BPHParking[kMuonMax];

  tree_new->Branch("HLT_Mu8p5_IP3p5",&_HLT_Mu8p5_IP3p5,"HLT_Mu8p5_IP3p5/O");
  tree_new->Branch("HLT_Mu10p5_IP3p5",&_HLT_Mu10p5_IP3p5,"HLT_Mu10p5_IP3p5/O");
  tree_new->Branch("HLT_Mu9_IP6",&_HLT_Mu9_IP6,"HLT_Mu9_IP6/O");
  tree_new->Branch("HLT_Mu8_IP3",&_HLT_Mu8_IP3,"HLT_Mu8_IP3/O");
  tree_new->Branch("HLT_BPHParking",&_HLT_BPHParking,"HLT_BPHParking/O");
  
  tree_new->Branch("Muon_isHLT_BPHParking",_Muon_isHLT_BPHParking,"Muon_isHLT_BPHParking[nMuon]/O");
  

  int nentries = tree->GetEntries();
  std::cout << " Nentries = " << nentries << std::endl;
  std::cout << " isMC = " << isMC << std::endl;
  
  for (int iEntry = 0; iEntry < nentries; ++iEntry){
    
    bool debug = false;
    //if(iEntry == 419) debug = true;

    int out = tree->GetEntry(iEntry);
    if(out<0){
      std::cout << " Error retrievieng entry #" << iEntry << std::endl;
      return -1;
    }

    if(iEntry%10000==0) std::cout << " Entry #" << iEntry << " " << int(100*float(iEntry)/nentries) << "%" << std::endl;

    _BToKstll_sel_index = -1;
    _BToKstll_order_index.clear();
    _Muon_sel_index = -1;  //tag muon
    _Muon_tag_index.clear();
    _Muon_probe_index = -1;

    _GenPart_BToKstll_index = -1;
    _GenPart_JPsiFromB_index = -1;
    _GenPart_KFromB_index = -1;
    _GenPart_lep1FromB_index = -1;
    _GenPart_lep2FromB_index = -1;
    _BToKstll_gen_index = -1; 
    _BToKstll_gendR_lep1FromB = -1; 
    _BToKstll_gendR_lep2FromB = -1; 
    _BToKstll_gendR_KFromB = -1; 
    _BToKstll_gen_llMass = -1;
    _BToKstll_gen_mass = -1;
    _BToKstll_gen_muonTag_hpT = -1;
    _Lep1_gen_index = -1;
    _Lep2_gen_index = -1;
    _KPFCand_gen_index = -1;

    _HLT_Mu8p5_IP3p5 = false;
    _HLT_Mu10p5_IP3p5 = false;
    _HLT_Mu9_IP6 = false;
    _HLT_Mu8_IP3 = false;
    _HLT_BPHParking = false;
    
    //for MC
    _HLT_Mu8_IP6 = false;
    _HLT_Mu8_IP5 = false;
    _HLT_Mu9_IP4 = false;
    _HLT_Mu7_IP4 = false;
    _HLT_Mu9_IP5 = false;
    _HLT_Mu12_IP6 = false;

    //look for trigger muon
    int nMuon = tree->nMuon;
    int nTriggeringMuons = 0;
    for(int i_mu=0; i_mu<nMuon; i_mu++){


      //Trigger selection + matching
      _HLT_Mu8p5_IP3p5 = tree->HLT_Mu8p5_IP3p5_part0
	  || tree->HLT_Mu8p5_IP3p5_part1
	  || tree->HLT_Mu8p5_IP3p5_part2
	  || tree->HLT_Mu8p5_IP3p5_part3
	  || tree->HLT_Mu8p5_IP3p5_part4
	  || tree->HLT_Mu8p5_IP3p5_part5;
	_HLT_Mu10p5_IP3p5 = tree->HLT_Mu10p5_IP3p5_part0
	  || tree->HLT_Mu10p5_IP3p5_part1
	  || tree->HLT_Mu10p5_IP3p5_part2
	  || tree->HLT_Mu10p5_IP3p5_part3
	  || tree->HLT_Mu10p5_IP3p5_part4
	  || tree->HLT_Mu10p5_IP3p5_part5;
	_HLT_Mu9_IP6 = tree->HLT_Mu9_IP6_part0
	  || tree->HLT_Mu9_IP6_part1
	  || tree->HLT_Mu9_IP6_part2
	  || tree->HLT_Mu9_IP6_part3
	  || tree->HLT_Mu9_IP6_part4
	  || tree->HLT_Mu9_IP6_part5;
	_HLT_Mu8_IP3 = tree->HLT_Mu8_IP3_part0
	  || tree->HLT_Mu8_IP3_part1
	  || tree->HLT_Mu8_IP3_part2
	  || tree->HLT_Mu8_IP3_part3
	  || tree->HLT_Mu8_IP3_part4
	  || tree->HLT_Mu8_IP3_part5;
	if(isMC != 0){
	  _HLT_Mu8_IP6 = tree->HLT_Mu8_IP6_part0
	    || tree->HLT_Mu8_IP6_part1
	    || tree->HLT_Mu8_IP6_part2
	    || tree->HLT_Mu8_IP6_part3
	    || tree->HLT_Mu8_IP6_part4
	    || tree->HLT_Mu8_IP6_part5;
	  _HLT_Mu8_IP5 = tree->HLT_Mu8_IP5_part0
	    || tree->HLT_Mu8_IP5_part1
	    || tree->HLT_Mu8_IP5_part2
	    || tree->HLT_Mu8_IP5_part3
	    || tree->HLT_Mu8_IP5_part4
	    || tree->HLT_Mu8_IP5_part5;
	  _HLT_Mu9_IP4 = tree->HLT_Mu9_IP4_part0
	    || tree->HLT_Mu9_IP4_part1
	    || tree->HLT_Mu9_IP4_part2
	    || tree->HLT_Mu9_IP4_part3
	    || tree->HLT_Mu9_IP4_part4
	    || tree->HLT_Mu9_IP4_part5;
	  _HLT_Mu7_IP4 = tree->HLT_Mu7_IP4_part0
	    || tree->HLT_Mu7_IP4_part1
	    || tree->HLT_Mu7_IP4_part2
	    || tree->HLT_Mu7_IP4_part3
	    || tree->HLT_Mu7_IP4_part4
	    || tree->HLT_Mu7_IP4_part5;
	  _HLT_Mu9_IP5 = tree->HLT_Mu9_IP5_part0
	    || tree->HLT_Mu9_IP5_part1
	    || tree->HLT_Mu9_IP5_part2
	    || tree->HLT_Mu9_IP5_part3
	    || tree->HLT_Mu9_IP5_part4
	    || tree->HLT_Mu9_IP5_part5;
	  _HLT_Mu12_IP6 = tree->HLT_Mu12_IP6_part0
	    || tree->HLT_Mu12_IP6_part1
	    || tree->HLT_Mu12_IP6_part2
	    || tree->HLT_Mu12_IP6_part3
	    || tree->HLT_Mu12_IP6_part4
	    || tree->HLT_Mu12_IP6_part5;
	}
	_HLT_BPHParking = _HLT_Mu8p5_IP3p5 || _HLT_Mu10p5_IP3p5 || _HLT_Mu9_IP6 || _HLT_Mu8_IP3;
	if(isMC != 0) 	_HLT_BPHParking = _HLT_Mu8p5_IP3p5 || _HLT_Mu10p5_IP3p5 || _HLT_Mu9_IP6 || _HLT_Mu8_IP3 || 
			 _HLT_Mu8_IP6 || _HLT_Mu8_IP5 || _HLT_Mu9_IP4 || _HLT_Mu7_IP4 || _HLT_Mu9_IP5 || _HLT_Mu12_IP6 ;

	TLorentzVector mu;
	//mu.SetPtEtaPhiM(tree->Muon_pt[i_mu],tree->Muon_eta[i_mu],tree->Muon_phi[i_mu],tree->Muon_mass[i_mu]);
	mu.SetPtEtaPhiM(tree->Muon_pt[i_mu],tree->Muon_eta[i_mu],tree->Muon_phi[i_mu],MuonMass_);

	bool isTrigMatched = false;
	int nTrigObj = tree->nTrigObj;
	for(int i_trig = 0; i_trig<nTrigObj; i_trig++){
	  if( tree->TrigObj_id[i_trig]==13 && ((tree->TrigObj_filterBits[i_trig])>>3)&1 ){
	    TLorentzVector trig;
	    trig.SetPtEtaPhiM(tree->TrigObj_pt[i_trig],tree->TrigObj_eta[i_trig],tree->TrigObj_phi[i_trig], MuonMass_);
	    float dR = mu.DeltaR(trig);
	    if(dR < 0.01){
	      isTrigMatched = true;
	      break;
	    }
	  }
	}
	if(isTrigMatched && debug) {
          ++nTriggeringMuons;
	  std::cout << " mu HLT matched = " << i_mu << std::endl;
        }
	_Muon_isHLT_BPHParking[i_mu] = isTrigMatched;
    }

    if(debug){
    int nRecoMuons = 0;
    int nRecoMuonsHLTMatched = 0;
    ///check reco muon collection
    for(int i_mu=0; i_mu<nMuon; i_mu++){
      //muon passing soft ID + trigger matched                                                                                                                          
      ++nRecoMuons;
      std::cout << " tree->Muon_pt[i_mu] = " << tree->Muon_pt[i_mu] << " tree->Muon_charge[i_mu] = " << tree->Muon_charge[i_mu] 
		<< " tree->Muon_eta[i_mu] = " << tree->Muon_eta[i_mu] << " tree->Muon_phi[i_mu] = " << tree->Muon_phi[i_mu] << std::endl;     
      if(!tree->Muon_softId[i_mu] || tree->Muon_pt[i_mu] < genMuPtCut) continue;
      if(isMC != 2 && !_Muon_isHLT_BPHParking[i_mu]) continue;
      ++nRecoMuonsHLTMatched;
    }
    std::cout << " nTriggeringMuons = " << nTriggeringMuons << " nRecoMuons = " << nRecoMuons << " nRecoMuonsHLTMatched = " << nRecoMuonsHLTMatched << std::endl;
    }

    //Select the BToKll candidate with reco criteria

    int nBinTree = tree->nBToKstll;
    float best_B_CL_vtx = -1.;
    std::vector<std::pair<int, float>> B_vtxCL_idx_val;
    _Muon_tag_index.resize(nBinTree);


    if(debug) std::cout << " >>> nBinTree = " << nBinTree << std::endl;
    int idxL1 = -1;
    int idxL2 = -1;
    for(int i_Btree=0; i_Btree<nBinTree; ++i_Btree){            

      //already applied in nanoAOD production
      if(tree->BToKstll_lep1_charge[i_Btree] * tree->BToKstll_lep2_charge[i_Btree] > 0.) continue;

      //consider triplet only if other trigger muon exists
      _Muon_sel_index = -1;
      
      //avoid rechecking extra muon match if same l1 and l2 as before
      if(i_Btree > 0 && isLeptonTrack == 0 &&
	 (idxL1 == tree->BToKstll_lep1_index[i_Btree] && idxL2 == tree->BToKstll_lep2_index[i_Btree])){
	if(debug) std::cout << " looping i_Btree = " << i_Btree << " idxL1 = " << idxL1 << " idxL2 = " << idxL2 
			    << " _Muon_tag_index[i_Btree-1] = " << _Muon_tag_index[i_Btree-1] << std::endl;
	_Muon_tag_index[i_Btree] = _Muon_tag_index[i_Btree-1];
	// idxL1 = tree->BToKstll_lep1_index[i_Btree];
	// idxL2 = tree->BToKstll_lep2_index[i_Btree];
	_Muon_sel_index = _Muon_tag_index[i_Btree-1];
      }
      else{
      for(int i_mu=0; i_mu<nMuon; i_mu++){

	//to avoid rechecking the exact same lepton lepton combination
	if(isLeptonTrack == 0){
	  idxL1 = tree->BToKstll_lep1_index[i_Btree];
	  idxL2 = tree->BToKstll_lep2_index[i_Btree];
	}

	if(debug) std::cout << " else looping i_Btree = " << i_Btree << " idxL1 = " << idxL1 << " idxL2 = " << idxL2 << std::endl;

	//muon passing soft ID + trigger matched
	if(!tree->Muon_softId[i_mu] || tree->Muon_pt[i_mu] < genMuPtCut) continue;
	if(isMC != 2 && !_Muon_isHLT_BPHParking[i_mu]) continue;


	//if lepton lepton just check index
	if(isLeptonTrack == 0 && isEleFinalState == 0 &&
	   (i_mu == tree->BToKstll_lep1_index[i_Btree] || i_mu == tree->BToKstll_lep2_index[i_Btree])) continue;
	


	//if tracks or electrons or is passed
	//check dR for leading and subleading
	
	TLorentzVector lep1_tlv;
	TLorentzVector lep2_tlv;
	TLorentzVector muHLT_tlv;
	
	lep1_tlv.SetPtEtaPhiM(tree->BToKstll_lep1_pt[i_Btree],
			      tree->BToKstll_lep1_eta[i_Btree],
			      tree->BToKstll_lep1_phi[i_Btree],
			    (isEleFinalState == 1) ? ElectronMass_ : MuonMass_);
	lep2_tlv.SetPtEtaPhiM(tree->BToKstll_lep2_pt[i_Btree],
			      tree->BToKstll_lep2_eta[i_Btree],
			      tree->BToKstll_lep2_phi[i_Btree],
			      (isEleFinalState == 1) ? ElectronMass_ : MuonMass_);
	muHLT_tlv.SetPtEtaPhiM(tree->Muon_pt[i_mu],
			       tree->Muon_eta[i_mu],
			       tree->Muon_phi[i_mu],
			       MuonMass_);
	
	float dR_lep1FromHLT = lep1_tlv.DeltaR(muHLT_tlv);
	float dR_lep2FromHLT = lep2_tlv.DeltaR(muHLT_tlv);

	if(dR_lep1FromHLT < 0.02) continue;
	if(dR_lep2FromHLT < 0.02) continue;

	if(debug) std::cout << " >>> l1idx = " << tree->BToKstll_lep1_index[i_Btree]
                            << " >>> l2idx = " << tree->BToKstll_lep2_index[i_Btree]
                            << " muIdx = " << i_mu
                            << std::endl;

	//enough to find 1 extra muon matched to the trigger    
	_Muon_sel_index = i_mu;
	_Muon_tag_index[i_Btree] = i_mu;

	break;
	
      }//loop over muons
      }
      /////////////////////
      //not found other reco muon matched to trigger for this triplet
      if(_Muon_sel_index == -1) {
	_Muon_tag_index[i_Btree] = -1;
      }


      float B_CL_vtx = tree->BToKstll_B_CL_vtx[i_Btree];
      B_vtxCL_idx_val.push_back(std::pair<int, float>(i_Btree, B_CL_vtx));
      
      if( best_B_CL_vtx < 0. || B_CL_vtx>best_B_CL_vtx ){
	best_B_CL_vtx = B_CL_vtx;
	_BToKstll_sel_index = i_Btree;
      }
    }//loop over triplets

    std::sort(B_vtxCL_idx_val.begin(), B_vtxCL_idx_val.end(), comparePairs);
    int ipCount = 0;
    _BToKstll_order_index.resize(B_vtxCL_idx_val.size());
    for(auto ip : B_vtxCL_idx_val){
      _BToKstll_order_index[ip.first] = ipCount;
      ++ipCount;
    }

    if(debug) std::cout << " >>> _BToKstll_sel_index = " << _BToKstll_sel_index << std::endl;


    //Take as tag muon leading soft ID muon + trigger-matched
    if(isMC == 0 && _BToKstll_sel_index<0) continue;

    //re-assign proper Muon_sel_index for chosen triplet
    if(_BToKstll_sel_index>=0){
      _Muon_sel_index = _Muon_tag_index[_BToKstll_sel_index]; 
    }
    if(debug)    std::cout << " iEntry = " << iEntry << " _Muon_sel_index = " << _Muon_sel_index << std::endl;


    //!!! Can have selected B->Kstll even without additional tag muons
    //Needs to ask both BToKstll_sel_index>=0 && Muon_sel_index>=0 for analysis on probe side
    //BToKstll_sel_index not reset to -1 for potential analysis of the probe side
    

    //Select the BToKstll candidate based on gen matching

    if(isMC != 0){
      int nGenPart = tree->nGenPart;

      int leptonID = (isEleFinalState == 1) ? 11 : 13;

      if(isResonant){
	for(int i_Bu=0; i_Bu<nGenPart; i_Bu++){

	  _GenPart_JPsiFromB_index = -1;
	  _GenPart_lep1FromB_index = -1;
	  _GenPart_lep2FromB_index = -1;
	  _GenPart_KFromB_index = -1;

	  if(abs(tree->GenPart_pdgId[i_Bu])==521){
	    for(int i_gen=0; i_gen<nGenPart; i_gen++){

	      int pdgId = tree->GenPart_pdgId[i_gen];
	      int mother_index = tree->GenPart_genPartIdxMother[i_gen];

	      int lep1_index = -1;
	      int lep2_index = -1;

	      if(abs(pdgId)==443 && mother_index == i_Bu){

		for(int j_gen=0; j_gen<nGenPart; j_gen++){
		  int pdgId = tree->GenPart_pdgId[j_gen];
		  float partPt = tree->GenPart_pt[i_gen];
		  float partEta = tree->GenPart_eta[i_gen];
		  if(partPt < minPtacceptance_) continue;
		  if(std::abs(partEta) > maxEtacceptance_) continue;
		  int mother_index = tree->GenPart_genPartIdxMother[j_gen];
		  if(abs(pdgId)==leptonID && mother_index == i_gen && lep1_index < 0)
		    lep1_index = j_gen;
		  else if(abs(pdgId)==leptonID && mother_index == i_gen)
		    lep2_index = j_gen;
		  if(lep1_index >= 0 && lep2_index >= 0) break;
		}

		if(lep1_index >= 0 && lep2_index >= 0){
		  _GenPart_JPsiFromB_index = i_gen;
		  _GenPart_lep1FromB_index = lep1_index;
		  _GenPart_lep2FromB_index = lep2_index;
		}
		else break;

	      }

	      else if(abs(pdgId)==321 && mother_index == i_Bu) _GenPart_KFromB_index = i_gen;

	      else if(mother_index == i_Bu) break; //Additional B decay products
	    }
	  }

	  if(_GenPart_JPsiFromB_index >= 0 && _GenPart_KFromB_index >=0){
	    _GenPart_BToKstll_index = i_Bu;
	    break;
	  }
	}	
      }//resonant
      else{ 

	for(int i_Bu=0; i_Bu<nGenPart; i_Bu++){
	  _GenPart_JPsiFromB_index = -1;
	  _GenPart_lep1FromB_index = -1;
	  _GenPart_lep2FromB_index = -1;
	  _GenPart_KFromB_index = -1;

	  if(abs(tree->GenPart_pdgId[i_Bu]) == 521){

	    for(int i_gen=0; i_gen<nGenPart; i_gen++){

	      int pdgId = tree->GenPart_pdgId[i_gen];
	      float partPt = tree->GenPart_pt[i_gen];
	      float partEta = tree->GenPart_eta[i_gen];
	      if(partPt < minPtacceptance_) continue;
	      if(std::abs(partEta) > maxEtacceptance_) continue;
	      int mother_index = tree->GenPart_genPartIdxMother[i_gen];
	      if(abs(pdgId) == leptonID && mother_index == i_Bu && _GenPart_lep1FromB_index < 0)
		_GenPart_lep1FromB_index = i_gen;
	      else if(abs(pdgId)==leptonID && mother_index == i_Bu)
		_GenPart_lep2FromB_index = i_gen;
	      else if(abs(pdgId)==321 && mother_index == i_Bu)
		_GenPart_KFromB_index = i_gen;
	      else if(mother_index == i_Bu) break;
	    }
	  }//if B

	  if(_GenPart_lep1FromB_index >= 0 && _GenPart_lep2FromB_index >= 0 && _GenPart_KFromB_index >= 0){
	    _GenPart_BToKstll_index = i_Bu;
	    break;
	  }
	} //loop over B
      }// non resonant
    
      if(_GenPart_BToKstll_index >= 0){

	//lep1FromB stored a leading daughter
	if(tree->GenPart_pt[_GenPart_lep2FromB_index] > tree->GenPart_pt[_GenPart_lep1FromB_index]){
	  int i_temp = _GenPart_lep1FromB_index;
	  _GenPart_lep1FromB_index = _GenPart_lep2FromB_index;
	  _GenPart_lep2FromB_index = i_temp;
	}

	TLorentzVector gen_KFromB_tlv;
	TLorentzVector gen_lep1FromB_tlv;
	TLorentzVector gen_lep2FromB_tlv;
	
	gen_KFromB_tlv.SetPtEtaPhiM(tree->GenPart_pt[_GenPart_KFromB_index],
				    tree->GenPart_eta[_GenPart_KFromB_index],
				    tree->GenPart_phi[_GenPart_KFromB_index],
				    KaonMass_);
	gen_lep1FromB_tlv.SetPtEtaPhiM(tree->GenPart_pt[_GenPart_lep1FromB_index],
				      tree->GenPart_eta[_GenPart_lep1FromB_index],
				      tree->GenPart_phi[_GenPart_lep1FromB_index],
				       (isEleFinalState == 1) ? ElectronMass_ : MuonMass_);
	gen_lep2FromB_tlv.SetPtEtaPhiM(tree->GenPart_pt[_GenPart_lep2FromB_index],
				       tree->GenPart_eta[_GenPart_lep2FromB_index],
				       tree->GenPart_phi[_GenPart_lep2FromB_index],
				       (isEleFinalState == 1) ? ElectronMass_ : MuonMass_);

	_BToKstll_gen_llMass = (gen_lep1FromB_tlv+gen_lep2FromB_tlv).Mag();
	_BToKstll_gen_mass = (gen_lep1FromB_tlv+gen_lep2FromB_tlv+gen_KFromB_tlv).Mag();

	float best_dR = -1.;
	
	for(int i_Btree=0; i_Btree<nBinTree; ++i_Btree){

	  TLorentzVector kaon_tlv;
	  TLorentzVector lep1_tlv;
	  TLorentzVector lep2_tlv;

	  kaon_tlv.SetPtEtaPhiM(tree->BToKstll_kaon_pt[i_Btree],
				tree->BToKstll_kaon_eta[i_Btree],
				tree->BToKstll_kaon_phi[i_Btree],
				KaonMass_);
	  lep1_tlv.SetPtEtaPhiM(tree->BToKstll_lep1_pt[i_Btree],
				tree->BToKstll_lep1_eta[i_Btree],
				tree->BToKstll_lep1_phi[i_Btree],
				(isEleFinalState == 1) ? ElectronMass_ : MuonMass_);
	  lep2_tlv.SetPtEtaPhiM(tree->BToKstll_lep2_pt[i_Btree],
				tree->BToKstll_lep2_eta[i_Btree],
				tree->BToKstll_lep2_phi[i_Btree],
				(isEleFinalState == 1) ? ElectronMass_ : MuonMass_);
	  
	  float dR_KFromB = kaon_tlv.DeltaR(gen_KFromB_tlv);
	  float dR_lep1FromB = min(lep1_tlv.DeltaR(gen_lep1FromB_tlv), lep2_tlv.DeltaR(gen_lep1FromB_tlv));
	  float dR_lep2FromB = min(lep1_tlv.DeltaR(gen_lep2FromB_tlv), lep2_tlv.DeltaR(gen_lep2FromB_tlv));
	  //Should check that same objects not selected twice

	  float dR_tot = dR_KFromB + dR_lep1FromB + dR_lep2FromB; //In case several BToKmumu matches, take the closest one in dR_tot

	  //if( dR_lep1FromB <0.1 && dR_lep2FromB <0.1
	  //  && 
	  if(best_dR <0. || dR_tot < best_dR){
	    best_dR = dR_tot;
	    _BToKstll_gen_index = i_Btree;
	    _BToKstll_gendR_lep1FromB = dR_lep1FromB;
	    _BToKstll_gendR_lep2FromB = dR_lep2FromB;
	    _BToKstll_gendR_KFromB = dR_KFromB;
	  }
	}

	/*
	float best_dR_lep1FromB = -1.;
	float best_dR_lep2FromB = -1.;
	float best_dR_KFromB = -1.;
	*/

	//for the moment just check best reco closest to gen
	//for lepton-lepton combination only
	/*
	if(!isLeptonTrack){  
	for(int i_mu=0; i_mu<nMuon; i_mu++){

	  TLorentzVector mu_tlv;
	  mu_tlv.SetPtEtaPhiM(tree->Muon_pt[i_mu],
			      tree->Muon_eta[i_mu],
			      tree->Muon_phi[i_mu],
			      MuonMass_);

	  float dR_mu1FromB = mu_tlv.DeltaR(gen_mu1FromB_tlv);
	  float dR_mu2FromB = mu_tlv.DeltaR(gen_mu2FromB_tlv);

	  if(dR_mu1FromB<0.1 && (best_dR_mu1FromB<0. || dR_mu1FromB<best_dR_mu1FromB)){
	    _Muon_mu1FromB_index = i_mu;
	    best_dR_mu1FromB = dR_mu1FromB;
	  }
	  if(dR_mu2FromB<0.1 && (best_dR_mu2FromB<0. || dR_mu2FromB<best_dR_mu2FromB)){
	    _Muon_mu2FromB_index = i_mu;
	    best_dR_mu2FromB = dR_mu2FromB;
	  }
	}
	}//end of lepton-lepton comb

	int nPFCand = tree->nPFCand;
	for(int i_pf=0;i_pf<nPFCand;i_pf++){

	  if(abs(tree->PFCand_pdgId[i_pf])!=211) continue;

	  TLorentzVector PFCand_tlv;
	  PFCand_tlv.SetPtEtaPhiM(tree->PFCand_pt[i_pf],tree->PFCand_eta[i_pf],tree->PFCand_phi[i_pf],tree->PFCand_mass[i_pf]);

	  float dR_KFromB = PFCand_tlv.DeltaR(gen_KFromB_tlv);

	  if(dR_KFromB<0.1 && (best_dR_KFromB<0. || dR_KFromB<best_dR_KFromB)){
	    _PFCand_genKFromB_index = i_pf;
	    best_dR_KFromB = dR_KFromB;
	  }
	}
	*/

	//Gen muon pt filter on tag
	bool isTagMuonHighPt = false;
	for(int i_gen=0; i_gen<nGenPart; i_gen++){

	  int pdgId = tree->GenPart_pdgId[i_gen];
	  if(abs(pdgId)!=13) continue;

	  //exclude the gen muons that are from the chosen B
	  if(!isEleFinalState && 
	     (i_gen == _GenPart_lep1FromB_index || i_gen == _GenPart_lep2FromB_index)) continue;
	  	  
	  TLorentzVector gen_tagMu_tlv;
	  gen_tagMu_tlv.SetPtEtaPhiM(tree->GenPart_pt[i_gen],
				     tree->GenPart_eta[i_gen],
				     tree->GenPart_phi[i_gen],
				     MuonMass_);

	  //In case there are several copies of same muon (with FSR for instance)
	  if(gen_tagMu_tlv.DeltaR(gen_lep1FromB_tlv) > 0.01 && gen_tagMu_tlv.DeltaR(gen_lep2FromB_tlv) > 0.01 &&
	     //	     gen_tagMu_tlv.Pt() > 5.){	     
	     (gen_tagMu_tlv.Pt() > _BToKstll_gen_muonTag_hpT || _BToKstll_gen_muonTag_hpT == -1)){
	    _BToKstll_gen_muonTag_hpT = gen_tagMu_tlv.Pt();
	    isTagMuonHighPt = true;
	    //break;
	  }
	}
	if(!isTagMuonHighPt) continue; //Skip events where the gen filter in MC has been applied to the probe muon
      }


      
      //Require tag muon non matched to gen muons that are from B
      /*    
      if(!isEleFianlState && !isLeptonTrack){
	bool isTagMuonSoftID = false;
	for(int i_mu=0; i_mu<nMuon; ++i_mu){
	  if(i_mu == _Lep_mu1FromB_index || i_mu==_Muon_mu2FromB_index) continue;
	  if(tree->Muon_softId[i_mu] && tree->Muon_pt[i_mu] > 8.){
	    isProbeMuonSoftID = true;
	    _Muon_selgen_index = i_mu;
	    break;
	  }
	  
	}
	if(!isProbeMuonSoftID) continue; //Skip events where there is no probe muon passing the soft ID
      }
      */
    }



    if(isMC == 0){
      _Muon_probe_index = _Muon_sel_index;
    }
    
    tree_new->Fill();
  }


  f_new->cd();
  if(!saveFullNanoAOD) tree_new->AddFriend("Events",input.c_str());

  tree_new->Write();
  f_new->Close();
  return 0;

}



