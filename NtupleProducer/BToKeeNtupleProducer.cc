// Author: T. Strebler (IC)
// Date:   31 May 2018
//
// Add new variables from NanoAOD tree for BToKpipi analysis
// Can be either included in a interpreted macro or compiled in c++
// (use `root-config --glibs --cflags`)
//



#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "NanoAODTree.h"




float JPsiMass_ = 3.0969;
float BuMass_ = 5.279;
float KaonMass_ = 0.493677;
float ElectronMass_ = 0.5109989e-3;

float min_CL = 0.1;
float min_BcosAlpha = 0.9;
float min_KDCA = 6;
float min_BIP = 1;


int main(int argc, char** argv) {

  string status_sample = *(argv + 1);
  bool isMC = false;
  bool isData = false;
  if (status_sample.compare("mc") == 0) isMC = true;
  if (status_sample.compare("data") == 0) isData = true;

  bool isBPHParking = false;
  for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "BPHParking") {
      isBPHParking = true;
      break;
    }
  }


  bool isResonant = false;
  for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "Resonant") {
      isResonant = true;
      break;
    }
  }


  string output;
  for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "--output") {
      if (i + 1 < argc) {
	output = argv[i+1];
	break;
      } else {
	std::cerr << "--output option requires one argument." << std::endl;
	return 1;
      }      
    }  
  }
  if(output==""){
    std::cerr << "--output argument required" << std::endl;
    return 1;
  }
    
  

  string input;
  string listFilesTXT;
  bool inputTXT = false;
  for (int i = 1; i < argc; ++i) {
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

  if(listFilesTXT == "" && input == ""){
    std::cerr << "--input argument required" << std::endl;
    return 1;
  }


  bool overwrite;
  for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "--overwrite") {
      overwrite = true;
      break;
    }
  }


  //Always saveFullNanoAOD info because we are skimming
  bool saveFullNanoAOD = true;
  /*for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "--saveFullNanoAOD") {
      saveFullNanoAOD = true;
      break;
    }
    }*/
 

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
      std::cout << " Adding " << reader << std::endl;
      oldtree->Add(reader.c_str());
    }
  }
  else   oldtree->Add(input.c_str());
  NanoAODTree* tree = new NanoAODTree(oldtree);

  TTree* tree_new=new TTree("BToKeeTree","BToKeeTree");
  if(saveFullNanoAOD)
    tree_new=tree->GetTree()->CloneTree(0);


  //New branches

  int _Muon_sel_index = -1;
  int _BToKee_sel_index = -1;

  tree_new->Branch("Muon_sel_index",&_Muon_sel_index,"Muon_sel_index/I");
  tree_new->Branch("BToKee_sel_index",&_BToKee_sel_index,"BToKee_sel_index/I");

  int _GenPart_BToKee_index = -1;
  int _GenPart_JPsiFromB_index = -1;
  int _GenPart_KFromB_index = -1;
  int _GenPart_e1FromJPsi_index = -1;
  int _GenPart_e2FromJPsi_index = -1;
  int _BToKee_gen_index = -1;

  if(isMC){

    tree_new->Branch("GenPart_BToKee_index",&_GenPart_BToKee_index,"GenPart_BToKee_index/I");
    tree_new->Branch("GenPart_JPsiFromB_index",&_GenPart_JPsiFromB_index,"GenPart_JPsiFromB_index/I");
    tree_new->Branch("GenPart_KFromB_index",&_GenPart_KFromB_index,"GenPart_KFromB_index/I");
    tree_new->Branch("GenPart_e1FromJPsi_index",&_GenPart_e1FromJPsi_index,"GenPart_e1FromJPsi_index/I");
    tree_new->Branch("GenPart_e2FromJPsi_index",&_GenPart_e2FromJPsi_index,"GenPart_e2FromJPsi_index/I");
    tree_new->Branch("BToKee_gen_index",&_BToKee_gen_index,"BToKee_gen_index/I");

  }

  bool _HLT_Mu8p5_IP3p5 = false;
  bool _HLT_Mu10p5_IP3p5 = false;
  bool _HLT_Mu9_IP6 = false;
  bool _HLT_Mu8_IP3 = false;
  bool _HLT_BPHParking = false;

  if(isBPHParking){

    tree_new->Branch("HLT_Mu8p5_IP3p5",&_HLT_Mu8p5_IP3p5,"HLT_Mu8p5_IP3p5/O");
    tree_new->Branch("HLT_Mu10p5_IP3p5",&_HLT_Mu10p5_IP3p5,"HLT_Mu10p5_IP3p5/O");
    tree_new->Branch("HLT_Mu9_IP6",&_HLT_Mu9_IP6,"HLT_Mu9_IP6/O");
    tree_new->Branch("HLT_Mu8_IP3",&_HLT_Mu8_IP3,"HLT_Mu8_IP3/O");
    tree_new->Branch("HLT_BPHParking",&_HLT_BPHParking,"HLT_BPHParking/O");

  }


  int nentries = tree->GetEntries();
  cout<<"Nentries="<<nentries<<endl;
  cout<<"isMC="<<isMC<<endl;

  for (int iEntry = 0; iEntry < nentries ; iEntry++){

    tree->GetEntry(iEntry);

    if(iEntry%10000==0) cout<<"Entry #"<<iEntry<<" "<< int(100*float(iEntry)/nentries)<<"%"<<endl;

    _Muon_sel_index = -1;
    _BToKee_sel_index = -1;

    _GenPart_BToKee_index = -1;
    _GenPart_JPsiFromB_index = -1;
    _GenPart_KFromB_index = -1;
    _GenPart_e1FromJPsi_index = -1;
    _GenPart_e2FromJPsi_index = -1;
    _BToKee_gen_index = -1;

    _HLT_Mu8p5_IP3p5 = false;
    _HLT_Mu10p5_IP3p5 = false;
    _HLT_Mu9_IP6 = false;
    _HLT_Mu8_IP3 = false;
    _HLT_BPHParking = false;

    //Select the muon

    int nMuon = tree->nMuon;

    for(int i_mu=0; i_mu<nMuon; i_mu++){

      bool isMuSel = tree->Muon_softId[i_mu];

      //Trigger selection + matching

      //Only trigger for isBPHParking for now
      if(isBPHParking){

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
	_HLT_BPHParking = _HLT_Mu8p5_IP3p5 || _HLT_Mu10p5_IP3p5 || _HLT_Mu9_IP6 || _HLT_Mu8_IP3;

	TLorentzVector mu;
	mu.SetPtEtaPhiM(tree->Muon_pt[i_mu],tree->Muon_eta[i_mu],tree->Muon_phi[i_mu],tree->Muon_mass[i_mu]);

	bool isTrigMatched = false;
	int nTrigObj = tree->nTrigObj;
	for(unsigned int i_trig = 0; i_trig<nTrigObj; i_trig++){

	  if( tree->TrigObj_id[i_trig]==13 && ((tree->TrigObj_filterBits[i_trig])>>3)&1 ){
	    TLorentzVector trig;
	    trig.SetPtEtaPhiM(tree->TrigObj_pt[i_trig],tree->TrigObj_eta[i_trig],tree->TrigObj_phi[i_trig],0);
	    float dR = mu.DeltaR(trig);
	    if(dR<0.1){
	      isTrigMatched = true;
	      break;
	    }
	  }

	}

	isMuSel &= (_HLT_BPHParking && isTrigMatched);

      }//is Bparking

      //Should implement something to avoid overlap with BToKmm final state (nMuon>1 won't work because only muons with pT>3 GeV are stored in default NanoAOD)

      if( isMuSel ){
	_Muon_sel_index = i_mu;
	break; //Take leading muon passing the selections (muons are pt-ordered)
      }

    }//loop over muons

    //when MC with trigger is available re-enable for everything
    if(_Muon_sel_index <0){
      //Let's skim events which do not have a tag muon (including trigger)
      //tree_new->Fill();
      continue;
    }
  
    //Select the BToKee candidate with reco criteria

    int nBToKee = tree->nBToKee;
    float best_JPsi_mass = -1.;
    float best_Bu_mass = -1.;

    for(int i_BToKee=0; i_BToKee<nBToKee; i_BToKee++){            

      //check for dataset with trigger info
      if(tree->BToKee_kaon_charge[i_BToKee]*tree->Muon_charge[_Muon_sel_index]>0) continue; //Only consider BToKee with opposite charge to muon

      //can be too restrictive for the moment
      /*
      if(BToKee_ele1_charge[i_BToKee] != Electron_charge[BToKee_ele1_charge[i_BToKee]] ||
	 BToKee_ele1_charge[i_BToKee] == -1 || BToKee_ele1_charge[i_BToKee] > nElectron) continue;
      if(BToKee_ele2_charge[i_BToKee] != Electron_charge[BToKee_ele2_charge[i_BToKee]] ||
	 BToKee_ele2_charge[i_BToKee] == -1 || BToKee_ele2_charge[i_BToKee] > nElectron) continue;
      */

      //should study effect of charge flip...
      if(tree->BToKee_ele1_charge[i_BToKee] * tree->BToKee_ele2_charge[i_BToKee] > 0.) continue;


      float ee_mass = tree->BToKee_ee_mass[i_BToKee];
      float ee_CL_vtx = tree->BToKee_ee_CL_vtx[i_BToKee];


      //do not want to bias the ee choice, remove this
      /*
	//JPsi selection
	if( !(best_JPsi_mass < 0. 
	      || abs(best_JPsi_mass-ee_mass)<1e-3 //Several BToKee can share the same JPsi->ee
	      || abs(ee_mass-JPsiMass_) < abs(best_JPsi_mass-JPsiMass_)) )
	  continue;

	//if( ee_CL_vtx < min_CL) continue; //cut on ee vtx refitting
	*/

      float B_mass = tree->BToKee_mass[i_BToKee];
      float B_CL_vtx = tree->BToKee_CL_vtx[i_BToKee];
      
      float BcosAlpha = tree->BToKee_cosAlpha[i_BToKee];
      float B_K_DCAsig = tree->BToKee_kaon_DCASig[i_BToKee];
      float B_VtxIP = tree->BToKee_Lxy[i_BToKee];

      if( !(best_Bu_mass < 0. 
	    || abs(B_mass-BuMass_) < abs(best_Bu_mass-BuMass_)) )       
	continue;

      //can check all of this at next stage on ntuples
      /*
      if( B_CL_vtx < min_CL) continue; //cut on vtx refitting
      if(BcosAlpha < min_BcosAlpha) continue;
      if(B_K_DCAsig < min_KDCA) continue;
      if(B_VtxIP < min_BIP) continue;
      */
      best_JPsi_mass = ee_mass;
      best_Bu_mass = B_mass;
      _BToKee_sel_index = i_BToKee;

    }


    
    //Select the BToKee candidate based on gen matching

    if(isMC){
      
      int nGenPart = tree->nGenPart;

      if(isResonant){
      for(int i_Bu=0; i_Bu<nGenPart; i_Bu++){

	if(abs(tree->GenPart_pdgId[i_Bu])==521){

	  for(int i_gen=0; i_gen<nGenPart; i_gen++){
	    int pdgId = tree->GenPart_pdgId[i_gen];
	    int mother_index = tree->GenPart_genPartIdxMother[i_gen];
	    if(abs(pdgId)==443 && mother_index == i_Bu)
	      _GenPart_JPsiFromB_index = i_gen;
	    else if(abs(pdgId)==321 && mother_index == i_Bu)
	      _GenPart_KFromB_index = i_gen;
	    if(_GenPart_JPsiFromB_index>=0 && _GenPart_KFromB_index>=0){
	      _GenPart_BToKee_index = i_Bu;
	      break;
	    }
	  }	  

	}
	if(_GenPart_BToKee_index>=0) break;
      }

      for(int i_gen=0; i_gen<nGenPart; i_gen++){
	int pdgId = tree->GenPart_pdgId[i_gen];
	int mother_index = tree->GenPart_genPartIdxMother[i_gen];
	if(abs(pdgId)==11 && mother_index == _GenPart_JPsiFromB_index && _GenPart_e1FromJPsi_index<0)
	  _GenPart_e1FromJPsi_index = i_gen;
	else if(abs(pdgId)==11 && mother_index == _GenPart_JPsiFromB_index)
	  _GenPart_e2FromJPsi_index = i_gen;
	if(_GenPart_e1FromJPsi_index>=0 && _GenPart_e2FromJPsi_index>=0) break;

      }
      }//resonant
      else{
	for(int i_Bu=0; i_Bu<nGenPart; i_Bu++){

	  if(abs(tree->GenPart_pdgId[i_Bu])==521){

	    for(int i_gen=0; i_gen<nGenPart; i_gen++){
	      int pdgId = tree->GenPart_pdgId[i_gen];
	      int mother_index = tree->GenPart_genPartIdxMother[i_gen];
	      if(abs(pdgId)==11 && mother_index == i_Bu && _GenPart_e1FromJPsi_index<0)
		_GenPart_e1FromJPsi_index = i_gen;
	      else if(abs(pdgId)==11 && mother_index == i_Bu)
		_GenPart_e2FromJPsi_index = i_gen;
	      else if(abs(pdgId)==321 && mother_index == i_Bu)
		_GenPart_KFromB_index = i_gen;
	      if(_GenPart_e1FromJPsi_index>=0 && _GenPart_e2FromJPsi_index>=0 && _GenPart_KFromB_index>=0){
		_GenPart_BToKee_index = i_Bu;
		break;
	      }
	     }
	  }//if B
	  if(_GenPart_BToKee_index>=0) break;
	}//loop over B
      }


      //e1FromJPsi stored a leading daughter
      if(tree->GenPart_pt[_GenPart_e2FromJPsi_index]>tree->GenPart_pt[_GenPart_e1FromJPsi_index]){
	int i_temp = _GenPart_e1FromJPsi_index;
	_GenPart_e1FromJPsi_index = _GenPart_e2FromJPsi_index;
	_GenPart_e2FromJPsi_index = i_temp;
      }


      TLorentzVector gen_KFromB_tlv;
      TLorentzVector gen_e1FromJPsi_tlv;
      TLorentzVector gen_e2FromJPsi_tlv;
      
      gen_KFromB_tlv.SetPtEtaPhiM(tree->GenPart_pt[_GenPart_KFromB_index],
				   tree->GenPart_eta[_GenPart_KFromB_index],
				   tree->GenPart_phi[_GenPart_KFromB_index],
				   KaonMass_);
      gen_e1FromJPsi_tlv.SetPtEtaPhiM(tree->GenPart_pt[_GenPart_e1FromJPsi_index],
				   tree->GenPart_eta[_GenPart_e1FromJPsi_index],
				   tree->GenPart_phi[_GenPart_e1FromJPsi_index],
				   ElectronMass_);
      gen_e2FromJPsi_tlv.SetPtEtaPhiM(tree->GenPart_pt[_GenPart_e2FromJPsi_index],
				    tree->GenPart_eta[_GenPart_e2FromJPsi_index],
				    tree->GenPart_phi[_GenPart_e2FromJPsi_index],
				    ElectronMass_);

      float best_dR = -1.;

      for(int i_BToKee=0; i_BToKee<nBToKee; i_BToKee++){

	TLorentzVector kaon_tlv;
	TLorentzVector ele1_tlv;
	TLorentzVector ele2_tlv;
	
	kaon_tlv.SetPtEtaPhiM(tree->BToKee_kaon_pt[i_BToKee],
			      tree->BToKee_kaon_eta[i_BToKee],
			      tree->BToKee_kaon_phi[i_BToKee],
			      KaonMass_);
	ele1_tlv.SetPtEtaPhiM(tree->BToKee_ele1_pt[i_BToKee],
			      tree->BToKee_ele1_eta[i_BToKee],
			      tree->BToKee_ele1_phi[i_BToKee],
			      ElectronMass_);
	ele2_tlv.SetPtEtaPhiM(tree->BToKee_ele2_pt[i_BToKee],
			      tree->BToKee_ele2_eta[i_BToKee],
			      tree->BToKee_ele2_phi[i_BToKee],
			      ElectronMass_);

	float dR_KFromB = kaon_tlv.DeltaR(gen_KFromB_tlv);
	float dR_e1FromJPsi = min(ele1_tlv.DeltaR(gen_e1FromJPsi_tlv),ele2_tlv.DeltaR(gen_e1FromJPsi_tlv));
	float dR_e2FromJPsi = min(ele1_tlv.DeltaR(gen_e2FromJPsi_tlv),ele2_tlv.DeltaR(gen_e2FromJPsi_tlv));
	//Should check that same objects not selected twice

	float dR_tot = dR_KFromB + dR_e1FromJPsi + dR_e2FromJPsi; //In case several BToKee matches, take the closest one in dR_tot

	if( dR_KFromB<0.1 && dR_e1FromJPsi<0.1 && dR_e2FromJPsi<0.1
	    && (best_dR<0. || dR_tot<best_dR) ){
	  best_dR = dR_tot;
	  _BToKee_gen_index = i_BToKee;	  
	}

      }

    }

    tree_new->Fill();

  }


  f_new->cd();
  if(!saveFullNanoAOD) tree_new->AddFriend("Events",input.c_str());

  tree_new->Write();
  f_new->Close();
  return 0;

}



