//g++ -Wall -o analyzeCharged_fastDATA `root-config --cflags --glibs` -lRooFitCore analyzeCharged_fastDATA.cpp

//./analyzeCharged_fastDATA --isEle (0,1) --dataset (-1, runA, runB, MC) --run (1,2,3,...) --typeSelection (tightCB, NN_BkgR, NN_SigEff) --ntupleList (list.txt) --JOBid (1,2..) --outputFolder ("outfolder") --nMaxEvents (-1, N) --saveSelectedNTU (1,0) --outSelectedNTU (path for selected ntuples)


#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>

#include "TApplication.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TCut.h"
#include "TError.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TKey.h"
#include "TLegend.h"
#include "TMath.h"
#include "TPad.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

#include "RooAddPdf.h"
#include "RooCBShape.h"
#include "RooChebychev.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooGaussian.h"
#include "RooMsgService.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"

#include <TLorentzVector.h>

using namespace RooFit;

const int kBToKllMax = 50000;
const int kMuonMax = 100;
const int kPFCandMax = 10000;
const int kGenPartMax = 10000;
const float ElectronMass = 0.5109989e-3;
const float MuonMass = 0.10565837;
const float KaonMass = 0.493677;
const float PionMass = 0.139570;

int main(int argc, char **argv){

  if(argc < 2) {
    std::cout << " Missing arguments " << std::endl;
    return -1;
  }
  int isEleFinalState = -1;
  std::string dataset = "-1";
  std::string BPHRun = "-1";
  std::string typeSelection = "-1";
  std::string ntupleList = "-1";
  std::string JOBid = "-1";
  std::string outputFolder = "-1";
  int nMaxEvents = -1;
  int saveOUTntu = 0;
  std::string outSelectedNTU = "-1";
  for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "--isEle") {
      if (i + 1 < argc) {
	isEleFinalState = atoi(argv[i+1]);
	break;
      } else {
	std::cerr << " --isEle option requires one argument " << std::endl;
	return 1;
      }
    }
  }  for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "--dataset") {
      if (i + 1 < argc) {
        dataset = argv[i+1];
        break;
      } else {
	std::cerr << " --dataset option requires one argument " << std::endl;
        return 1;
      }
    }
  }  for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "--run") {
      if (i + 1 < argc) {
        BPHRun = argv[i+1];
        break;
      } else {
	std::cerr << " --run option requires one argument " << std::endl;
        return 1;
      }
    }
  }  for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "--typeSelection") {
      if (i + 1 < argc) {
        typeSelection = argv[i+1];
        break;
      } else {
	std::cerr << " --typeSelection option requires one argument " << std::endl;
        return 1;
      }
    }
  }  for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "--ntupleList") {
      if (i + 1 < argc) {
        ntupleList = argv[i+1];
        break;
      } else {
	std::cerr << " --ntupleList option requires one argument " << std::endl;
        return 1;
      }
    }
  }  for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "--JOBid") {
      if (i + 1 < argc) {
        JOBid = argv[i+1];
        break;
      } else {
	std::cerr << " --JOBid option requires one argument " << std::endl;
        return 1;
      }
    }
  }  for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "--outputFolder") {
      if (i + 1 < argc) {
        outputFolder = argv[i+1];
        break;
      } else {
	std::cerr << " --outputFolder option requires one argument " << std::endl;
        return 1;
      }
    }
  }  for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "--nMaxEvents") {
      if (i + 1 < argc) {
        nMaxEvents = atoi(argv[i+1]);
        break;
      } else {
	std::cerr << " --nMaxEvents option requires one argument " << std::endl;
        return 1;
      }
    }
  }  for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "--saveSelectedNTU") {
      if (i + 1 < argc) {
        saveOUTntu = atoi(argv[i+1]);
        break;
      } else {
	std::cerr << " --saveSelectedNTU option requires one argument " << std::endl;
        return 1;
      }
    }
  }  for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "--outSelectedNTU") {
      if (i + 1 < argc) {
        outSelectedNTU = argv[i+1];
        break;
      } else {
	std::cerr << " --outSelectedNTU option requires one argument " << std::endl;
        return 1;
      }
    }
  }


  if(ntupleList != "-1" && (JOBid == "-1" || outputFolder == "-1")){
    std::cout << " configuration ERROR => splitting file based but missing JOBid and output folder " << std::endl;
    return -1;
  }

  if(saveOUTntu != 0 && outSelectedNTU == "-1"){
    std::cout << " configuration ERROR => missing output folder for final trees " << std::endl;
    return -1;
  }

  std::cout << " isEleFinalState = " << isEleFinalState << " dataset = " << dataset << " BPHRun = " << BPHRun << " typeSelection = " << typeSelection
	    << " ntupleList = " << ntupleList << " JOBid = " << JOBid << " outputFolder = " << outputFolder
	    << " nMaxEvents = " << nMaxEvents << " saveOUTntu = " << saveOUTntu << " outSelectedNTU = " << outSelectedNTU << std::endl;


  gROOT->Reset();
  gROOT->Macro("./setStyle.C");
  gSystem->Load("libRooFit") ;
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TChain* t1 = new TChain("Events");
  //new prod Kee
  if(isEleFinalState){
    if((dataset == "runA" && BPHRun == "1") || dataset == "-1") t1->Add("/vols/cms/tstreble/BPH/BToKee_ntuple/BPHParking1_2018A_18_09_07_elechargefix/*root");
    if((dataset == "runA" && BPHRun == "2") || dataset == "-1") t1->Add("/vols/cms/tstreble/BPH/BToKee_ntuple/BPHParking2_2018A_18_09_07_elechargefix/*root");
    if((dataset == "runA" && BPHRun == "3") || dataset == "-1") t1->Add("/vols/cms/tstreble/BPH/BToKee_ntuple/BPHParking3_2018A_18_09_07_elechargefix/*root");
    if((dataset == "runA" && BPHRun == "4") || dataset == "-1") t1->Add("/vols/cms/tstreble/BPH/BToKee_ntuple/BPHParking4_2018A_18_09_07_elechargefix/*root");
    if((dataset == "runA" && BPHRun == "5") || dataset == "-1") t1->Add("/vols/cms/tstreble/BPH/BToKee_ntuple/BPHParking5_2018A_18_09_07_elechargefix/*root");
    if((dataset == "runA" && BPHRun == "6") || dataset == "-1") t1->Add("/vols/cms/tstreble/BPH/BToKee_ntuple/BPHParking6_2018A_18_09_07_elechargefix/*root");

    if((dataset == "runB" && BPHRun == "1") || dataset == "-1") t1->Add("/vols/cms/amartell/BParking/ntuPROD/BPHrun2018B/BToKee_2018B_BPH1/*root");
    if((dataset == "runB" && BPHRun == "2") || dataset == "-1") t1->Add("/vols/cms/amartell/BParking/ntuPROD/BPHrun2018B/BToKee_2018B_BPH2/*root");
    if((dataset == "runB" && BPHRun == "3") || dataset == "-1") t1->Add("/vols/cms/amartell/BParking/ntuPROD/BPHrun2018B/BToKee_2018B_BPH3/*root");
    if((dataset == "runB" && BPHRun == "4") || dataset == "-1") t1->Add("/vols/cms/amartell/BParking/ntuPROD/BPHrun2018B/BToKee_2018B_BPH4/*root");
    if((dataset == "runB" && BPHRun == "5") || dataset == "-1") t1->Add("/vols/cms/amartell/BParking/ntuPROD/BPHrun2018B/BToKee_2018B_BPH5/*root");

    if(dataset == "MC" && BPHRun == "nnResonant"){
      t1->Add("/vols/cms/amartell/BParking/ntuPROD/newNANO_20Aug/ntu_BToKee_18_09_10.root");
    }
    if(dataset == "MC" && BPHRun == "Resonant"){
      t1->Add("/vols/cms/amartell/BParking/ntuPROD/newNANO_20Aug/ntu_BToKJPsiee_18_09_10.root");
    }
  }
  else{
    //round1 ntuples
    // if((dataset == "runA" && BPHRun == "1") || dataset == "-1") t1->Add("/vols/cms/tstreble/BPH/BToKmumu_ntuple/BPHParking1_2018A_18_08_14_new/*root");
    // if((dataset == "runA" && BPHRun == "2") || dataset == "-1") t1->Add("/vols/cms/tstreble/BPH/BToKmumu_ntuple/BPHParking2_2018A_18_08_14_new/*root");
    // if((dataset == "runA" && BPHRun == "3") || dataset == "-1") t1->Add("/vols/cms/tstreble/BPH/BToKmumu_ntuple/BPHParking3_2018A_18_08_14_new/*root");
    // if((dataset == "runA" && BPHRun == "4") || dataset == "-1") t1->Add("/vols/cms/tstreble/BPH/BToKmumu_ntuple/BPHParking4_2018A_18_08_14_new/*root");
    // if((dataset == "runA" && BPHRun == "5") || dataset == "-1") t1->Add("/vols/cms/tstreble/BPH/BToKmumu_ntuple/BPHParking5_2018A_18_08_14_new/*root");
    // if((dataset == "runA" && BPHRun == "6") || dataset == "-1") t1->Add("/vols/cms/tstreble/BPH/BToKmumu_ntuple/BPHParking6_2018A_18_08_14_new/*root");

    // if((dataset == "runB" && BPHRun == "1") || dataset == "-1") t1->Add("/vols/cms/tstreble/BPH/BToKmumu_ntuple/BPHParking1_2018B_18_08_14_new/*root");
    // if((dataset == "runB" && BPHRun == "2") || dataset == "-1") t1->Add("/vols/cms/tstreble/BPH/BToKmumu_ntuple/BPHParking2_2018B_18_08_14_new/*root");
    // if((dataset == "runB" && BPHRun == "3") || dataset == "-1") t1->Add("/vols/cms/tstreble/BPH/BToKmumu_ntuple/BPHParking3_2018B_18_08_14_new/*root");
    // if((dataset == "runB" && BPHRun == "4") || dataset == "-1") t1->Add("/vols/cms/tstreble/BPH/BToKmumu_ntuple/BPHParking4_2018B_18_08_14_new/*root");
    // if((dataset == "runB" && BPHRun == "5") || dataset == "-1") t1->Add("/vols/cms/tstreble/BPH/BToKmumu_ntuple/BPHParking5_2018B_18_08_14_new/*root");

    //need update wrt new production => check https://docs.google.com/spreadsheets/d/1Kdtaw0nGNXZ_O5-e7DR5WGOI0AOpUhccMfn4WIieGYU/edit#gid=0
    if(ntupleList == "-1"){
      if((dataset == "runA" && BPHRun == "1") || dataset == "-1") t1->Add("/vols/cms/vc1116/BParking/ntuPROD/data_BToKmumuNtuple/old_BToKmumu/A1/*root");
      if((dataset == "runA" && BPHRun == "2") || dataset == "-1") t1->Add("/vols/cms/vc1116/BParking/ntuPROD/data_BToKmumuNtuple/old_BToKmumu/A2/*root");
      if((dataset == "runA" && BPHRun == "3") || dataset == "-1") t1->Add("/vols/cms/vc1116/BParking/ntuPROD/data_BToKmumuNtuple/old_BToKmumu/A3/*root");
      if((dataset == "runA" && BPHRun == "4") || dataset == "-1") t1->Add("/vols/cms/vc1116/BParking/ntuPROD/data_BToKmumuNtuple/old_BToKmumu/A4/*root");
      if((dataset == "runA" && BPHRun == "5") || dataset == "-1") t1->Add("/vols/cms/vc1116/BParking/ntuPROD/data_BToKmumuNtuple/old_BToKmumu/A5/*root");
      if((dataset == "runA" && BPHRun == "6") || dataset == "-1") t1->Add("/vols/cms/vc1116/BParking/ntuPROD/data_BToKmumuNtuple/old_BToKmumu/A6/*root");

      if((dataset == "runB" && BPHRun == "1") || dataset == "-1") t1->Add("/vols/cms/amartell/BParking/ntuPROD/BPHrun2018B/BToKmumu_2018B_BPH1_NN/*root");
      if((dataset == "runB" && BPHRun == "2") || dataset == "-1") t1->Add("/vols/cms/amartell/BParking/ntuPROD/BPHrun2018B/BToKmumu_2018B_BPH2_NN/*root");
      if((dataset == "runB" && BPHRun == "3") || dataset == "-1") t1->Add("/vols/cms/amartell/BParking/ntuPROD/BPHrun2018B/BToKmumu_2018B_BPH3_NN/*root");
      if((dataset == "runB" && BPHRun == "4") || dataset == "-1") t1->Add("/vols/cms/amartell/BParking/ntuPROD/BPHrun2018B/BToKmumu_2018B_BPH4_NN/*root");
      if((dataset == "runB" && BPHRun == "5") || dataset == "-1") t1->Add("/vols/cms/amartell/BParking/ntuPROD/BPHrun2018B/BToKmumu_2018B_BPH5_NN/*root");
    }

    if(ntupleList != "-1"){
      std::string rootFileName;
      std::ifstream inFileLong;
      inFileLong.open(ntupleList.c_str(), std::ios::in);
      while(!inFileLong.eof()){
	if(inFileLong >> rootFileName){
	  t1->Add(rootFileName.c_str());
	  std::cout << " adding " << rootFileName << std::endl;
	}
      }
    }

    if(dataset == "MC" && BPHRun == "nnResonant"){
      t1->Add("/vols/cms/tstreble/BPH/BToKmumu_ntuple/BToKmumu_18_08_14_new/*root");
    }
    if(dataset == "MC" && BPHRun == "Resonant"){
      t1->Add("/vols/cms/tstreble/BPH/BToKmumu_ntuple/BToKJPsimumu_18_08_14_new/*root");
    }
  }

  int nEvts = t1->GetEntries();
  std::cout << " #initial n. events: " << nEvts << std::endl;

  std::string outNtuName = "selectedEvents_Kee_"+dataset+"_BPHRun"+BPHRun;
  if(!isEleFinalState) outNtuName = "selectedEvents_Kmumu_"+dataset+"_BPHRun"+BPHRun;
  if(JOBid != "-1") outNtuName = outSelectedNTU+"/"+outNtuName+"_JOB_"+JOBid;
  outNtuName += ".root";
  gROOT->cd();
  TFile* newFile;
  TTree* newT;

  Float_t Kll_cosAlpha;
  Float_t Kll_CL_vtx;
  Float_t Kll_kaon_DCASig;
  Float_t Kll_Lxy;
  Float_t Kll_ctxy;
  Float_t Kll_llmass;
  Float_t Kll_llRefitmass;
  Float_t Kll_llPrefitmass;
  Float_t Kll_mass;
  Float_t Kll_pt;
  Float_t Kll_kaon_pt;
  Float_t Kll_kaon_eta;
  Float_t Kll_kaon_phi;
  Int_t Kll_kaon_charge;
  Float_t Kll_lep1_pt;
  Float_t Kll_lep2_pt;
  Float_t Kll_lep1_eta;
  Float_t Kll_lep2_eta;
  Float_t Kll_lep1_phi;
  Float_t Kll_lep2_phi;
  Float_t Kll_lep1_preFpt;
  Float_t Kll_lep2_preFpt;
  Float_t Kll_lep1_preFeta;
  Float_t Kll_lep2_preFeta;
  Float_t Kll_lep1_preFphi;
  Float_t Kll_lep2_preFphi;
  Float_t Kll_lep1_pfRelIso03;
  Float_t Kll_lep2_pfRelIso03;
  Float_t Kll_lep1_pfRelIso04;
  Float_t Kll_lep2_pfRelIso04;
  Float_t Kll_lep1_dxy;
  Float_t Kll_lep2_dxy;
  Float_t Kll_lep1_dz;
  Float_t Kll_lep2_dz;
  Float_t Kll_kaon_dxy;
  Float_t Kll_kaon_dz;

  if(saveOUTntu){
    newFile = new TFile(outNtuName.c_str(),"recreate");
    newT = new TTree("selectedEvents", "");
    newT->Branch("Kll_cosAlpha", &Kll_cosAlpha, "Kll_cosAlpha/F");
    newT->Branch("Kll_CL_vtx", &Kll_CL_vtx, "Kll_CL_vtx/F");
    newT->Branch("Kll_kaon_DCASig", &Kll_kaon_DCASig, "Kll_kaon_DCASig/F");       
    newT->Branch("Kll_Lxy", &Kll_Lxy, "Kll_Lxy/F");
    newT->Branch("Kll_ctxy", &Kll_ctxy, "Kll_ctxy/F");
    newT->Branch("Kll_llRefitmass", &Kll_llRefitmass, "Kll_llRefitmass/F");
    newT->Branch("Kll_llPrefitmass", &Kll_llPrefitmass, "Kll_llPrefitmass/F");
    newT->Branch("Kll_llmass", &Kll_llmass, "Kll_llmass/F");
    newT->Branch("Kll_mass", &Kll_mass, "Kll_mass/F");
    newT->Branch("Kll_pt", &Kll_pt, "Kll_pt/F");                
    newT->Branch("Kll_kaon_pt", &Kll_kaon_pt, "Kll_kaon_pt/F");
    newT->Branch("Kll_kaon_eta", &Kll_kaon_eta, "Kll_kaon_eta/F");
    newT->Branch("Kll_kaon_phi", &Kll_kaon_phi, "Kll_kaon_phi/F");
    newT->Branch("Kll_kaon_charge", &Kll_kaon_charge, "Kll_kaon_charge/I");
    newT->Branch("Kll_lep1_pt", &Kll_lep1_pt, "Kll_lep1_pt/F");
    newT->Branch("Kll_lep2_pt", &Kll_lep2_pt, "Kll_lep2_pt/F");
    newT->Branch("Kll_lep1_eta", &Kll_lep1_eta, "Kll_lep1_eta/F");
    newT->Branch("Kll_lep2_eta", &Kll_lep2_eta, "Kll_lep2_eta/F");
    newT->Branch("Kll_lep1_phi", &Kll_lep1_phi, "Kll_lep1_phi/F");
    newT->Branch("Kll_lep2_phi", &Kll_lep2_phi, "Kll_lep2_phi/F");
    newT->Branch("Kll_lep1_preFpt", &Kll_lep1_preFpt, "Kll_lep1_preFpt/F");
    newT->Branch("Kll_lep2_preFpt", &Kll_lep2_preFpt, "Kll_lep2_preFpt/F");
    newT->Branch("Kll_lep1_preFeta", &Kll_lep1_preFeta, "Kll_lep1_preFeta/F");
    newT->Branch("Kll_lep2_preFeta", &Kll_lep2_preFeta, "Kll_lep2_preFeta/F");
    newT->Branch("Kll_lep1_preFphi", &Kll_lep1_preFphi, "Kll_lep1_preFphi/F");
    newT->Branch("Kll_lep2_preFphi", &Kll_lep2_preFphi, "Kll_lep2_preFphi/F");
    newT->Branch("Kll_lep1_pfRelIso03", &Kll_lep1_pfRelIso03, "Kll_lep1_pfRelIso03/F");
    newT->Branch("Kll_lep2_pfRelIso03", &Kll_lep2_pfRelIso03, "Kll_lep2_pfRelIso03/F");
    newT->Branch("Kll_lep1_pfRelIso04", &Kll_lep1_pfRelIso04, "Kll_lep1_pfRelIso04/F");
    newT->Branch("Kll_lep2_pfRelIso04", &Kll_lep2_pfRelIso04, "Kll_lep2_pfRelIso04/F");
    newT->Branch("Kll_lep1_dxy", &Kll_lep1_dxy, "Kll_lep1_dxy/F");
    newT->Branch("Kll_lep2_dxy", &Kll_lep2_dxy, "Kll_lep2_dxy/F");
    newT->Branch("Kll_lep1_dz", &Kll_lep1_dz, "Kll_lep1_dz/F");
    newT->Branch("Kll_lep2_dz", &Kll_lep2_dz, "Kll_lep2_dz/F");
    newT->Branch("Kll_kaon_dxy", &Kll_kaon_dxy, "Kll_kaon_dxy/F");
    newT->Branch("Kll_kaon_dz", &Kll_kaon_dz, "Kll_kaon_dz/F");
  }

  std::cout << " >>> qui ok " << std::endl;

  UInt_t run = 0;
  UInt_t lumi = 0;
  ULong64_t event = 0;

  float nnBMX = -1;
  int MuonTag_index = -1;
  int MuonRecoTag_index = -1;
  int BToKll_sel_index = -1;
  int BToKll_gen_index = -1;
  bool Lepton_softId[kMuonMax];
  bool Lepton_mediumId[kMuonMax];
  float BToKll_lep1_charge[kBToKllMax];
  float BToKll_lep2_charge[kBToKllMax];
  int BToKll_lep1_index[kBToKllMax];
  int BToKll_lep2_index[kBToKllMax];
  int BToKll_kaon_index[kBToKllMax];
  float BToKll_cosAlpha[kBToKllMax];
  float BToKll_CL_vtx[kBToKllMax];
  float BToKll_kaon_DCASig[kBToKllMax];
  float BToKll_Lxy[kBToKllMax];
  float BToKll_ctxy[kBToKllMax];
  float BToKll_llKFit_ll_mass[kBToKllMax];
  float BToKll_ll_mass[kBToKllMax];
  int passed_2trk[kBToKllMax];
  float BToKll_mass[kBToKllMax];
  float BToKll_pt[kBToKllMax];
  float BToKll_kaon_pt[kBToKllMax];
  float BToKll_kaon_eta[kBToKllMax];
  float BToKll_kaon_phi[kBToKllMax];
  int	BToKll_kaon_charge[kBToKllMax];
  float BToKll_lep1_pt[kBToKllMax];
  float BToKll_lep2_pt[kBToKllMax];
  float BToKll_lep1_eta[kBToKllMax];
  float BToKll_lep2_eta[kBToKllMax];
  float BToKll_lep1_phi[kBToKllMax];
  float BToKll_lep2_phi[kBToKllMax];
  float Lepton_pfRelIso03[kMuonMax];
  float Lepton_pfRelIso04[kMuonMax];
  float Lepton_pt[kMuonMax];
  float Lepton_eta[kMuonMax];
  float Lepton_phi[kMuonMax];
  float Lepton_mass[kMuonMax];
  int Lepton_charge[kMuonMax];
  float Lepton_dxy[kMuonMax];
  float Lepton_dz[kMuonMax];
  float PFCand_pt[kPFCandMax];
  float PFCand_eta[kPFCandMax];
  float PFCand_phi[kPFCandMax];
  float PFCand_mass[kPFCandMax];
  float PFCand_dxy[kPFCandMax];
  float PFCand_dz[kPFCandMax];
  float GenPart_pt[kGenPartMax];
  float GenPart_eta[kGenPartMax];
  float GenPart_phi[kGenPartMax];
  int GenPart_KFromB_index = -1;
  float BToKll_gen_mass[kBToKllMax];
  float BToKll_gen_llMass[kBToKllMax];
  
  t1->SetBranchStatus("*", 0);


  t1->SetBranchStatus("run", 1);                        t1->SetBranchAddress("run", &run);
  t1->SetBranchStatus("luminosityBlock", 1);            t1->SetBranchAddress("luminosityBlock", &lumi);
  t1->SetBranchStatus("event", 1);                      t1->SetBranchAddress("event", &event);

  if(isEleFinalState){
    t1->SetBranchStatus("Muon_sel_index", 1);            t1->SetBranchAddress("Muon_sel_index", &MuonTag_index);
    t1->SetBranchStatus("BToKee_sel_index", 1);          t1->SetBranchAddress("BToKee_sel_index", &BToKll_sel_index);
    if(dataset == "MC"){
      t1->SetBranchStatus("BToKee_gen_index", 1);          t1->SetBranchAddress("BToKee_gen_index", &BToKll_gen_index);
      t1->SetBranchStatus("GenPart_KFromB_index", 1);      t1->SetBranchAddress("GenPart_KFromB_index", &GenPart_KFromB_index);
      t1->SetBranchStatus("GenPart_pt", 1);                t1->SetBranchAddress("GenPart_pt", &GenPart_pt);
      t1->SetBranchStatus("GenPart_eta", 1);               t1->SetBranchAddress("GenPart_eta", &GenPart_eta);
      t1->SetBranchStatus("GenPart_phi", 1);               t1->SetBranchAddress("GenPart_phi", &GenPart_phi);
      t1->SetBranchStatus("BToKee_gen_mass", 1);           t1->SetBranchAddress("BToKee_gen_mass", &BToKll_gen_mass);
      t1->SetBranchStatus("BToKee_gen_eeMass", 1);         t1->SetBranchAddress("BToKee_gen_eeMass", &BToKll_gen_llMass);
    }
    t1->SetBranchStatus("BToKee_ele1_charge", 1);        t1->SetBranchAddress("BToKee_ele1_charge", &BToKll_lep1_charge);
    t1->SetBranchStatus("BToKee_ele2_charge", 1);        t1->SetBranchAddress("BToKee_ele2_charge", &BToKll_lep2_charge);
    t1->SetBranchStatus("BToKee_ele1_index", 1);         t1->SetBranchAddress("BToKee_ele1_index", &BToKll_lep1_index);
    t1->SetBranchStatus("BToKee_ele2_index", 1);         t1->SetBranchAddress("BToKee_ele2_index", &BToKll_lep2_index);
    t1->SetBranchStatus("BToKee_kaon_index", 1);         t1->SetBranchAddress("BToKee_kaon_index", &BToKll_kaon_index);
    t1->SetBranchStatus("BToKee_cosAlpha", 1);           t1->SetBranchAddress("BToKee_cosAlpha", &BToKll_cosAlpha);
    t1->SetBranchStatus("BToKee_CL_vtx", 1);             t1->SetBranchAddress("BToKee_CL_vtx", &BToKll_CL_vtx);
    t1->SetBranchStatus("BToKee_kaon_DCASig", 1);        t1->SetBranchAddress("BToKee_kaon_DCASig", &BToKll_kaon_DCASig);
    t1->SetBranchStatus("BToKee_Lxy", 1);                t1->SetBranchAddress("BToKee_Lxy", &BToKll_Lxy);
    t1->SetBranchStatus("BToKee_ctxy", 1);               t1->SetBranchAddress("BToKee_ctxy", &BToKll_ctxy);
    t1->SetBranchStatus("BToKee_eeKFit_ee_mass", 1);     t1->SetBranchAddress("BToKee_eeKFit_ee_mass", &BToKll_llKFit_ll_mass);
    t1->SetBranchStatus("BToKee_ee_mass", 1);            t1->SetBranchAddress("BToKee_ee_mass", &BToKll_ll_mass);
    t1->SetBranchStatus("BToKee_eeRefit", 1);            t1->SetBranchAddress("BToKee_eeRefit", &passed_2trk);
    t1->SetBranchStatus("BToKee_mass", 1);               t1->SetBranchAddress("BToKee_mass", &BToKll_mass);
    t1->SetBranchStatus("BToKee_pt", 1);                 t1->SetBranchAddress("BToKee_pt", &BToKll_pt);
    t1->SetBranchStatus("BToKee_kaon_pt", 1);            t1->SetBranchAddress("BToKee_kaon_pt", &BToKll_kaon_pt);
    t1->SetBranchStatus("BToKee_kaon_eta", 1);           t1->SetBranchAddress("BToKee_kaon_eta", &BToKll_kaon_eta);
    t1->SetBranchStatus("BToKee_kaon_phi", 1);           t1->SetBranchAddress("BToKee_kaon_phi", &BToKll_kaon_phi);
    t1->SetBranchStatus("BToKee_kaon_charge", 1);        t1->SetBranchAddress("BToKee_kaon_charge", &BToKll_kaon_charge);
    t1->SetBranchStatus("BToKee_ele1_pt", 1);            t1->SetBranchAddress("BToKee_ele1_pt", &BToKll_lep1_pt);
    t1->SetBranchStatus("BToKee_ele2_pt", 1);            t1->SetBranchAddress("BToKee_ele2_pt", &BToKll_lep2_pt);
    t1->SetBranchStatus("BToKee_ele1_eta", 1);            t1->SetBranchAddress("BToKee_ele1_eta", &BToKll_lep1_eta);
    t1->SetBranchStatus("BToKee_ele2_eta", 1);            t1->SetBranchAddress("BToKee_ele2_eta", &BToKll_lep2_eta);
    t1->SetBranchStatus("BToKee_ele1_phi", 1);            t1->SetBranchAddress("BToKee_ele1_phi", &BToKll_lep1_phi);
    t1->SetBranchStatus("BToKee_ele2_phi", 1);            t1->SetBranchAddress("BToKee_ele2_phi", &BToKll_lep2_phi);
    t1->SetBranchStatus("Electron_pfRelIso03_all", 1);    t1->SetBranchAddress("Electron_pfRelIso03_all", &Lepton_pfRelIso03);
    t1->SetBranchStatus("Electron_pt", 1);                t1->SetBranchAddress("Electron_pt", &Lepton_pt);
    t1->SetBranchStatus("Electron_eta", 1);               t1->SetBranchAddress("Electron_eta", &Lepton_eta);
    t1->SetBranchStatus("Electron_phi", 1);               t1->SetBranchAddress("Electron_phi", &Lepton_phi);
    t1->SetBranchStatus("Electron_mass", 1);              t1->SetBranchAddress("Electron_mass", &Lepton_mass);
    t1->SetBranchStatus("Electron_charge", 1);            t1->SetBranchAddress("Electron_charge", &Lepton_charge);
    t1->SetBranchStatus("Electron_dxy", 1);               t1->SetBranchAddress("Electron_dxy", &Lepton_dxy);
    t1->SetBranchStatus("Electron_dz", 1);                t1->SetBranchAddress("Electron_dz", &Lepton_dz);
    t1->SetBranchStatus("PFCand_pt", 1);                  t1->SetBranchAddress("PFCand_pt", &PFCand_pt);
    t1->SetBranchStatus("PFCand_eta", 1);                 t1->SetBranchAddress("PFCand_eta", &PFCand_eta);
    t1->SetBranchStatus("PFCand_phi", 1);                 t1->SetBranchAddress("PFCand_phi", &PFCand_phi);
    t1->SetBranchStatus("PFCand_mass", 1);                t1->SetBranchAddress("PFCand_mass", &PFCand_mass);
    t1->SetBranchStatus("PFCand_dxy", 1);                 t1->SetBranchAddress("PFCand_dxy", &PFCand_dxy);
    t1->SetBranchStatus("PFCand_dz", 1);                  t1->SetBranchAddress("PFCand_dz", &PFCand_dz);
  }
  else{
    t1->SetBranchStatus("Muon_probe_index", 1);            t1->SetBranchAddress("Muon_probe_index", &MuonTag_index);
    t1->SetBranchStatus("Muon_sel_index", 1);              t1->SetBranchAddress("Muon_sel_index", &MuonRecoTag_index);
    t1->SetBranchStatus("Muon_softId", 1);                 t1->SetBranchAddress("Muon_softId", &Lepton_softId);
    t1->SetBranchStatus("Muon_mediumId", 1);               t1->SetBranchAddress("Muon_mediumId", &Lepton_mediumId);
    t1->SetBranchStatus("BToKmumu_sel_index", 1);          t1->SetBranchAddress("BToKmumu_sel_index", &BToKll_sel_index);
    if(typeSelection == "NN_BkgR" || typeSelection == "NN_SigEff"){
      t1->SetBranchStatus("nnBMX", 1);                      t1->SetBranchAddress("nnBMX", &nnBMX);
    }
    if(dataset == "MC"){
      t1->SetBranchStatus("BToKmumu_gen_index", 1);          t1->SetBranchAddress("BToKmumu_gen_index", &BToKll_gen_index);
      t1->SetBranchStatus("GenPart_KFromB_index", 1);      t1->SetBranchAddress("GenPart_KFromB_index", &GenPart_KFromB_index);
      t1->SetBranchStatus("GenPart_pt", 1);                t1->SetBranchAddress("GenPart_pt", &GenPart_pt);
      t1->SetBranchStatus("GenPart_eta", 1);               t1->SetBranchAddress("GenPart_eta", &GenPart_eta);
      t1->SetBranchStatus("GenPart_phi", 1);               t1->SetBranchAddress("GenPart_phi", &GenPart_phi);
      t1->SetBranchStatus("BToKmumu_gen_mass", 1);         t1->SetBranchAddress("BToKmumu_gen_mass", BToKll_gen_mass);
      t1->SetBranchStatus("BToKmumu_gen_mumuMass", 1);     t1->SetBranchAddress("BToKmumu_gen_mumuMass", BToKll_gen_llMass);
    }
    t1->SetBranchStatus("BToKmumu_mu1_charge", 1);         t1->SetBranchAddress("BToKmumu_mu1_charge", &BToKll_lep1_charge);
    t1->SetBranchStatus("BToKmumu_mu2_charge", 1);         t1->SetBranchAddress("BToKmumu_mu2_charge", &BToKll_lep2_charge);
    t1->SetBranchStatus("BToKmumu_mu1_index", 1);         t1->SetBranchAddress("BToKmumu_mu1_index", &BToKll_lep1_index);
    t1->SetBranchStatus("BToKmumu_mu2_index", 1);         t1->SetBranchAddress("BToKmumu_mu2_index", &BToKll_lep2_index);
    t1->SetBranchStatus("BToKmumu_kaon_index", 1);         t1->SetBranchAddress("BToKmumu_kaon_index", &BToKll_kaon_index);
    t1->SetBranchStatus("BToKmumu_cosAlpha", 1);           t1->SetBranchAddress("BToKmumu_cosAlpha", &BToKll_cosAlpha);
    t1->SetBranchStatus("BToKmumu_CL_vtx", 1);             t1->SetBranchAddress("BToKmumu_CL_vtx", &BToKll_CL_vtx);
    t1->SetBranchStatus("BToKmumu_kaon_DCASig", 1);        t1->SetBranchAddress("BToKmumu_kaon_DCASig", &BToKll_kaon_DCASig);
    t1->SetBranchStatus("BToKmumu_Lxy", 1);                t1->SetBranchAddress("BToKmumu_Lxy", &BToKll_Lxy);
    t1->SetBranchStatus("BToKmumu_ctxy", 1);               t1->SetBranchAddress("BToKmumu_ctxy", &BToKll_ctxy);
    t1->SetBranchStatus("BToKmumu_mumuKFit_mumu_mass", 1); t1->SetBranchAddress("BToKmumu_mumuKFit_mumu_mass", &BToKll_llKFit_ll_mass);
    t1->SetBranchStatus("BToKmumu_mumu_mass", 1);          t1->SetBranchAddress("BToKmumu_mumu_mass", &BToKll_ll_mass);
    t1->SetBranchStatus("BToKmumu_mumuRefit", 1);          t1->SetBranchAddress("BToKmumu_mumuRefit", &passed_2trk);
    t1->SetBranchStatus("BToKmumu_mass", 1);               t1->SetBranchAddress("BToKmumu_mass", &BToKll_mass);
    t1->SetBranchStatus("BToKmumu_pt", 1);                 t1->SetBranchAddress("BToKmumu_pt", &BToKll_pt);
    t1->SetBranchStatus("BToKmumu_kaon_pt", 1);            t1->SetBranchAddress("BToKmumu_kaon_pt", &BToKll_kaon_pt);
    t1->SetBranchStatus("BToKmumu_kaon_eta", 1);           t1->SetBranchAddress("BToKmumu_kaon_eta", &BToKll_kaon_eta);
    t1->SetBranchStatus("BToKmumu_kaon_phi", 1);           t1->SetBranchAddress("BToKmumu_kaon_phi", &BToKll_kaon_phi);
    t1->SetBranchStatus("BToKmumu_kaon_charge", 1);        t1->SetBranchAddress("BToKmumu_kaon_charge", &BToKll_kaon_charge);
    t1->SetBranchStatus("BToKmumu_mu1_pt", 1);           t1->SetBranchAddress("BToKmumu_mu1_pt", &BToKll_lep1_pt);
    t1->SetBranchStatus("BToKmumu_mu2_pt", 1);           t1->SetBranchAddress("BToKmumu_mu2_pt", &BToKll_lep2_pt);
    t1->SetBranchStatus("BToKmumu_mu1_eta", 1);            t1->SetBranchAddress("BToKmumu_mu1_eta", &BToKll_lep1_eta);
    t1->SetBranchStatus("BToKmumu_mu2_eta", 1);            t1->SetBranchAddress("BToKmumu_mu2_eta", &BToKll_lep2_eta);
    t1->SetBranchStatus("BToKmumu_mu1_phi", 1);            t1->SetBranchAddress("BToKmumu_mu1_phi", &BToKll_lep1_phi);
    t1->SetBranchStatus("BToKmumu_mu2_phi", 1);            t1->SetBranchAddress("BToKmumu_mu2_phi", &BToKll_lep2_phi);
    t1->SetBranchStatus("Muon_pfRelIso03_all", 1);         t1->SetBranchAddress("Muon_pfRelIso03_all", &Lepton_pfRelIso03);
    t1->SetBranchStatus("Muon_pfRelIso04_all", 1);         t1->SetBranchAddress("Muon_pfRelIso04_all", &Lepton_pfRelIso04);
    t1->SetBranchStatus("Muon_pt", 1);                     t1->SetBranchAddress("Muon_pt", &Lepton_pt);
    t1->SetBranchStatus("Muon_eta", 1);                    t1->SetBranchAddress("Muon_eta", &Lepton_eta);
    t1->SetBranchStatus("Muon_phi", 1);                    t1->SetBranchAddress("Muon_phi", &Lepton_phi);
    t1->SetBranchStatus("Muon_mass", 1);                    t1->SetBranchAddress("Muon_mass", &Lepton_mass);
    t1->SetBranchStatus("Muon_charge", 1);                 t1->SetBranchAddress("Muon_charge", &Lepton_charge);
    t1->SetBranchStatus("Muon_dxy", 1);                    t1->SetBranchAddress("Muon_dxy", &Lepton_dxy);
    t1->SetBranchStatus("Muon_dz", 1);                     t1->SetBranchAddress("Muon_dz", &Lepton_dz);
    t1->SetBranchStatus("PFCand_pt", 1);                t1->SetBranchAddress("PFCand_pt", &PFCand_pt);
    t1->SetBranchStatus("PFCand_eta", 1);                t1->SetBranchAddress("PFCand_eta", &PFCand_eta);
    t1->SetBranchStatus("PFCand_phi", 1);                t1->SetBranchAddress("PFCand_phi", &PFCand_phi);
    t1->SetBranchStatus("PFCand_mass", 1);                t1->SetBranchAddress("PFCand_mass", &PFCand_mass);
    t1->SetBranchStatus("PFCand_dxy", 1);                t1->SetBranchAddress("PFCand_dxy", &PFCand_dxy);
    t1->SetBranchStatus("PFCand_dz", 1);                t1->SetBranchAddress("PFCand_dz", &PFCand_dz);
  }



  std::vector<float> llMassBoundary;
  llMassBoundary.push_back(0.);
  llMassBoundary.push_back(1.);
  llMassBoundary.push_back(2.5);
  llMassBoundary.push_back(2.9);
  llMassBoundary.push_back(3.3);
  llMassBoundary.push_back(3.58);
  llMassBoundary.push_back(100.);


  std::string outName = "outMassHistos_Kee_"+dataset+"_BPHRun"+BPHRun;
  if(!isEleFinalState) outName = "outMassHistos_Kmumu_"+dataset+"_BPHRun"+BPHRun;
  if(JOBid != "-1") outName = outputFolder; // + "/" +outName + "_JOB_"+JOBid;
  else  outName += ".root";
  TFile outMassHistos(outName.c_str(), "recreate");

  ///histos: 1 per bin plus inclusive
  TH1F* hAlpha[7];
  TH1F* hCLVtx[7];
  TH1F* hDCASig[7];
  TH1F* hLxy[7];
  TH1F* hctxy[7];
  TH1F* hKaonpt[7];
  TH1F* hLep1pt[7];
  TH1F* hLep2pt[7];
  TH1F* hLep1pt_EB[7];
  TH1F* hLep2pt_EB[7];
  TH1F* hLep1pt_EE[7];
  TH1F* hLep2pt_EE[7];
  TH1F* hBpt[7];
  TH1F* hllMass[7];
  TH1F* hllRefitMass[7];
  TH1F* hllPrefitMass[7];
  TH2F* hllMass_vs_Bmass[7];
  TH1F* hBmass[7];

  for(int ij=0; ij<7; ++ij){
    hAlpha[ij] = new TH1F(Form("hAlpha_bin%d", ij), "", 500, 0, 1.1);
    hAlpha[ij]->Sumw2();
    hAlpha[ij]->SetLineColor(kRed);
    hAlpha[ij]->SetLineWidth(2);

    hCLVtx[ij] = new TH1F(Form("hCLVtx_%d", ij), "", 100, 0., 1.);
    hCLVtx[ij]->Sumw2();
    hCLVtx[ij]->SetLineColor(kRed);
    hCLVtx[ij]->SetLineWidth(2);
 
    hDCASig[ij] = new TH1F(Form("hDCASig_%d", ij), "", 100, -50., 50.);
    hDCASig[ij]->Sumw2();
    hDCASig[ij]->SetLineColor(kRed);
    hDCASig[ij]->SetLineWidth(2);

    hLxy[ij] = new TH1F(Form("hLxy_%d", ij), "", 100, 0., 100.);
    hLxy[ij]->Sumw2();
    hLxy[ij]->SetLineColor(kRed);
    hLxy[ij]->SetLineWidth(2);

    hllRefitMass[ij] = new TH1F(Form("hllRefitMass_%d", ij), "", 750, 0., 15.);
    hllRefitMass[ij]->Sumw2();
    hllRefitMass[ij]->SetLineColor(kRed);
    hllRefitMass[ij]->SetLineWidth(2);

    hllPrefitMass[ij] = new TH1F(Form("hllPrefitMass_%d", ij), "", 750, 0., 15.);
    hllPrefitMass[ij]->Sumw2();
    hllPrefitMass[ij]->SetLineColor(kRed);
    hllPrefitMass[ij]->SetLineWidth(2);

    hllMass[ij] = new TH1F(Form("hllMass_%d", ij), "", 750, 0., 15.);
    hllMass[ij]->Sumw2();
    hllMass[ij]->SetLineColor(kRed);
    hllMass[ij]->SetLineWidth(2);

    hllMass_vs_Bmass[ij] = new TH2F(Form("hllMass_vs_Bmass_%d", ij), "", 500, 0., 15., 500, 0., 15.);
    hllMass_vs_Bmass[ij]->Sumw2();
    hllMass_vs_Bmass[ij]->SetMarkerColor(kRed);
    hllMass_vs_Bmass[ij]->SetMarkerStyle(20);

    hBmass[ij] = new TH1F(Form("Bmass_%d", ij), "", 750, 0., 15.); // 75, 4.5, 6.);
    hBmass[ij]->Sumw2();
    hBmass[ij]->SetLineColor(kRed);
    hBmass[ij]->SetLineWidth(2);

    hctxy[ij] = new TH1F(Form("hctxy_%d", ij), "", 1000, 0., 10.);
    hctxy[ij]->Sumw2();
    hctxy[ij]->SetLineColor(kRed);
    hctxy[ij]->SetLineWidth(2);

    hKaonpt[ij] = new TH1F(Form("hKaonpt_%d", ij), "", 100, 0., 10.);
    hKaonpt[ij]->Sumw2();
    hKaonpt[ij]->SetLineColor(kRed);
    hKaonpt[ij]->SetLineWidth(2);

    hLep1pt[ij] = new TH1F(Form("hLep1pt_%d", ij), "", 100, 0., 10.);
    hLep1pt[ij]->Sumw2();
    hLep1pt[ij]->SetLineColor(kRed);
    hLep1pt[ij]->SetLineWidth(2);

    hLep2pt[ij] = new TH1F(Form("hLep2pt_%d", ij), "", 100, 0., 10.);
    hLep2pt[ij]->Sumw2();
    hLep2pt[ij]->SetLineColor(kRed);
    hLep2pt[ij]->SetLineWidth(2);

    //
    hLep1pt_EB[ij] = new TH1F(Form("hLep1pt_EB_%d", ij), "", 100, 0., 10.);
    hLep1pt_EB[ij]->Sumw2();
    hLep1pt_EB[ij]->SetLineColor(kRed);
    hLep1pt_EB[ij]->SetLineWidth(2);

    hLep2pt_EB[ij] = new TH1F(Form("hLep2pt_EB_%d", ij), "", 100, 0., 10.);
    hLep2pt_EB[ij]->Sumw2();
    hLep2pt_EB[ij]->SetLineColor(kRed);
    hLep2pt_EB[ij]->SetLineWidth(2);

    //
    hLep1pt_EE[ij] = new TH1F(Form("hLep1pt_EE_%d", ij), "", 100, 0., 10.);
    hLep1pt_EE[ij]->Sumw2();
    hLep1pt_EE[ij]->SetLineColor(kRed);
    hLep1pt_EE[ij]->SetLineWidth(2);

    hLep2pt_EE[ij] = new TH1F(Form("hLep2pt_EE_%d", ij), "", 100, 0., 10.);
    hLep2pt_EE[ij]->Sumw2();
    hLep2pt_EE[ij]->SetLineColor(kRed);
    hLep2pt_EE[ij]->SetLineWidth(2);

    hBpt[ij] = new TH1F(Form("hBpt_%d", ij), "", 200, 0., 100.);
    hBpt[ij]->Sumw2();
    hBpt[ij]->SetLineColor(kRed);
    hBpt[ij]->SetLineWidth(2);
  }


  float nEv_muonTag[1] = {0.};
  float nEv_recoCand[1] = {0.};
  float nEv_chargeSel[1] = {0.};
  // float nEv_alphaEff[1] = {0.};
  // float nEv_vtxCLEff[1] = {0.};
  // float nEv_LxyEff[1] = {0.};
  float nEv_selected[7] = {0.};


  if(nMaxEvents == -1) nMaxEvents = nEvts;
  for(int iEvt = 0; iEvt<nMaxEvents; ++iEvt){
    if(iEvt%500000 == 0) std::cout << " >>> processing event " << iEvt << " " << 1.*iEvt/nEvts*100. << std::endl;

    t1->GetEntry(iEvt);

    if(MuonTag_index == -1) continue;
    if(!isEleFinalState && MuonRecoTag_index == -1) continue;
    ++nEv_muonTag[0];

    if(BToKll_sel_index == -1) continue; 
    if(dataset == "MC" && BPHRun == "nnResonant" && BToKll_sel_index != BToKll_gen_index) continue;
    ++nEv_recoCand[0];


    // if(!isEleFinalState && (Lepton_softId[BToKll_lep1_index[BToKll_sel_index]] != 1 ||
    // 			    Lepton_softId[BToKll_lep2_index[BToKll_sel_index]] != 1)) continue;
    // if(!isEleFinalState && (Lepton_mediumId[BToKll_lep1_index[BToKll_sel_index]] != 1 ||
    // 			    Lepton_mediumId[BToKll_lep2_index[BToKll_sel_index]] != 1)) continue;
 

    //opposite sign leptons
    if(BToKll_lep1_charge[BToKll_sel_index]*BToKll_lep2_charge[BToKll_sel_index] > 0.) continue;
    ++nEv_chargeSel[0];

    float llInvPrefitMass = 0.;
    if(isEleFinalState){
      TLorentzVector ele1cand;
      ele1cand.SetPtEtaPhiM(Lepton_pt[BToKll_lep1_index[BToKll_sel_index]], Lepton_eta[BToKll_lep1_index[BToKll_sel_index]], 
			    Lepton_phi[BToKll_lep1_index[BToKll_sel_index]], Lepton_mass[BToKll_lep1_index[BToKll_sel_index]]);
      TLorentzVector ele2cand;
      ele2cand.SetPtEtaPhiM(Lepton_pt[BToKll_lep2_index[BToKll_sel_index]], Lepton_eta[BToKll_lep2_index[BToKll_sel_index]], 
			    Lepton_phi[BToKll_lep2_index[BToKll_sel_index]], Lepton_mass[BToKll_lep2_index[BToKll_sel_index]]);
      /*
      std::cout << " ele1 pt = " << Lepton_pt[BToKll_lep1_index[BToKll_sel_index]] 
		<< " ele2 pt = " << Lepton_pt[BToKll_lep2_index[BToKll_sel_index]]
		<< " ele1 eta = " << Lepton_eta[BToKll_lep1_index[BToKll_sel_index]] 
		<< " ele2 eta = " << Lepton_eta[BToKll_lep2_index[BToKll_sel_index]] 
		<< " ele1 mass  = " << Lepton_mass[BToKll_lep1_index[BToKll_sel_index]]
		<< " ele2 mass = " << Lepton_mass[BToKll_lep2_index[BToKll_sel_index]] << std::endl;
      */

      llInvPrefitMass = (ele1cand+ele2cand).Mag();
    }
    //    std::cout << " llInvPrefitMass = " << llInvPrefitMass << std::endl;

    float llInvRefitMass = BToKll_llKFit_ll_mass[BToKll_sel_index];
    float llInvMass = BToKll_ll_mass[BToKll_sel_index];
    int massBin = -1;
    float massVar = llInvRefitMass;
    if(isEleFinalState) massVar = llInvPrefitMass;
    for(unsigned int kl=0; kl<llMassBoundary.size()-1; ++kl){
      if(massVar >= llMassBoundary[kl] && massVar < llMassBoundary[kl+1]){
	massBin = kl;
	break;
      }
    }

    //to synch. with Riccardo ~ tight selection
    if(typeSelection == "tightCB"){
      if((BToKll_kaon_pt[BToKll_sel_index] < 1.5 || BToKll_pt[BToKll_sel_index] < 10.)) continue;
      if(BToKll_cosAlpha[BToKll_sel_index] < 0.999) continue;
      if(BToKll_CL_vtx[BToKll_sel_index] < 0.1) continue;
      if(BToKll_Lxy[BToKll_sel_index] < 6.) continue;
    }

    // MVA selection to get same background rejection as with the cut base
    if(typeSelection == "NN_BkgR" && nnBMX < 0.993) continue;

    // MVA selection to get same signal efficiency as with the cut base
    if(typeSelection == "NN_SigEff" && nnBMX < 0.998) continue;


    if(massBin != -1) ++nEv_selected[massBin];
    else ++nEv_selected[6];

    //save output ntuple with selected events fro final plots
      if(saveOUTntu){
	Kll_cosAlpha = BToKll_cosAlpha[BToKll_sel_index];
	Kll_CL_vtx = BToKll_CL_vtx[BToKll_sel_index];
	Kll_kaon_DCASig = BToKll_kaon_DCASig[BToKll_sel_index];
	Kll_Lxy = BToKll_Lxy[BToKll_sel_index];
	Kll_ctxy = BToKll_ctxy[BToKll_sel_index];
	Kll_llmass = massVar;
	Kll_llRefitmass = llInvRefitMass;
	Kll_llPrefitmass = llInvPrefitMass;
	Kll_mass = BToKll_mass[BToKll_sel_index];
	Kll_pt = BToKll_pt[BToKll_sel_index];
	Kll_kaon_pt = BToKll_kaon_pt[BToKll_sel_index];
	Kll_kaon_eta = BToKll_kaon_eta[BToKll_sel_index];
	Kll_kaon_phi = BToKll_kaon_phi[BToKll_sel_index];
	Kll_kaon_charge = BToKll_kaon_charge[BToKll_sel_index];
	Kll_lep1_pt = BToKll_lep1_pt[BToKll_sel_index];
	Kll_lep2_pt = BToKll_lep2_pt[BToKll_sel_index];
	Kll_lep1_eta = BToKll_lep1_eta[BToKll_sel_index];
	Kll_lep2_eta = BToKll_lep2_eta[BToKll_sel_index];
	Kll_lep1_phi = BToKll_lep1_phi[BToKll_sel_index];
	Kll_lep2_phi = BToKll_lep2_phi[BToKll_sel_index];
	Kll_lep1_preFpt = Lepton_pt[BToKll_lep1_index[BToKll_sel_index]];
	Kll_lep2_preFpt = Lepton_pt[BToKll_lep2_index[BToKll_sel_index]];
	Kll_lep1_preFeta = Lepton_eta[BToKll_lep1_index[BToKll_sel_index]];
	Kll_lep2_preFeta = Lepton_eta[BToKll_lep2_index[BToKll_sel_index]];
	Kll_lep1_preFphi = Lepton_phi[BToKll_lep1_index[BToKll_sel_index]];
	Kll_lep2_preFphi = Lepton_phi[BToKll_lep2_index[BToKll_sel_index]];

	Kll_lep1_pfRelIso03 = Lepton_pfRelIso03[BToKll_lep1_index[BToKll_sel_index]];
	Kll_lep2_pfRelIso03 = Lepton_pfRelIso03[BToKll_lep2_index[BToKll_sel_index]];
	if(isEleFinalState){
	  Kll_lep1_pfRelIso04 = -1.;
	  Kll_lep2_pfRelIso04 = -1.;
	}
	else{
	Kll_lep1_pfRelIso04 = Lepton_pfRelIso04[BToKll_lep1_index[BToKll_sel_index]];
	Kll_lep2_pfRelIso04 = Lepton_pfRelIso04[BToKll_lep2_index[BToKll_sel_index]];
	}
	Kll_lep1_dxy = Lepton_dxy[BToKll_lep1_index[BToKll_sel_index]];
	Kll_lep2_dxy = Lepton_dxy[BToKll_lep2_index[BToKll_sel_index]];
	Kll_lep1_dz = Lepton_dz[BToKll_lep1_index[BToKll_sel_index]];
	Kll_lep2_dz = Lepton_dz[BToKll_lep2_index[BToKll_sel_index]];
	Kll_kaon_dxy = PFCand_dxy[BToKll_kaon_index[BToKll_sel_index]];
	Kll_kaon_dz = PFCand_dz[BToKll_kaon_index[BToKll_sel_index]];

	newT->Fill();
      }
      //end output Nutuple


    //histograms for each mass bin
    if(massBin != -1){
      hAlpha[massBin]->Fill(BToKll_cosAlpha[BToKll_sel_index]);
      hCLVtx[massBin]->Fill(BToKll_CL_vtx[BToKll_sel_index]);
      hDCASig[massBin]->Fill(BToKll_kaon_DCASig[BToKll_sel_index]);
      hLxy[massBin]->Fill(BToKll_Lxy[BToKll_sel_index]);
      hctxy[massBin]->Fill(BToKll_ctxy[BToKll_sel_index]);
      hKaonpt[massBin]->Fill(BToKll_kaon_pt[BToKll_sel_index]);
      hBpt[massBin]->Fill(BToKll_pt[BToKll_sel_index]);

      hLep1pt[massBin]->Fill(BToKll_lep1_pt[BToKll_sel_index]);
      hLep2pt[massBin]->Fill(BToKll_lep2_pt[BToKll_sel_index]);
      if(std::abs(BToKll_lep1_eta[BToKll_sel_index]) < 1.47) hLep1pt_EB[massBin]->Fill(BToKll_lep1_pt[BToKll_sel_index]);
      else hLep1pt_EE[massBin]->Fill(BToKll_lep1_pt[BToKll_sel_index]);
      if(std::abs(BToKll_lep2_eta[BToKll_sel_index]) < 1.47) hLep2pt_EB[massBin]->Fill(BToKll_lep2_pt[BToKll_sel_index]);
      else hLep2pt_EE[massBin]->Fill(BToKll_lep2_pt[BToKll_sel_index]);

      hllPrefitMass[massBin]->Fill(llInvPrefitMass);
      hllRefitMass[massBin]->Fill(llInvRefitMass);
      hllMass[massBin]->Fill(llInvMass);
      hllMass_vs_Bmass[massBin]->Fill(BToKll_mass[BToKll_sel_index], llInvPrefitMass);
      hBmass[massBin]->Fill(BToKll_mass[BToKll_sel_index]);
    }
    
    //histograms inclusive over all m(ll)
    hllPrefitMass[6]->Fill(llInvPrefitMass);
    hllRefitMass[6]->Fill(llInvRefitMass);
    hllMass[6]->Fill(llInvMass);
    hllMass_vs_Bmass[6]->Fill(BToKll_mass[BToKll_sel_index], llInvPrefitMass);
    hBmass[6]->Fill(BToKll_mass[BToKll_sel_index]);
    
    hAlpha[6]->Fill(BToKll_cosAlpha[BToKll_sel_index]);
    hCLVtx[6]->Fill(BToKll_CL_vtx[BToKll_sel_index]);
    hDCASig[6]->Fill(BToKll_kaon_DCASig[BToKll_sel_index]);
    hLxy[6]->Fill(BToKll_Lxy[BToKll_sel_index]);
    hctxy[6]->Fill(BToKll_ctxy[BToKll_sel_index]);
    hKaonpt[6]->Fill(BToKll_kaon_pt[BToKll_sel_index]);
    hBpt[6]->Fill(BToKll_pt[BToKll_sel_index]);

    hLep1pt[6]->Fill(BToKll_lep1_pt[BToKll_sel_index]);
    hLep2pt[6]->Fill(BToKll_lep2_pt[BToKll_sel_index]);
    if(std::abs(BToKll_lep1_eta[BToKll_sel_index]) < 1.47) hLep1pt_EB[6]->Fill(BToKll_lep1_pt[BToKll_sel_index]);
    else hLep1pt_EE[6]->Fill(BToKll_lep1_pt[BToKll_sel_index]);
    if(std::abs(BToKll_lep2_eta[BToKll_sel_index]) < 1.47) hLep2pt_EB[6]->Fill(BToKll_lep2_pt[BToKll_sel_index]);
    else hLep2pt_EE[6]->Fill(BToKll_lep2_pt[BToKll_sel_index]);

    /*
    if(massBin == 3 && BToKll_mass[BToKll_sel_index] < 6.)
      std::cout << run << " " << lumi << " " << event << std::endl;
    */
  }//loop over events

  
  if(saveOUTntu){
    newFile->cd();
    newT->Write();
    newFile->Close(); 
  }

  outMassHistos.cd();
  
  std::cout << " ***** summary ***** "<< std::endl;
  for(int ij=0; ij<7; ++ij){
    hAlpha[ij]->Write(hAlpha[ij]->GetName());
    hCLVtx[ij]->Write(hCLVtx[ij]->GetName());
    hDCASig[ij]->Write(hDCASig[ij]->GetName());
    hLxy[ij]->Write(hLxy[ij]->GetName());
    hctxy[ij]->Write(hctxy[ij]->GetName());
    hKaonpt[ij]->Write(hKaonpt[ij]->GetName());
    hLep1pt[ij]->Write(hLep1pt[ij]->GetName());
    hLep2pt[ij]->Write(hLep2pt[ij]->GetName());
    hLep1pt_EB[ij]->Write(hLep1pt_EB[ij]->GetName());
    hLep1pt_EE[ij]->Write(hLep1pt_EE[ij]->GetName());
    hLep2pt_EB[ij]->Write(hLep2pt_EB[ij]->GetName());
    hLep2pt_EE[ij]->Write(hLep2pt_EE[ij]->GetName());

    hBpt[ij]->Write(hBpt[ij]->GetName());
    std::cout << " >>> hBpt[ij]->GetName() = " << hBpt[ij]->GetName() << " entries = " << hBpt[ij]->GetEntries() << std::endl;

    hllMass[ij]->Write(hllMass[ij]->GetName());
    hllRefitMass[ij]->Write(hllRefitMass[ij]->GetName());
    hllPrefitMass[ij]->Write(hllPrefitMass[ij]->GetName());
    hllMass_vs_Bmass[ij]->Write(hllMass_vs_Bmass[ij]->GetName());
    hBmass[ij]->Write(hBmass[ij]->GetName());

    if(ij > 5) continue;
    std::cout << "\n massBin: " << llMassBoundary[ij] << " - " << llMassBoundary[ij+1]
	      << " \n \t muonTag = " << nEv_muonTag[0] << " \t recoEvts = " << nEv_recoCand[0]
              << " \n \t chargeSel = " << nEv_chargeSel[0] << " \t selected Events = " << nEv_selected[ij] << std::endl;

  }
  outMassHistos.Close();

}
