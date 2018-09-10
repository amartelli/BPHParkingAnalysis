//g++ -Wall -o analyzeCharged_fastDATA `root-config --cflags --glibs` -lRooFitCore analyzeCharged_fastDATA.cpp
//
//to run on ele all dataset 
// ./analyzeCharged_fastDATA 1 -1
//to run on muon runA
// ./analyzeCharged_fastDATA 0  runA
// options are: isEleFinalState (1, 0)   dataset (-1, runA, runB, MC)  run (1,2,3,... ) nMaxEvents (-1, N)   saveOUTntu

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

const float ElectronMass = 0.5109989e-3;

int main(int argc, char *argv[]){

  if(argc < 2) {
    std::cout << " Missing arguments " << std::endl;
    return -1;
  }
  int isEleFinalState = atoi(argv[1]);
  std::string dataset = "-1";
  std::string BPHRun = "-1";
  int nMaxEvents = -1;
  int saveOUTntu = 0;
  if(argc > 2) dataset = argv[2];
  if(argc > 3) BPHRun = argv[3];
  if(argc > 4) nMaxEvents = atoi(argv[4]);
  if(argc > 5) saveOUTntu = atoi(argv[5]);


  std::cout << " isEleFinalState = " << isEleFinalState << " dataset = " << dataset << " BPHRun = " << BPHRun << " nMaxEvents = " << nMaxEvents << " saveOUTntu = " << saveOUTntu << std::endl;


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
    if((dataset == "runA" && BPHRun == "1") || dataset == "-1") t1->Add("/vols/cms/tstreble/BPH/BToKmumu_ntuple/BPHParking1_2018A_18_08_14_new/*root");
    if((dataset == "runA" && BPHRun == "2") || dataset == "-1") t1->Add("/vols/cms/tstreble/BPH/BToKmumu_ntuple/BPHParking2_2018A_18_08_14_new/*root");
    if((dataset == "runA" && BPHRun == "3") || dataset == "-1") t1->Add("/vols/cms/tstreble/BPH/BToKmumu_ntuple/BPHParking3_2018A_18_08_14_new/*root");
    if((dataset == "runA" && BPHRun == "4") || dataset == "-1") t1->Add("/vols/cms/tstreble/BPH/BToKmumu_ntuple/BPHParking4_2018A_18_08_14_new/*root");
    if((dataset == "runA" && BPHRun == "5") || dataset == "-1") t1->Add("/vols/cms/tstreble/BPH/BToKmumu_ntuple/BPHParking5_2018A_18_08_14_new/*root");
    if((dataset == "runA" && BPHRun == "6") || dataset == "-1") t1->Add("/vols/cms/tstreble/BPH/BToKmumu_ntuple/BPHParking6_2018A_18_08_14_new/*root");

    if((dataset == "runB" && BPHRun == "1") || dataset == "-1") t1->Add("/vols/cms/tstreble/BPH/BToKmumu_ntuple/BPHParking1_2018B_18_08_14_new/*root");
    if((dataset == "runB" && BPHRun == "2") || dataset == "-1") t1->Add("/vols/cms/tstreble/BPH/BToKmumu_ntuple/BPHParking2_2018B_18_08_14_new/*root");
    if((dataset == "runB" && BPHRun == "3") || dataset == "-1") t1->Add("/vols/cms/tstreble/BPH/BToKmumu_ntuple/BPHParking3_2018B_18_08_14_new/*root");
    if((dataset == "runB" && BPHRun == "4") || dataset == "-1") t1->Add("/vols/cms/tstreble/BPH/BToKmumu_ntuple/BPHParking4_2018B_18_08_14_new/*root");
    if((dataset == "runB" && BPHRun == "5") || dataset == "-1") t1->Add("/vols/cms/tstreble/BPH/BToKmumu_ntuple/BPHParking5_2018B_18_08_14_new/*root");

    if(dataset == "MC" && BPHRun == "nnResonant"){
      t1->Add("/vols/cms/tstreble/BPH/BToKmumu_ntuple/BToKmumu_18_08_14_new/*root");
    }
    if(dataset == "MC" && BPHRun == "Resonant"){
      t1->Add("/vols/cms/tstreble/BPH/BToKmumu_ntuple/BToKJPsimumu_18_08_14_new/*root");
    }
  }

  int nEvts = t1->GetEntries();
  std::cout << " #initial n. events: " << nEvts << std::endl;

  std::string outNtuName = "/vols/cms/amartell/BParking/data_allStat/selectedEvents_Kee_"+dataset+"_BPHRun"+BPHRun+".root";
  if(!isEleFinalState) outNtuName = "/vols/cms/amartell/BParking/data_allStat/selectedEvents_Kmumu_"+dataset+"_BPHRun"+BPHRun+".root";
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
  }

  std::cout << " >>> qui ok " << std::endl;

  UInt_t run = 0;
  UInt_t lumi = 0;
  ULong64_t event = 0;

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

  
  t1->SetBranchStatus("*", 0);


  t1->SetBranchStatus("run", 1);                        t1->SetBranchAddress("run", &run);
  t1->SetBranchStatus("luminosityBlock", 1);            t1->SetBranchAddress("luminosityBlock", &lumi);
  t1->SetBranchStatus("event", 1);                      t1->SetBranchAddress("event", &event);

  if(isEleFinalState){
    t1->SetBranchStatus("Muon_sel_index", 1);            t1->SetBranchAddress("Muon_sel_index", &MuonTag_index);
    t1->SetBranchStatus("BToKee_sel_index", 1);          t1->SetBranchAddress("BToKee_sel_index", &BToKll_sel_index);
    if(dataset == "MC"){
      t1->SetBranchStatus("BToKee_gen_index", 1);          t1->SetBranchAddress("BToKee_gen_index", &BToKll_gen_index);
    }
    t1->SetBranchStatus("BToKee_ele1_charge", 1);        t1->SetBranchAddress("BToKee_ele1_charge", &BToKll_lep1_charge);
    t1->SetBranchStatus("BToKee_ele2_charge", 1);        t1->SetBranchAddress("BToKee_ele2_charge", &BToKll_lep2_charge);
    t1->SetBranchStatus("BToKee_ele1_index", 1);         t1->SetBranchAddress("BToKee_ele1_index", &BToKll_lep1_index);
    t1->SetBranchStatus("BToKee_ele2_index", 1);         t1->SetBranchAddress("BToKee_ele2_index", &BToKll_lep2_index);
    t1->SetBranchStatus("BToKee_cosAlpha", 1);           t1->SetBranchAddress("BToKee_cosAlpha", &BToKll_cosAlpha);
    t1->SetBranchStatus("BToKee_CL_vtx", 1);             t1->SetBranchAddress("BToKee_CL_vtx", &BToKll_CL_vtx);
    t1->SetBranchStatus("BToKee_kaon_DCASig", 1);        t1->SetBranchAddress("BToKee_kaon_DCASig", &BToKll_kaon_DCASig);
    t1->SetBranchStatus("BToKee_Lxy", 1);                t1->SetBranchAddress("BToKee_Lxy", &BToKll_Lxy);
    t1->SetBranchStatus("BToKee_ctxy", 1);                t1->SetBranchAddress("BToKee_ctxy", &BToKll_ctxy);
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
    t1->SetBranchStatus("Electron_eta", 1);                t1->SetBranchAddress("Electron_eta", &Lepton_eta);
    t1->SetBranchStatus("Electron_phi", 1);                t1->SetBranchAddress("Electron_phi", &Lepton_phi);
    t1->SetBranchStatus("Electron_mass", 1);                t1->SetBranchAddress("Electron_mass", &Lepton_mass);
    t1->SetBranchStatus("Electron_charge", 1);             t1->SetBranchAddress("Electron_charge", &Lepton_charge);
    t1->SetBranchStatus("Electron_dxy", 1);                t1->SetBranchAddress("Electron_dxy", &Lepton_dxy);
    t1->SetBranchStatus("Electron_dz", 1);                t1->SetBranchAddress("Electron_dz", &Lepton_dz);

  }
  else{
    t1->SetBranchStatus("Muon_probe_index", 1);            t1->SetBranchAddress("Muon_probe_index", &MuonTag_index);
    t1->SetBranchStatus("Muon_sel_index", 1);              t1->SetBranchAddress("Muon_sel_index", &MuonRecoTag_index);
    t1->SetBranchStatus("Muon_softId", 1);                 t1->SetBranchAddress("Muon_softId", &Lepton_softId);
    t1->SetBranchStatus("Muon_mediumId", 1);               t1->SetBranchAddress("Muon_mediumId", &Lepton_mediumId);
    t1->SetBranchStatus("BToKmumu_sel_index", 1);          t1->SetBranchAddress("BToKmumu_sel_index", &BToKll_sel_index);
    if(dataset == "MC"){
      t1->SetBranchStatus("BToKmumu_gen_index", 1);          t1->SetBranchAddress("BToKmumu_gen_index", &BToKll_gen_index);
    }
    t1->SetBranchStatus("BToKmumu_mu1_charge", 1);         t1->SetBranchAddress("BToKmumu_mu1_charge", &BToKll_lep1_charge);
    t1->SetBranchStatus("BToKmumu_mu2_charge", 1);         t1->SetBranchAddress("BToKmumu_mu2_charge", &BToKll_lep2_charge);
    t1->SetBranchStatus("BToKmumu_mu1_index", 1);         t1->SetBranchAddress("BToKmumu_mu1_index", &BToKll_lep1_index);
    t1->SetBranchStatus("BToKmumu_mu2_index", 1);         t1->SetBranchAddress("BToKmumu_mu2_index", &BToKll_lep2_index);
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

  }



  std::vector<float> llMassBoundary;
  llMassBoundary.push_back(0.);
  llMassBoundary.push_back(1.);
  llMassBoundary.push_back(2.5);
  llMassBoundary.push_back(2.9);
  llMassBoundary.push_back(3.2);
  llMassBoundary.push_back(3.58);



  std::string outName = "outMassHistos_Kee_"+dataset+"_BPHRun"+BPHRun+".root";
  if(!isEleFinalState) outName = "outMassHistos_Kmumu_"+dataset+"_BPHRun"+BPHRun+".root";
  TFile outMassHistos(outName.c_str(), "recreate");

  ///histos: 1 per bin plus inclusive
  TH1F* hAlpha[6];
  TH1F* hCLVtx[6];
  TH1F* hDCASig[6];
  TH1F* hLxy[6];
  TH1F* hctxy[6];
  TH1F* hKaonpt[6];
  TH1F* hLep1pt[6];
  TH1F* hLep2pt[6];
  TH1F* hLep1pt_EB[6];
  TH1F* hLep2pt_EB[6];
  TH1F* hLep1pt_EE[6];
  TH1F* hLep2pt_EE[6];
  TH1F* hBpt[6];
  TH1F* hllMass[6];
  TH1F* hllRefitMass[6];
  TH1F* hllPrefitMass[6];
  TH2F* hllMass_vs_Bmass[6];
  TH1F* hBmass[6];

  for(int ij=0; ij<6; ++ij){
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


  float nEv_muonTag[5] = {0.};
  float nEv_recoCand[5] = {0.};
  float nEv_chargeEff[5] = {0.};
  float nEv_alphaEff[5] = {0.};
  float nEv_vtxCLEff[5] = {0.};
  float nEv_LxyEff[5] = {0.};

  int nEvents_hm_5p6 = 0;
  int nEvents_hm_5p7 = 0;

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
 

    if(BToKll_lep1_charge[BToKll_sel_index]*BToKll_lep2_charge[BToKll_sel_index] > 0.) continue;


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



    if(dataset == "MC" || BToKll_mass[BToKll_sel_index] > 5.7){
      ++nEvents_hm_5p7;

      if(saveOUTntu){
	Kll_cosAlpha = BToKll_cosAlpha[BToKll_sel_index];
	Kll_CL_vtx = BToKll_CL_vtx[BToKll_sel_index];
	Kll_kaon_DCASig = BToKll_kaon_DCASig[BToKll_sel_index];
	Kll_Lxy = BToKll_Lxy[BToKll_sel_index];
	Kll_ctxy = BToKll_ctxy[BToKll_sel_index];
	Kll_llmass = BToKll_ll_mass[BToKll_sel_index];
	Kll_llRefitmass = llInvRefitMass;
	Kll_llPrefitmass = llInvPrefitMass;
	Kll_mass = llInvMass;
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

	newT->Fill();
      }
    }

    if(massBin != -1) ++nEv_chargeEff[massBin];

    if((BToKll_kaon_pt[BToKll_sel_index] < 1.5 || BToKll_pt[BToKll_sel_index] < 10.)) continue;

    if(BToKll_cosAlpha[BToKll_sel_index] < 0.999) continue;
    //if(BToKll_cosAlpha[BToKll_sel_index] < 0.99) continue;

    if(massBin != -1) ++nEv_alphaEff[massBin];
    if(BToKll_CL_vtx[BToKll_sel_index] < 0.1) continue;
    if(massBin != -1) ++nEv_vtxCLEff[massBin];    
    if(BToKll_Lxy[BToKll_sel_index] < 6.) continue;
    if(massBin != -1) ++nEv_LxyEff[massBin];



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
    
    hllPrefitMass[5]->Fill(llInvPrefitMass);
    hllRefitMass[5]->Fill(llInvRefitMass);
    hllMass[5]->Fill(llInvMass);
    hllMass_vs_Bmass[5]->Fill(BToKll_mass[BToKll_sel_index], llInvPrefitMass);
    hBmass[5]->Fill(BToKll_mass[BToKll_sel_index]);
    
    hAlpha[5]->Fill(BToKll_cosAlpha[BToKll_sel_index]);
    hCLVtx[5]->Fill(BToKll_CL_vtx[BToKll_sel_index]);
    hDCASig[5]->Fill(BToKll_kaon_DCASig[BToKll_sel_index]);
    hLxy[5]->Fill(BToKll_Lxy[BToKll_sel_index]);
    hctxy[5]->Fill(BToKll_ctxy[BToKll_sel_index]);
    hKaonpt[5]->Fill(BToKll_kaon_pt[BToKll_sel_index]);
    hBpt[5]->Fill(BToKll_pt[BToKll_sel_index]);

    hLep1pt[5]->Fill(BToKll_lep1_pt[BToKll_sel_index]);
    hLep2pt[5]->Fill(BToKll_lep2_pt[BToKll_sel_index]);      
    if(std::abs(BToKll_lep1_eta[BToKll_sel_index]) < 1.47) hLep1pt_EB[5]->Fill(BToKll_lep1_pt[BToKll_sel_index]);
    else hLep1pt_EE[5]->Fill(BToKll_lep1_pt[BToKll_sel_index]);
    if(std::abs(BToKll_lep2_eta[BToKll_sel_index]) < 1.47) hLep2pt_EB[5]->Fill(BToKll_lep2_pt[BToKll_sel_index]);
    else hLep2pt_EE[5]->Fill(BToKll_lep2_pt[BToKll_sel_index]);

    if(massBin == 3 && BToKll_mass[BToKll_sel_index] < 6.)
      std::cout << run << " " << lumi << " " << event << std::endl;
  }//loop over events

  
  if(saveOUTntu){
    newFile->cd();
    newT->Write();
    newFile->Close(); 
  }

  
  // std::string outName = "outMassHistos_Kee_"+dataset+"_BPHRun"+BPHRun+".root";
  // if(!isEleFinalState) outName = "outMassHistos_Kmumu_"+dataset+"_BPHRun"+BPHRun+".root";
  // TFile outMassHistos(outName.c_str(), "recreate");
  outMassHistos.cd();
  
  std::cout << " ***** summary ***** "<< std::endl;
  for(int ij=0; ij<6; ++ij){
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

    if(ij > 4) continue;
    std::cout << "\n massBin: " << llMassBoundary[ij] << " - " << llMassBoundary[ij+1]
              << " \n \t recoEvts = " << nEv_chargeEff[ij] << " \t alphaCut = " << nEv_alphaEff[ij]
	      << " \t vtxCLCut = " << nEv_vtxCLEff[ij] << " \t LxyCut = " << nEv_LxyEff[ij] << std::endl;

  }
  outMassHistos.Close();

  return 100;

  TLegend *legM = new TLegend(0.70,0.70,0.98,0.95,NULL,"brNDC");
  legM->SetTextFont(42);
  legM->SetFillColor(kWhite);
  legM->SetLineColor(kWhite);
  legM->SetShadowColor(kWhite);
  legM->SetFillStyle(0);
  legM->SetTextSize(0.05);
  legM->AddEntry(hAlpha[3], "DATA", "l");
  

  TCanvas* c1 = new TCanvas();
  c1->cd();
  gPad->SetLogy();
  hAlpha[3]->GetXaxis()->SetTitle("cosAlpha");
  hAlpha[3]->Draw("hist");
  legM->Draw("same");
  if(isEleFinalState){
    c1->Print("plots/Kee_cosAlpha_DA.png", "png");
    c1->Print("plots/Kee_cosAlpha_DA.pdf", "pdf");
    c1->Print("plots/Kee_cosAlpha_DA.root", "root");
  }
  else{
    c1->Print("plots/Kmumu_cosAlpha_DA.png", "png");
    c1->Print("plots/Kmumu_cosAlpha_DA.pdf", "pdf");
    c1->Print("plots/Kmumu_cosAlpha_DA.root", "root");
  }

  TCanvas* c1b = new TCanvas();
  c1b->cd();
  gPad->SetLogy();
  hCLVtx[3]->GetXaxis()->SetTitle("CL_vtx");
  hCLVtx[3]->Draw("hist");
  legM->Draw("same");
  if(isEleFinalState){
    c1b->Print("plots/Kee_CL_vtx_DA.png", "png");
    c1b->Print("plots/Kee_CL_vtx_DA.pdf", "pdf");
    c1b->Print("plots/Kee_CL_vtx_DA.root", "root");
  }
  else{
    c1b->Print("plots/Kmumu_CL_vtx_DA.png", "png");
    c1b->Print("plots/Kmumu_CL_vtx_DA.pdf", "pdf");
    c1b->Print("plots/Kmumu_CL_vtx_DA.root", "root");
  }

  TCanvas* c1c = new TCanvas();
  hDCASig[3]->GetXaxis()->SetTitle("kaon DCA SIP");
  hDCASig[3]->Draw("hist");
  legM->Draw("same");
  if(isEleFinalState){
    c1c->Print("plots/Kee_kaonDCA_DA.png", "png");
    c1c->Print("plots/Kee_kaonDCA_DA.pdf", "pdf");
    c1c->Print("plots/Kee_kaonDCA_DA.root", "root");
  }
  else{
    c1c->Print("plots/Kmumu_kaonDCA_DA.png", "png");
    c1c->Print("plots/Kmumu_kaonDCA_DA.pdf", "pdf");
    c1c->Print("plots/Kmumu_kaonDCA_DA.root", "root");
  }
  
  TCanvas* c1d = new TCanvas();
  gPad->SetLogy();
  hLxy[3]->GetXaxis()->SetTitle("Lxy IP(cm)");
  hLxy[3]->Draw("hist");
  legM->Draw("same");
  if(isEleFinalState){
    c1d->Print("plots/Kee_Lxy_DA.png", "png");
    c1d->Print("plots/Kee_Lxy_DA.pdf", "pdf");
    c1d->Print("plots/Kee_Lxy_DA.root", "root");
  }
  else{
    c1d->Print("plots/Kmumu_Lxy_DA.png", "png");
    c1d->Print("plots/Kmumu_Lxy_DA.pdf", "pdf");
    c1d->Print("plots/Kmumu_Lxy_DA.root", "root");
  }

  TCanvas* c1e = new TCanvas();
  gPad->SetLogy();
  if(isEleFinalState)
    hllMass[3]->GetXaxis()->SetTitle("Mee mass");
  else hllMass[3]->GetXaxis()->SetTitle("Mmumu mass");
  hllMass[3]->Draw("hist");
  legM->Draw("same");
  if(isEleFinalState){
    c1e->Print("plots/Kee_Mee_DA.png", "png");
    c1e->Print("plots/Kee_Mee_DA.pdf", "pdf");
    c1e->Print("plots/Kee_Mee_DA.root", "root");
  }
  else{
    c1e->Print("plots/Kmumu_Mee_DA.png", "png");
    c1e->Print("plots/Kmumu_Mee_DA.pdf", "pdf");
    c1e->Print("plots/Kmumu_Mee_DA.root", "root");
  }

  TCanvas* c1f = new TCanvas();
  if(isEleFinalState){
    hllMass_vs_Bmass[3]->GetXaxis()->SetTitle("Kee mass");
    hllMass_vs_Bmass[3]->GetYaxis()->SetTitle("ee mass");
  }
  else{
    hllMass_vs_Bmass[3]->GetXaxis()->SetTitle("Kmumu mass");
    hllMass_vs_Bmass[3]->GetYaxis()->SetTitle("mumu mass");
  }
  hllMass_vs_Bmass[3]->Draw("");
  legM->Draw("same");
  if(isEleFinalState){
    c1f->Print("plots/Kee_Mee_vsBmass_DA.png", "png");
    c1f->Print("plots/Kee_Mee_vsBmass_DA.pdf", "pdf");
    c1f->Print("plots/Kee_Mee_vsBmass_DA.root", "root");
  }
  else{
    c1f->Print("plots/Kmumu_Mmumu_vsBmass_DA.png", "png");
    c1f->Print("plots/Kmumu_Mmumu_vsBmass_DA.pdf", "pdf");
    c1f->Print("plots/Kmumu_Mmumu_vsBmass_DA.root", "root");
  }

  TCanvas* c1g = new TCanvas();
  gPad->SetLogy();
  if(isEleFinalState)  hBmass[3]->GetXaxis()->SetTitle("B(Kee) mass");
  else hBmass[3]->GetXaxis()->SetTitle("B(Kmumu) mass");
  hBmass[3]->Draw("hist");
  legM->Draw("same");
  if(isEleFinalState){
    c1g->Print("plots/Kee_mass_DA.png", "png");
    c1g->Print("plots/Kee_mass_DA.pdf", "pdf");
    c1g->Print("plots/Kee_mass_DA.root", "root");
  }
  else{
    c1g->Print("plots/Kmumu_mass_DA.png", "png");
    c1g->Print("plots/Kmumu_mass_DA.pdf", "pdf");
    c1g->Print("plots/Kmumu_mass_DA.root", "root");
  }



  //comment out for the moment
  ///////////////////////////////////////////////////////////
  //now fitting

  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

  float nEv_postFit[5] = {0.};
  float nEvError_postFit[5] = {0.};

  RooWorkspace w("w");
  w.factory("x[0, 10]");
  //w.factory("x[4.5, 6.]");                                                                                                                                           

  w.factory("nbackground[10000, 0, 10000]");
  w.factory("nsignal[100, 0.0, 10000.0]");

  for(int ij=0; ij<5; ++ij){
    //w.factory("Gaussian::smodel(x,mu[5.3,4.5,6],sigma[0.05,0,0.2])");                                                                                                
    //w.factory("RooCBShape::smodel(x,m[5.3,4.5,6],s[0.1,0.,1.],a[1.2,0.,3.],n[1,0.1,6.])");
    w.factory("RooCBShape::smodel(x,m[5.3,0.,10.],s[0.1,0.,1.],a[1.2,0.,3.],n[1,0.1,6.])");
    //w.factory("RooCBShape::CBall(x[0,15], mean[11000,13000], sigma[5000,200000], alpha[0,10000],n[0,100000])");                                                      
    RooAbsPdf * smodel = w.pdf("smodel");

    w.factory("Exponential::bmodel(x,tau[-2,-3,0])");
    RooAbsPdf * bmodel = w.pdf("bmodel");

    w.factory("SUM::model(nbackground*bmodel, nsignal*smodel)");
    RooAbsPdf * model = w.pdf("model");

    RooDataHist hBMass("hBMass", "hBMass", *w.var("x"), Import(*(hBmass[ij])));

    RooFitResult * r = model->fitTo(hBMass, Minimizer("Minuit2"),Save(true));

    RooPlot * plot = w.var("x")->frame();
    if(isEleFinalState){
      plot->SetXTitle("Kee mass (GeV)");
    }
    else{
      plot->SetXTitle("K#mu#mu mass (GeV)");
    }
    plot->SetTitle("");
    plot->SetAxisRange(4.,6);
    hBMass.plotOn(plot);
    model->plotOn(plot);
    model->plotOn(plot, Components("bmodel"),LineStyle(kDashed));
    model->plotOn(plot, Components("smodel"),LineColor(kRed));

    TCanvas * cc = new TCanvas();
    cc->SetLogy(0);
    plot->Draw();
    if(isEleFinalState) cc->Print(Form("plots/Bmass_DATA/Kee_%s.png",hBmass[ij]->GetName()), "png");
    else cc->Print(Form("plots/Bmass_DATA/Kmumu_%s.png",hBmass[ij]->GetName()), "png");

    RooRealVar* parS = (RooRealVar*) r->floatParsFinal().find("nsignal");
    RooRealVar* parB = (RooRealVar*) r->floatParsFinal().find("nbackground");
    nEv_postFit[ij] = parS->getValV();
    nEvError_postFit[ij] = parS->getError();

    std::cout << " selection signal events = \t " << parS->getValV() << " error = " << parS->getError()
              << " bkg events = " << parB->getValV() << " error = " << parB->getError() << std::endl;
  }

  std::cout << " ***** summary ***** "<< std::endl;
  for(int ij=0; ij<5; ++ij){

    std::cout << "\n category = " << hBmass[ij]->GetName()
              << " \n \t recoEvts = " << nEv_LxyEff[ij] << " postFit " << nEv_postFit[ij] << "+/-" << nEvError_postFit[ij] << std::endl;
  }


  std::cout << " nEvents_hm_5p6 = " << nEvents_hm_5p6 << " nEvents_hm_5p7 " << nEvents_hm_5p7 << std::endl;

}
