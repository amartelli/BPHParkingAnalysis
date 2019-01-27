//compute double ratio nnReso/JPsi  event counts following all levels of selections
//muonTag + 
//gen Acceptance +
//gen Efficiency (reco final state matched to gen level) +
//selections for the analysis in data  

//should be ok for electron final state 
//example to run
//electron final state
//root -l computeAE_charged_fitBmass_Kstll.C'("/path-to-the-Kee-directory/ntu_MC_PR38_BToKee.root", "/path-to-the-KJPsiee-directory/ntu_MC_PR38_BToKJPsiee.root", 1, true)' 
//muon final state
//root -l computeAE_charged_fitBmass_Kstll.C'("/path-to-the-Kmumu-directory/ntu_MC_PR38_BToKmumu.root", "/path-to-the-KJPsimumu-directory/ntu_MC_PR38_BToKJPsimumu.root", 0, true)' 


#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>

#include "TROOT.h"
#include "TSystem.h"
#include "TKey.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TString.h"
#include "TCut.h"
#include "TMath.h"
#include "TApplication.h"
#include "TError.h"
#include "TGraphErrors.h"
#include "TPad.h"
#include "TMultiGraph.h"
#include "TStyle.h"
#include "TChain.h"

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooCBShape.h"
#include "RooAddPdf.h"
#include "RooWorkspace.h"
#include "RooFitResult.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"

using namespace RooFit;


void computeAE_charged_fitBmass_Kstll(std::string nonResonantFile, std::string ResonantFile, int isEleFinalState, bool foldGenMassBin=true){

  gROOT->Reset();
  gROOT->Macro("setStyle.C");

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);


  std::cout << " isEleFinalState = " << isEleFinalState << std::endl;

  TChain* t1 = new TChain("Events");
  TChain* t2 = new TChain("Events");

  t1->Add(nonResonantFile.c_str());
  t2->Add(ResonantFile.c_str());


  //muon tag with soft ID https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonIdRun2#Soft_Muon
  int nMuonTag_nonReso = t1->GetEntries();
  int nMuonTag_Reso = t2->GetEntries();

  std::cout << " #muon tag events: nnreso = " << nMuonTag_nonReso << " resonant = " << nMuonTag_Reso << std::endl;

  //if false => move to reco ee invariant mass bins

   //selections
  std::string cut_muonTag = "Muon_sel_index != -1";
  if(!isEleFinalState){
    cut_muonTag = "Muon_probe_index != -1";
  }
  std::string cut_Recocandidate = "BToKstll_sel_index != -1 && BToKstll_gen_index != -1 && BToKstll_sel_index == BToKstll_gen_index";
  std::string cut_chargeEff = "BToKstll_lep1_charge[BToKstll_sel_index]*BToKstll_lep2_charge[BToKstll_sel_index] < 0.";
  std::string cut_pTEff = "BToKstll_B_pt[BToKstll_sel_index] >= 10. && BToKstll_kaon_pt[BToKstll_sel_index] >= 1.5";
  std::string cut_alphaEff = "BToKstll_B_cosAlpha[BToKstll_sel_index] >= 0.999";
  std::string cut_vtxCLEff = "BToKstll_B_CL_vtx[BToKstll_sel_index] >= 0.1";  
  std::string cut_LxyEff  = "BToKstll_B_Lxy[BToKstll_sel_index] >= 6"; 

  std::vector<std::string> llMassCut;
  std::vector<std::string> llGenMassCut;

  std::vector<float> llMassBoundary;
  llMassBoundary.push_back(0.);
  llMassBoundary.push_back(1.);
  llMassBoundary.push_back(2.5);
  llMassBoundary.push_back(2.9);
  llMassBoundary.push_back(3.3);
  llMassBoundary.push_back(3.58);
  llMassBoundary.push_back(100.);

  
  for(int ij=0; ij<6; ++ij){
    
    std::string cut = Form("BToKstll_ll_mass[BToKstll_sel_index] > %.2f && BToKstll_ll_mass[BToKstll_sel_index] < %.2f",            
    llMassBoundary.at(ij), llMassBoundary.at(ij+1));
    
    std::string gencut = Form("BToKstll_gen_llMass > %.2f && BToKstll_gen_llMass < %.2f", llMassBoundary.at(ij), llMassBoundary.at(ij+1));
    
    llMassCut.push_back(cut);
    if(foldGenMassBin) llGenMassCut.push_back(gencut);
    else llGenMassCut.push_back(" 1 == 1");  
      
  }

  
  //JPsi first - non resonant last
  float nEv_muonTag[7] = {0.};
  float nEv_cuts[7] = {0.};
  

  TH1F* h_Bmass_llbin[7];
  h_Bmass_llbin[0] = new TH1F(Form("h_Bmass_JPsi_%.2f-%.2f", llMassBoundary.at(3), llMassBoundary.at(3+1)), "", 1000, 0., 10.);
    
  nEv_muonTag[0] = t2->Draw("Muon_sel_index", (cut_muonTag+" && "+llGenMassCut.at(3)).c_str(), "goff");
  
  nEv_cuts[0] = t2->Draw(Form("BToKstll_B_mass[BToKstll_sel_index] >> %s", h_Bmass_llbin[0]->GetName()), (cut_muonTag+" && "+cut_Recocandidate+" && "+cut_chargeEff+" && "+cut_pTEff+" && "+cut_alphaEff+" && "+cut_vtxCLEff+" && "+cut_LxyEff+" && "+llMassCut.at(3)+" && "+llGenMassCut.at(3)).c_str());
  
  std::cout << " tot JPsi_MC muonTag events = " << nEv_muonTag[0] << std::endl;
  std::cout << " tot JPsi_MC endSelection Events = " << nEv_cuts[0] << std::endl;
  

  for(int ij=0; ij<6; ++ij){
    h_Bmass_llbin[ij+1] = new TH1F(Form("h_Bmass_bin_%.2f-%.2f", llMassBoundary.at(ij), llMassBoundary.at(ij+1)), "", 1000, 0., 10.);

    nEv_muonTag[ij+1] = t1->Draw("Muon_sel_index", (cut_muonTag+" && "+llGenMassCut.at(ij)).c_str(), "goff");
    
    nEv_cuts[ij+1] = t1->Draw(Form("BToKstll_B_mass[BToKstll_sel_index] >> %s", h_Bmass_llbin[ij+1]->GetName()), (cut_muonTag+" && "+cut_Recocandidate+" && "+cut_chargeEff+" && "+cut_pTEff+" && "+cut_alphaEff+" && "+cut_vtxCLEff+" && "+cut_LxyEff+" && "+llMassCut.at(ij)+" && "+llGenMassCut.at(ij)).c_str());

    std::cout << " tot nonResonant_MC muonTag events = " << nEv_muonTag[ij+1] << " in mass bin " << llMassCut.at(ij) << std::endl;
    std::cout << " tot nonResonant_MC endSelection Events = " << nEv_cuts[ij+1] << std::endl;
  }//loop


  std::string outName = "outMassHistos_MC_ee.root";
  if(!isEleFinalState) outName ="outMassHistos_MC_mumu.root";
  TFile outMassHistos(outName.c_str(), "recreate");
  outMassHistos.cd();
  for(int ij=0; ij<7; ++ij)
    h_Bmass_llbin[ij]->Write(h_Bmass_llbin[ij]->GetName());
  outMassHistos.Close();



  //now fitting
  float nEv_postFit[7] = {0.};
  float nEvError_postFit[7] = {0.};

  RooWorkspace w("w");    
  w.factory("x[0, 10]");  
  //w.factory("x[4.5, 6.]");  

  w.factory("nbackground[10000, 0, 10000]");   
  w.factory("nsignal[100, 0.0, 10000.0]");

  for(int ij=0; ij<7; ++ij){

    h_Bmass_llbin[ij]->Rebin(2);

    w.factory("Gaussian::smodel(x,mu[5.3,4.5,6],sigma[0.05,0,0.2])");
    //w.factory("Gaussian::smodel(x,mu[5.3,4.5,6],sigma[0.05,0,0.2])");
    //w.factory("RooCBShape::smodel(x,m[5.3,4.5,6],s[0.1,0.,1.],a[1.2,0.,3.],n[1,0.1,6.])");
    //w.factory("RooCBShape::CBall(x[0,15], mean[11000,13000], sigma[5000,200000], alpha[0,10000],n[0,100000])");
    RooAbsPdf * smodel = w.pdf("smodel");
    
    w.factory("Exponential::bmodel(x,tau[-2,-3,0])");
    RooAbsPdf * bmodel = w.pdf("bmodel");
    
    w.factory("SUM::model(nbackground*bmodel, nsignal*smodel)");
    RooAbsPdf * model = w.pdf("model");
    
    RooDataHist hBMass("hBMass", "hBMass", *w.var("x"), Import(*(h_Bmass_llbin[ij])));

    RooFitResult * r = model->fitTo(hBMass, Minimizer("Minuit2"),Save(true));

    RooPlot * plot = w.var("x")->frame();
    if(isEleFinalState){
      if(ij == 0) plot->SetXTitle("K(JPsi)ee mass (GeV)");
      else plot->SetXTitle("Kee mass (GeV)");
    }
    else{
      if(ij == 0) plot->SetXTitle("K(JPsi)#mu#mu mass (GeV)");
      else plot->SetXTitle("K#mu#mu mass (GeV)");
    }
    plot->SetTitle("");
    plot->SetAxisRange(4.5,6);
    hBMass.plotOn(plot);
    model->plotOn(plot);
    model->plotOn(plot, Components("bmodel"),LineStyle(kDashed));
    model->plotOn(plot, Components("smodel"),LineColor(kRed));

    TCanvas * cc = new TCanvas();
    cc->SetLogy(0);
    plot->Draw();
    if(isEleFinalState) cc->Print(Form("plots_computeAE_ee_PR38/computeAE_ee_PR38_%s.png",h_Bmass_llbin[ij]->GetName()), "png");
    else cc->Print(Form("plots_computeAE_mumu_PR38/computeAE_mumu_PR38_%s.png",h_Bmass_llbin[ij]->GetName()), "png");

    RooRealVar* parS = (RooRealVar*) r->floatParsFinal().find("nsignal");
    RooRealVar* parB = (RooRealVar*) r->floatParsFinal().find("nbackground");
    nEv_postFit[ij] = parS->getValV();
    nEvError_postFit[ij] = parS->getError();

    std::cout << " JPsi selection signal events = \t " << parS->getValV() << " error = " << parS->getError() 
	      << " bkg events = " << parB->getValV() << " error = " << parB->getError() << std::endl;
  }
  
  float errb = pow(nEvError_postFit[0]/nEv_postFit[0], 2);

  std::cout << " ***** summary ***** "<< std::endl;
  for(int ij=0; ij<7; ++ij){

    float errRatio = pow(nEvError_postFit[ij]/nEv_postFit[ij], 2) + errb;
    errRatio = sqrt(errRatio) * nEvError_postFit[ij]/nEvError_postFit[0];

    std::cout << "\n \n  category: " << h_Bmass_llbin[ij]->GetName()
              << " \t integral muonTag = " << nEv_muonTag[ij]
              << " evtCount " << nEv_cuts[ij] << " withFit " << nEv_postFit[ij] <<"+/-"<<nEvError_postFit[ij]
              << "\n \t\t\t eff(/muonTag)    evtCount " << nEv_cuts[ij]/nEv_muonTag[ij]
              << " withFit " << nEv_postFit[ij]/nEv_muonTag[ij] << "+/-" << nEvError_postFit[ij]/nEv_muonTag[ij]
              << "\n  \t\t\t double ratio wrt JPsi    evtCount " << (nEv_cuts[ij]/nEv_muonTag[ij])/(nEv_cuts[0]/nEv_muonTag[0])
              << " withFit " << (nEv_postFit[ij]/nEv_muonTag[ij])/(nEv_postFit[0]/nEv_muonTag[0]) << "+/-"<<errRatio<< std::endl;
  }

} 
