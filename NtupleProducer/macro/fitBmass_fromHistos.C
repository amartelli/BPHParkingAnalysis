//plot and fit BToKll 

//example to run
//root fitBmass_fromHistos.C'(1, "/vols/cms/amartell/BParking/data_allStat/outMassHistos_Kee_tightSel_KeePrefitMass.root")'
//root fitBmass_fromHistos.C'(0, "/vols/cms/amartell/BParking/data_allStat/outMassHistos_Kmumu_tightSel.root")'


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
#include "TText.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"

using namespace RooFit;


void fitBmass_fromHistos(int isEleFinalState, std::string inFile){

  gROOT->Reset();
  gROOT->Macro("setStyle.C");

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  std::cout << " isEleFinalState = " << isEleFinalState << std::endl;


  TFile* inF = TFile::Open(inFile.c_str());

  std::vector<float> llMassBoundary;
  llMassBoundary.push_back(0.);
  llMassBoundary.push_back(1.);
  llMassBoundary.push_back(2.5);
  llMassBoundary.push_back(2.9);
  llMassBoundary.push_back(3.2);
  llMassBoundary.push_back(3.58);
  llMassBoundary.push_back(100.);


  TH1F* h_Bmass[7];
  TH1F* h_Bmass_llt[7];
  TH1F* h_Bmass_ltt[7];
  for(int ij=0; ij<7; ++ij){
    h_Bmass[ij] = (TH1F*)inF->Get(Form("Bmass_%d", ij))->Clone(Form("h_Bmass_%d", ij));
    h_Bmass_llt[ij] = (TH1F*)inF->Get(Form("Bmass_llt_%d", ij))->Clone(Form("h_Bmass_llt_%d", ij));
    h_Bmass_ltt[ij] = (TH1F*)inF->Get(Form("Bmass_ltt_%d", ij))->Clone(Form("h_Bmass_ltt_%d", ij));
    //    h_Bmass[ij]->GetXaxis()->SetRangeUser(4.5, 6.);
  }//loop


  //now fitting
  float nEv_postFit[7] = {0.};
  float nEvError_postFit[7] = {0.};
  float nBkg_postFit[7] = {0.};
  float nBkgError_postFit[7] = {0.};
  float nBkgInt_postFit[7] = {0.};
  float nBkgIntError_postFit[7] = {0.};
  float chi2[7] = {0.};

  for(int ij=0; ij<7; ++ij){
    RooWorkspace w("w");    
    w.factory("x[4.5, 6.]");  

    w.factory("nbackground[10000, 0, 100000]");   
    w.factory("nbackgroundR[10, 0, 100]");   
    w.factory("nsignal[100, 0.0, 10000]");
    
    w.factory("Gaussian::smodel(x,mu[5.3,4.5,6],sigma[0.05,0,0.15])");
    //w.factory("Gaussian::smodel(x,mu[5.3,4.5,6],sigma[0.02,0,0.03])");
    //w.factory("RooCBShape::smodel(x,m[5.3,4.5,6],s[0.1,0.,1.],a[1.2,0.,3.],n[1,0.1,6.])");
    //w.factory("RooCBShape::CBall(x[0,15], mean[11000,13000], sigma[5000,200000], alpha[0,10000],n[0,100000])");
    RooAbsPdf * smodel = w.pdf("smodel");
    
    w.factory("Exponential::bmodel1(x,tau[-2,-3,0])");
    RooAbsPdf * bmodel1 = w.pdf("bmodel1");

    w.factory("Gaussian::bmodel2(x,mub[5.3,4.5,6],sigmab[0.15,0.05,2.])");
    RooAbsPdf * bmodel2 = w.pdf("bmodel2");
    
    w.factory("SUM::modelb(nbackground * bmodel1, -nbackgroundR * bmodel1, nbackgroundR * bmodel2)");
    RooAbsPdf * modelb = w.pdf("modelb");


    w.factory("SUM::model(nbackground * modelb, nsignal * smodel)");
    RooAbsPdf * model = w.pdf("model");
    
    RooDataHist hBMass("hBMass", "hBMass", *w.var("x"), Import(*(h_Bmass[ij])));
    RooDataHist hBMass_llt("hBMass_llt", "hBMass_llt", *w.var("x"), Import(*(h_Bmass_llt[ij])));
    RooDataHist hBMass_ltt("hBMass_ltt", "hBMass_ltt", *w.var("x"), Import(*(h_Bmass_ltt[ij])));
    w.Print();

    RooFitResult * r = model->fitTo(hBMass, Minimizer("Minuit2"),Save(true));
    std::cout << " fit status = " << r->status() << std::endl;

    RooPlot * plot = w.var("x")->frame();
    if(isEleFinalState){
      if(ij == 3) plot->SetXTitle("K(JPsi)ee mass (GeV)");
      else plot->SetXTitle("Kee mass (GeV)");
    }
    else{
      if(ij == 3) plot->SetXTitle("K(JPsi)#mu#mu mass (GeV)");
      else plot->SetXTitle("K#mu#mu mass (GeV)");
    }
    plot->SetTitle("");
    plot->SetAxisRange(4.,6);
    hBMass.plotOn(plot);
    model->plotOn(plot);
    model->plotOn(plot, Components("modelb"),LineStyle(kDashed));
    model->plotOn(plot, Components("smodel"),LineColor(kRed));
    hBMass_llt.plotOn(plot,LineColor(kGreen+2),MarkerColor(kGreen+2)) ;
    hBMass_ltt.plotOn(plot,LineColor(kViolet),MarkerColor(kViolet)) ;
    chi2[ij] = plot->chiSquare();

    RooRealVar* parS = (RooRealVar*) r->floatParsFinal().find("nsignal");
    RooRealVar* parB = (RooRealVar*) r->floatParsFinal().find("nbackground");
    nEv_postFit[ij] = parS->getValV();
    nEvError_postFit[ij] = parS->getError();
    nBkg_postFit[ij] = parB->getValV();
    nBkgError_postFit[ij] = parB->getError();

    std::cout << " **** JPsi selection signal events = \t " << parS->getValV() << " error = " << parS->getError() 
	      << " bkg events = " << parB->getValV() << " error = " << parB->getError() << std::endl;


    RooRealVar* parMean = (RooRealVar*) r->floatParsFinal().find("mu");
    RooRealVar* parSigma = (RooRealVar*) r->floatParsFinal().find("sigma");

    float meanVal = parMean->getValV();
    float sigmaVal = parSigma->getValV();

    std::cout << "\n  parMean = " << parMean->getValV() << " parSigma = " << parSigma->getValV() << std::endl;
    
    w.var("x")->setRange("signalRange", meanVal - 3.*sigmaVal, meanVal + 3.*sigmaVal);
    //RooAbsReal* bkgIntegral = w.pdf("bmodel")->createIntegral(*(w.var("x")), NormSet(*(w.var("x"))), Range("signalRange")) ;
    RooAbsReal* bkgIntegral = w.pdf("modelb")->createIntegral(*(w.var("x")), NormSet(*(w.var("x"))), Range("signalRange")) ;
    // w.var("x")->setRange("signalRange", 4.5, 6.);
    // RooAbsReal* bkgIntegral = w.pdf("bmodel")->createIntegral(*(w.var("x")), NormSet(*(w.var("x"))), Range("signalRange")) ;
    std::cout << "\n  bkgIntegral = " << bkgIntegral->getVal() << std::endl;
    nBkgInt_postFit[ij] = bkgIntegral->getVal() * nBkg_postFit[ij];
    nBkgIntError_postFit[ij] = nBkgError_postFit[ij] * bkgIntegral->getVal();

    // Add text to frame
    /*
    TText* txt = new TText(0.8,0.9, Form("Signal %.1f +/- %1.f",nEv_postFit[ij], nEvError_postFit[ij]));
    TText* txt2 = new TText(0.8,0.8, Form("Bkg %.1f +/- %1.f", nBkgInt_postFit[ij], nBkgIntError_postFit[ij]));
    txt->SetTextSize(0.04);
    txt->SetTextColor(kRed);
    txt2->SetTextSize(0.04);
    txt2->SetTextColor(kRed);
    plot->addObject(txt);
    plot->addObject(txt2);
    */

    TCanvas * cc = new TCanvas();
    cc->SetLogy(0);
    plot->Draw();

    TLatex tL;
    tL.SetNDC();
    tL.SetTextSize(0.05);
    tL.SetTextFont(42);
    tL.DrawLatex(0.65,0.9, Form("S %.1f +/- %1.f",nEv_postFit[ij], nEvError_postFit[ij]));
    TLatex tL2;
    tL2.SetNDC();
    tL2.SetTextSize(0.05);
    tL2.SetTextFont(42);
    tL2.DrawLatex(0.65,0.85, Form("B %.1f +/- %1.f",nBkgInt_postFit[ij], nBkgIntError_postFit[ij]));


    if(isEleFinalState) cc->Print(Form("plots/Bmass_DATA/Kee_%s.png",h_Bmass[ij]->GetName()), "png");
    else cc->Print(Form("plots/Bmass_DATA/Kmumu_%s.png",h_Bmass[ij]->GetName()), "png");
  }
  

  std::cout << " ***** summary ***** "<< std::endl;
  for(int ij=0; ij<7; ++ij){

    std::cout << "\n \n  category: " << h_Bmass[ij]->GetName()
	      << "\n \t signal = " << nEv_postFit[ij] << "+/-" << nEvError_postFit[ij]
              << "\n \t bkg = " << nBkg_postFit[ij] << "+/-" << nBkgError_postFit[ij] 
	      << "\n \t bkg in +/- 3sigma = " << nBkgInt_postFit[ij] << " +/- " << nBkgIntError_postFit[ij]
	      << " chi2 = " << chi2[ij] << std::endl;
  }

}

