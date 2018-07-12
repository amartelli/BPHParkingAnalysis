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
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TPad.h"
#include "TMultiGraph.h"


void analyzeCharged(){

  gROOT->Reset();
  gROOT->Macro("~/setStyle.C");
  gROOT->Macro("~/setStyle.C");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TChain* t1 = new TChain("Events");
  TChain* t2 = new TChain("Events");

  t1->Add("/vols/cms/amartell/BParking/ntuPROD/ntu_BToKee_v18_03_22_and_21.root");
  t2->Add("/vols/cms/amartell/BParking/ntuPROD/ntu_BToKJPsiee_v18_06_4.root");


  int nNNR = t1->GetEntries();
  int nReso = t2->GetEntries();

  std::cout << " #initial n. events: nnreso = " << nNNR << " resonant = " << nReso << std::endl;


  TH1F* h1 = new TH1F("h1", "", 500, 0, 1.1);
  TH1F* h2 = new TH1F("h2", "", 500, 0, 1.1);

  h1->Sumw2();
  h2->Sumw2();

  h1->SetLineColor(kRed);
  h2->SetLineColor(kBlue);
  h1->SetLineWidth(2);
  h2->SetLineWidth(2);
  
  std:: string cut_gen_charge = "BToKee_gen_index != -1 && BToKee_ele1_charge[BToKee_gen_index]*BToKee_ele2_charge[BToKee_gen_index] < 0.";

  t1->Draw("BToKee_cosAlpha[BToKee_gen_index] >> h1", cut_gen_charge.c_str());
  t2->Draw("BToKee_cosAlpha[BToKee_gen_index] >> h2", cut_gen_charge.c_str());

  int genSelectedOpposietCharge[2];
  genSelectedOpposietCharge[0] = h1->GetEntries();
  genSelectedOpposietCharge[1] = h2->GetEntries();

  std::cout << cut_gen_charge << " #survived: nnreso = " << genSelectedOpposietCharge[0] << " resonant = " << genSelectedOpposietCharge[1] << std::endl;

  h1->Scale(1./nNNR);
  h2->Scale(1./nReso);

  TLegend *legM = new TLegend(0.70,0.70,0.98,0.95,NULL,"brNDC");
  legM->SetTextFont(42);
  legM->SetFillColor(kWhite);
  legM->SetLineColor(kWhite);
  legM->SetShadowColor(kWhite);
  legM->SetFillStyle(0);
  legM->SetTextSize(0.05);
  legM->AddEntry(h1, "BToKee", "l");
  legM->AddEntry(h2, "BToJPsiKee", "l");
  

  TCanvas* c1 = new TCanvas();
  c1->cd();
  gPad->SetLogy();
  h1->GetXaxis()->SetTitle("cosAlpha");
  h1->Draw("hist");
  h2->Draw("hist, same");
  legM->Draw("same");
  c1->Print("plots/Kee_cosAlpha.png", "png");
  c1->Print("plots/Kee_cosAlpha.pdf", "pdf");
  c1->Print("plots/Kee_cosAlpha.root", "root");


  ////
  std::string cut_alpha = " && BToKee_cosAlpha[BToKee_gen_index] > 0.99";

  TH1F* h1b = new TH1F("h1b", "", 100, 0., 1.);
  TH1F* h2b = new TH1F("h2b", "", 100, 0., 1.);
  h1b->Sumw2();
  h2b->Sumw2();

  h1b->SetLineColor(kRed);
  h2b->SetLineColor(kBlue);
  h1b->SetLineWidth(2);
  h2b->SetLineWidth(2);


  int cosAlpha[2];
  cosAlpha[0] = t1->Draw("BToKee_CL_vtx[BToKee_gen_index] >> h1b", (cut_gen_charge+cut_alpha).c_str());
  cosAlpha[1] = t2->Draw("BToKee_CL_vtx[BToKee_gen_index] >> h2b", (cut_gen_charge+cut_alpha).c_str());


  std::cout << cut_alpha << " #survived: nnreso = " << cosAlpha[0] << " resonant = " << cosAlpha[1] << std::endl;

  h1b->Scale(1./nNNR);
  h2b->Scale(1./nReso);

  TCanvas* c1b = new TCanvas();
  c1b->cd();
  gPad->SetLogy();
  h1b->GetXaxis()->SetTitle("CL_vtx");
  h1b->Draw("hist");
  h2b->Draw("hist, same");
  legM->Draw("same");
  c1b->Print("plots/Kee_CL_vtx.png", "png");
  c1b->Print("plots/Kee_CL_vtx.pdf", "pdf");
  c1b->Print("plots/Kee_CL_vtx.root", "root");

  //  return;
  /////////////////

  std::string cut_CL = " && BToKee_CL_vtx[BToKee_gen_index] > 0.1";

  TH1F* h1c = new TH1F("h1c", "", 100, -50., 50.);
  TH1F* h2c = new TH1F("h2c", "", 100, -50., 50.);
  h1c->Sumw2();
  h2c->Sumw2();

  h1c->SetLineColor(kRed);
  h2c->SetLineColor(kBlue);
  h1c->SetLineWidth(2);
  h2c->SetLineWidth(2);

  int vtxCL[2];
  vtxCL[0] = t1->Draw("BToKee_kaon_DCASig[BToKee_gen_index] >> h1c", (cut_gen_charge+cut_alpha+cut_CL).c_str());
  vtxCL[1] = t2->Draw("BToKee_kaon_DCASig[BToKee_gen_index] >> h2c", (cut_gen_charge+cut_alpha+cut_CL).c_str());

  std::cout << cut_CL << " #survived: nnreso = " << vtxCL[0] << " resonant = " << vtxCL[1] << std::endl;


  h1c->Scale(1./nNNR);
  h2c->Scale(1./nReso);

  TCanvas* c1c = new TCanvas();

  h1c->GetXaxis()->SetTitle("kaon DCA SIP");
  h1c->Draw("hist");
  h2c->Draw("hist, same");
  legM->Draw("same");
  c1c->Print("plots/Kee_kaonDCA.png", "png");
  c1c->Print("plots/Kee_kaonDCA.pdf", "pdf");
  c1c->Print("plots/Kee_kaonDCA.root", "root");

  
  //////
  /////////////////
  
  std::string cut_DCA  = " && (BToKee_kaon_DCASig[BToKee_gen_index] > 6 || BToKee_kaon_DCASig[BToKee_gen_index] < -6)";

  TH1F* h1d = new TH1F("h1d", "", 100, 0., 100.);
  TH1F* h2d = new TH1F("h2d", "", 100, 0., 100.);
  h1d->Sumw2();
  h2d->Sumw2();

  h1d->SetLineColor(kRed);
  h2d->SetLineColor(kBlue);
  h1d->SetLineWidth(2);
  h2d->SetLineWidth(2);


  int kaonDCA[2];
  kaonDCA[0] = t1->Draw("BToKee_Lxy[BToKee_gen_index] >> h1d", (cut_gen_charge+cut_alpha+cut_CL+cut_DCA).c_str());
  kaonDCA[1] = t2->Draw("BToKee_Lxy[BToKee_gen_index] >> h2d", (cut_gen_charge+cut_alpha+cut_CL+cut_DCA).c_str());

  std::cout << cut_DCA << " #survived: nnreso = " << kaonDCA[0] << " resonant = " << kaonDCA[1] << std::endl;


  h1d->Scale(1./nNNR);
  h2d->Scale(1./nReso);

  TCanvas* c1d = new TCanvas();
  gPad->SetLogy();
  h1d->GetXaxis()->SetTitle("Lxy IP(cm)");
  h1d->Draw("hist");
  h2d->Draw("hist, same");
  legM->Draw("same");
  c1d->Print("plots/Kee_Lxy.png", "png");
  c1d->Print("plots/Kee_Lxy.pdf", "pdf");
  c1d->Print("plots/Kee_Lxy.root", "root");


  /////////////////
  
  //  std::string cut_Lxy  = " ";

  TH1F* h1e = new TH1F("h1e", "", 100, 0., 15.);
  TH1F* h2e = new TH1F("h2e", "", 100, 0., 15.);
  h1e->Sumw2();
  h2e->Sumw2();

  h1e->SetLineColor(kRed);
  h2e->SetLineColor(kBlue);
  h1e->SetLineWidth(2);
  h2e->SetLineWidth(2);


  t1->Draw("BToKee_ee_mass[BToKee_gen_index] >> h1e", (cut_gen_charge+cut_alpha+cut_CL+cut_DCA).c_str());
  t2->Draw("BToKee_ee_mass[BToKee_gen_index] >> h2e", (cut_gen_charge+cut_alpha+cut_CL+cut_DCA).c_str());

  std::cout << cut_DCA << " #survived: nnreso = " << kaonDCA[0] << " resonant = " << kaonDCA[1] << std::endl;


  h1e->Scale(1./nNNR);
  h2e->Scale(1./nReso);

  TCanvas* c1e = new TCanvas();
  gPad->SetLogy();
  h1e->GetXaxis()->SetTitle("Mee mass");
  h1e->Draw("hist");
  h2e->Draw("hist, same");
  legM->Draw("same");
  c1e->Print("plots/Kee_Mee.png", "png");
  c1e->Print("plots/Kee_Mee.pdf", "pdf");
  c1e->Print("plots/Kee_Mee.root", "root");


  //////

  TH2F* h1f = new TH2F("h1f", "", 500, 0., 15., 500, 0., 15.);
  TH2F* h2f = new TH2F("h2f", "", 500, 0., 15., 500, 0., 15.);
  h1f->Sumw2();
  h2f->Sumw2();

  h1f->SetMarkerColor(kRed);
  h2f->SetMarkerColor(kBlue);
  h1f->SetMarkerStyle(20);
  h2f->SetMarkerStyle(25);


  t1->Draw("BToKee_ee_mass[BToKee_gen_index]:BToKee_mass[BToKee_gen_index] >> h1f", (cut_gen_charge+cut_alpha+cut_CL+cut_DCA).c_str());
  t2->Draw("BToKee_ee_mass[BToKee_gen_index]:BToKee_mass[BToKee_gen_index] >> h2f", (cut_gen_charge+cut_alpha+cut_CL+cut_DCA).c_str());

  std::cout << cut_DCA << " #survived: nnreso = " << kaonDCA[0] << " resonant = " << kaonDCA[1] << std::endl;


  // h1e->Scale(1./nNNR);
  // h2e->Scale(1./nReso);

  TCanvas* c1f = new TCanvas();
  h1f->GetXaxis()->SetTitle("Kee mass");
  h1f->GetYaxis()->SetTitle("ee mass");
  h1f->Draw("hist");
  h2f->Draw("hist, same");
  legM->Draw("same");
  c1f->Print("plots/Kee_Mee_vsBmass.png", "png");
  c1f->Print("plots/Kee_Mee_vsBmass.pdf", "pdf");
  c1f->Print("plots/Kee_Mee_vsBmass.root", "root");


  ///////////

  TH1F* h1g = new TH1F("h1g", "", 100, 0., 15.);
  TH1F* h2g = new TH1F("h2g", "", 100, 0., 15.);
  h1g->Sumw2();
  h2g->Sumw2();

  h1g->SetLineColor(kRed);
  h2g->SetLineColor(kBlue);
  h1g->SetLineWidth(2);
  h2g->SetLineWidth(2);


  t1->Draw("BToKee_mass[BToKee_gen_index] >> h1g", (cut_gen_charge+cut_alpha+cut_CL+cut_DCA).c_str());
  t2->Draw("BToKee_mass[BToKee_gen_index] >> h2g", (cut_gen_charge+cut_alpha+cut_CL+cut_DCA).c_str());

  std::cout << cut_DCA << " #survived: nnreso = " << kaonDCA[0] << " resonant = " << kaonDCA[1] << std::endl;


  h1g->Scale(1./nNNR);
  h2g->Scale(1./nReso);

  TCanvas* c1g = new TCanvas();
  gPad->SetLogy();
  h1g->GetXaxis()->SetTitle("Kee mass");
  h1g->Draw("hist");
  h2g->Draw("hist, same");
  legM->Draw("same");
  c1g->Print("plots/Kee_mass.png", "png");
  c1g->Print("plots/Kee_mass.pdf", "pdf");
  c1g->Print("plots/Kee_mass.root", "root");

  /////////////////compute acceptance in q bins


  //recompute ee mass


  TH1F* JPsiMass_bin = new TH1F("JPsiMass_bin", "", 1000, 0., 10.);
  TH1F* hMass_bins[5];
  TH1F* JPsiMass_bin_cut0 = new TH1F("JPsiMass_bin_cut0", "", 1000, 0., 10.);
  TH1F* hMass_bins_cut0[5];
  std::vector<std::string> eeMassRange;
  std::vector<float> eeMassBoundary;
  eeMassBoundary.push_back(0.);
  eeMassBoundary.push_back(1.);
  eeMassBoundary.push_back(2.5);
  eeMassBoundary.push_back(2.9);
  eeMassBoundary.push_back(3.3);
  eeMassBoundary.push_back(3.58);

  float eventCounts[5];
  float eventCounts_cut0[5];

  for(int ij=0; ij<5; ++ij){
    hMass_bins[ij] = new TH1F(Form("hMass_eeRange_%.2f-%.2f",eeMassBoundary.at(ij), eeMassBoundary.at(ij+1)), "", 1000, 0., 10.);
    hMass_bins_cut0[ij] = new TH1F(Form("hMass_eeRange_cut0_%.2f-%.2f",eeMassBoundary.at(ij), eeMassBoundary.at(ij+1)), "", 1000, 0., 10.);
    std::string cut = Form("&& BToKee_eeKFit_ee_mass[BToKee_gen_index] > %.2f && BToKee_eeKFit_ee_mass[BToKee_gen_index] < %.2f", 
			   eeMassBoundary.at(ij), eeMassBoundary.at(ij+1));
    eeMassRange.push_back(cut);

    eventCounts[ij] = t1->Draw(Form("BToKee_eeKFit_ee_mass[BToKee_gen_index] >> %s", hMass_bins[ij]->GetName()), 
			       (cut_gen_charge+cut_alpha+cut_CL+cut_DCA+eeMassRange.at(ij)).c_str());

    eventCounts_cut0[ij] = t1->Draw(Form("BToKee_eeKFit_ee_mass[BToKee_gen_index] >> %s", hMass_bins_cut0[ij]->GetName()), 
				    (cut_gen_charge+eeMassRange.at(ij)).c_str());

    std::cout << " >>> non-resonant ee mass bin " << eeMassBoundary.at(ij) << "-" << eeMassBoundary.at(ij+1) 
	      << " nEvents cut0 = " << eventCounts_cut0[ij] << " post selection = " << eventCounts[ij] << " in total MC events = " << nNNR << std::endl;
  }


  float resonantEvents = t1->Draw("BToKee_eeKFit_ee_mass[BToKee_gen_index] >> JPsiMass_bin",
				  (cut_gen_charge+cut_alpha+cut_CL+cut_DCA+eeMassRange.at(3)).c_str());
  float resonantEvents_cut0 = t1->Draw("BToKee_eeKFit_ee_mass[BToKee_gen_index] >> JPsiMass_bin_cut0",
				       (cut_gen_charge+eeMassRange.at(3)).c_str());
  
  std::cout << " >>> resonant JPsi mass bin = " << eeMassRange.at(3) 
	    << " nEvents cut0 = "  << resonantEvents_cut0 << " post selection = " << resonantEvents << " in total MC events = " << nReso << std::endl;

  
  TFile outMassHistos("outMassHistos.root", "recreate");
  outMassHistos.cd();
  for(int ij=0; ij<5; ++ij){
    hMass_bins[ij]->Write(Form("hMass_eeRange_%.2f-%.2f",eeMassBoundary.at(ij), eeMassBoundary.at(ij+1)));
    hMass_bins_cut0[ij]->Write(Form("hMass_eeRange_cut0_%.2f-%.2f",eeMassBoundary.at(ij), eeMassBoundary.at(ij+1)));
  }
  JPsiMass_bin->Write();
  JPsiMass_bin_cut0->Write();
  outMassHistos.Close();


}
