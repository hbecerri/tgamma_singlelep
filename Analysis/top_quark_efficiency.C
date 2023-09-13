#include <algorithm>
#include <iterator>
#include <TROOT.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLatex.h>
#include "TCanvas.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH1F.h"
#include "THStack.h"
#include "TRandom.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TFrame.h"
#include "TPaveLabel.h"
#include "TPad.h"
#include "TLegend.h"
#include "TRandom3.h"
#include <math.h> 

void top_quark_efficiency()
{

  TFile *_file0 = TFile::Open("../output/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/0000.root");

  Int_t Lepton_charge[1];
  Double_t tophad_phi[1];
  Double_t toplep_phi[1];
  Double_t tophad_eta[1];
  Double_t toplep_eta[1];   
  Double_t tophad_m[1];
  Double_t toplep_m[1];

  
  Float_t GenTopPos_pt_phi_eta_mass_phi;
  Float_t GenTopNeg_pt_phi_eta_mass_phi;

  Float_t GenTopPos_pt_phi_eta_mass_eta;
  Float_t GenTopNeg_pt_phi_eta_mass_eta;


  Double_t delta_top;
  Double_t delta_atop;  

  TTree *t1 = (TTree *) _file0->Get("Events"); 


  t1->SetBranchAddress("Lepton_charge",&Lepton_charge[0]);
  t1->SetBranchAddress("tophad_phi",&tophad_phi[0]);
  t1->SetBranchAddress("toplep_phi",&toplep_phi[0]);
  t1->SetBranchAddress("tophad_eta",&tophad_eta[0]);
  t1->SetBranchAddress("toplep_eta",&toplep_eta[0]);
  t1->SetBranchAddress("tophad_m",&tophad_m[0]);
  t1->SetBranchAddress("toplep_m",&toplep_m[0]);


  t1->SetBranchAddress("GenTopPos_pt,phi,eta,mass_eta",&GenTopPos_pt_phi_eta_mass_eta);
  t1->SetBranchAddress("GenTopNeg_pt,phi,eta,mass_eta",&GenTopNeg_pt_phi_eta_mass_eta);

  t1->SetBranchAddress("GenTopPos_pt,phi,eta,mass_phi",&GenTopPos_pt_phi_eta_mass_phi);
  t1->SetBranchAddress("GenTopNeg_pt,phi,eta,mass_phi",&GenTopNeg_pt_phi_eta_mass_phi);



  TH1F *Delta_top  = new TH1F("Delta_top","",50,0,10);  
  TH1F *Delta_atop = new TH1F("Delta_atop","",50,0,10); 
 
  for(int i=0; i<t1->GetEntries(); i++){

    t1->GetEntry(i);

    if(toplep_m[0] < 20) continue;
    if(tophad_m[0] < 20) continue;
  
    if(Lepton_charge[0]>0){   
    
      delta_top  = sqrt((GenTopPos_pt_phi_eta_mass_phi - toplep_phi[0])*(GenTopPos_pt_phi_eta_mass_phi - toplep_phi[0]) + (GenTopPos_pt_phi_eta_mass_eta - toplep_eta[0])*(GenTopPos_pt_phi_eta_mass_eta - toplep_eta[0]));
      delta_atop = sqrt((GenTopNeg_pt_phi_eta_mass_phi - tophad_phi[0])*(GenTopNeg_pt_phi_eta_mass_phi - tophad_phi[0]) + (GenTopNeg_pt_phi_eta_mass_eta - tophad_eta[0])*(GenTopNeg_pt_phi_eta_mass_eta - tophad_eta[0]));
 
    }else{

      delta_top  = sqrt((GenTopPos_pt_phi_eta_mass_phi - tophad_phi[0])*(GenTopPos_pt_phi_eta_mass_phi - tophad_phi[0]) + (GenTopPos_pt_phi_eta_mass_eta - tophad_eta[0])*(GenTopPos_pt_phi_eta_mass_eta - tophad_eta[0]));
      delta_atop = sqrt((GenTopNeg_pt_phi_eta_mass_phi - toplep_phi[0])*(GenTopNeg_pt_phi_eta_mass_phi - toplep_phi[0]) + (GenTopNeg_pt_phi_eta_mass_eta - toplep_eta[0])*(GenTopNeg_pt_phi_eta_mass_eta - toplep_eta[0]));
   
    }

    Delta_top->Fill(delta_top);
    Delta_atop->Fill(delta_atop);

  }


auto c1    = new TCanvas("c1","c1",600,400);
c1->cd();  
Delta_top->Draw();
c1->Print("Anglar_separation_top.pdf");


auto c2    = new TCanvas("c2","c2",600,400);
c2->cd();
Delta_atop->Draw();
c2->Print("Anglar_separation_atop.pdf");

 

}

