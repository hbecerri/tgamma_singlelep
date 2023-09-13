//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Sep 13 10:37:40 2023 by ROOT version 6.24/08
// from TTree Events/
// found on file: ../output/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/0000.root
//////////////////////////////////////////////////////////

#ifndef Events_h
#define Events_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class Events {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           nLepton_pt_eta_phi_mass;
   Float_t         Lepton_pt_eta_phi_mass_pt[2];   //[nLepton_pt,eta,phi,mass]
   Float_t         Lepton_pt_eta_phi_mass_eta[2];   //[nLepton_pt,eta,phi,mass]
   Float_t         Lepton_pt_eta_phi_mass_phi[2];   //[nLepton_pt,eta,phi,mass]
   Float_t         Lepton_pt_eta_phi_mass_mass[2];   //[nLepton_pt,eta,phi,mass]
   Int_t           nLepton_charge;
   Int_t           Lepton_charge[2];   //[nLepton_charge]
   Bool_t          masktophad_m;
   Int_t           ntophad_m;
   Double_t        tophad_m[1];   //[ntophad_m]
   Bool_t          masktophad_pt;
   Int_t           ntophad_pt;
   Double_t        tophad_pt[1];   //[ntophad_pt]
   Bool_t          masktophad_eta;
   Int_t           ntophad_eta;
   Double_t        tophad_eta[1];   //[ntophad_eta]
   Bool_t          masktophad_phi;
   Int_t           ntophad_phi;
   Double_t        tophad_phi[1];   //[ntophad_phi]
   Bool_t          masktoplep_m;
   Int_t           ntoplep_m;
   Double_t        toplep_m[1];   //[ntoplep_m]
   Bool_t          masktoplep_pt;
   Int_t           ntoplep_pt;
   Double_t        toplep_pt[1];   //[ntoplep_pt]
   Bool_t          masktoplep_eta;
   Int_t           ntoplep_eta;
   Double_t        toplep_eta[1];   //[ntoplep_eta]
   Bool_t          masktoplep_phi;
   Int_t           ntoplep_phi;
   Double_t        toplep_phi[1];   //[ntoplep_phi]
   Double_t        chi2;
   Int_t           nJet_pt_eta_phi_mass;
   Float_t         Jet_pt_eta_phi_mass_pt[11];   //[nJet_pt,eta,phi,mass]
   Float_t         Jet_pt_eta_phi_mass_eta[11];   //[nJet_pt,eta,phi,mass]
   Float_t         Jet_pt_eta_phi_mass_phi[11];   //[nJet_pt,eta,phi,mass]
   Float_t         Jet_pt_eta_phi_mass_mass[11];   //[nJet_pt,eta,phi,mass]
   Int_t           nJet_btagged;
   Bool_t          Jet_btagged[11];   //[nJet_btagged]
   Float_t         MET_pt_phi_pt;
   Float_t         MET_pt_phi_phi;
   Int_t           nGenTopPos_pt_phi_eta_mass;
   Float_t         GenTopPos_pt_phi_eta_mass_pt[1];   //[nGenTopPos_pt,phi,eta,mass]
   Float_t         GenTopPos_pt_phi_eta_mass_phi[1];   //[nGenTopPos_pt,phi,eta,mass]
   Float_t         GenTopPos_pt_phi_eta_mass_eta[1];   //[nGenTopPos_pt,phi,eta,mass]
   Float_t         GenTopPos_pt_phi_eta_mass_mass[1];   //[nGenTopPos_pt,phi,eta,mass]
   Int_t           nGenTopNeg_pt_phi_eta_mass;
   Float_t         GenTopNeg_pt_phi_eta_mass_pt[1];   //[nGenTopNeg_pt,phi,eta,mass]
   Float_t         GenTopNeg_pt_phi_eta_mass_phi[1];   //[nGenTopNeg_pt,phi,eta,mass]
   Float_t         GenTopNeg_pt_phi_eta_mass_eta[1];   //[nGenTopNeg_pt,phi,eta,mass]
   Float_t         GenTopNeg_pt_phi_eta_mass_mass[1];   //[nGenTopNeg_pt,phi,eta,mass]
   Int_t           nNeutrino1_pt_phi_eta_mass;
   Double_t        Neutrino1_pt_phi_eta_mass_pt[2];   //[nNeutrino1_pt,phi,eta,mass]
   Double_t        Neutrino1_pt_phi_eta_mass_phi[2];   //[nNeutrino1_pt,phi,eta,mass]
   Double_t        Neutrino1_pt_phi_eta_mass_eta[2];   //[nNeutrino1_pt,phi,eta,mass]
   Double_t        Neutrino1_pt_phi_eta_mass_mass[2];   //[nNeutrino1_pt,phi,eta,mass]
   Double_t        weight;
   Double_t        systematics_weight;

   // List of branches
   TBranch        *b_nLepton_pt_eta_phi_mass;   //!
   TBranch        *b_Lepton_pt_eta_phi_mass_pt;   //!
   TBranch        *b_Lepton_pt_eta_phi_mass_eta;   //!
   TBranch        *b_Lepton_pt_eta_phi_mass_phi;   //!
   TBranch        *b_Lepton_pt_eta_phi_mass_mass;   //!
   TBranch        *b_nLepton_charge;   //!
   TBranch        *b_Lepton_charge;   //!
   TBranch        *b_masktophad_m;   //!
   TBranch        *b_ntophad_m;   //!
   TBranch        *b_tophad_m;   //!
   TBranch        *b_masktophad_pt;   //!
   TBranch        *b_ntophad_pt;   //!
   TBranch        *b_tophad_pt;   //!
   TBranch        *b_masktophad_eta;   //!
   TBranch        *b_ntophad_eta;   //!
   TBranch        *b_tophad_eta;   //!
   TBranch        *b_masktophad_phi;   //!
   TBranch        *b_ntophad_phi;   //!
   TBranch        *b_tophad_phi;   //!
   TBranch        *b_masktoplep_m;   //!
   TBranch        *b_ntoplep_m;   //!
   TBranch        *b_toplep_m;   //!
   TBranch        *b_masktoplep_pt;   //!
   TBranch        *b_ntoplep_pt;   //!
   TBranch        *b_toplep_pt;   //!
   TBranch        *b_masktoplep_eta;   //!
   TBranch        *b_ntoplep_eta;   //!
   TBranch        *b_toplep_eta;   //!
   TBranch        *b_masktoplep_phi;   //!
   TBranch        *b_ntoplep_phi;   //!
   TBranch        *b_toplep_phi;   //!
   TBranch        *b_chi2;   //!
   TBranch        *b_nJet_pt_eta_phi_mass;   //!
   TBranch        *b_Jet_pt_eta_phi_mass_pt;   //!
   TBranch        *b_Jet_pt_eta_phi_mass_eta;   //!
   TBranch        *b_Jet_pt_eta_phi_mass_phi;   //!
   TBranch        *b_Jet_pt_eta_phi_mass_mass;   //!
   TBranch        *b_nJet_btagged;   //!
   TBranch        *b_Jet_btagged;   //!
   TBranch        *b_MET_pt_phi_pt;   //!
   TBranch        *b_MET_pt_phi_phi;   //!
   TBranch        *b_nGenTopPos_pt_phi_eta_mass;   //!
   TBranch        *b_GenTopPos_pt_phi_eta_mass_pt;   //!
   TBranch        *b_GenTopPos_pt_phi_eta_mass_phi;   //!
   TBranch        *b_GenTopPos_pt_phi_eta_mass_eta;   //!
   TBranch        *b_GenTopPos_pt_phi_eta_mass_mass;   //!
   TBranch        *b_nGenTopNeg_pt_phi_eta_mass;   //!
   TBranch        *b_GenTopNeg_pt_phi_eta_mass_pt;   //!
   TBranch        *b_GenTopNeg_pt_phi_eta_mass_phi;   //!
   TBranch        *b_GenTopNeg_pt_phi_eta_mass_eta;   //!
   TBranch        *b_GenTopNeg_pt_phi_eta_mass_mass;   //!
   TBranch        *b_nNeutrino1_pt_phi_eta_mass;   //!
   TBranch        *b_Neutrino1_pt_phi_eta_mass_pt;   //!
   TBranch        *b_Neutrino1_pt_phi_eta_mass_phi;   //!
   TBranch        *b_Neutrino1_pt_phi_eta_mass_eta;   //!
   TBranch        *b_Neutrino1_pt_phi_eta_mass_mass;   //!
   TBranch        *b_weight;   //!
   TBranch        *b_systematics_weight;   //!

   Events(TTree *tree=0);
   virtual ~Events();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Events_cxx
Events::Events(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../output/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/0000.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../output/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/0000.root");
      }
      f->GetObject("Events",tree);

   }
   Init(tree);
}

Events::~Events()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Events::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Events::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Events::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("nLepton_pt,eta,phi,mass", &nLepton_pt_eta_phi_mass, &b_nLepton_pt_eta_phi_mass);
   fChain->SetBranchAddress("Lepton_pt,eta,phi,mass_pt", Lepton_pt_eta_phi_mass_pt, &b_Lepton_pt_eta_phi_mass_pt);
   fChain->SetBranchAddress("Lepton_pt,eta,phi,mass_eta", Lepton_pt_eta_phi_mass_eta, &b_Lepton_pt_eta_phi_mass_eta);
   fChain->SetBranchAddress("Lepton_pt,eta,phi,mass_phi", Lepton_pt_eta_phi_mass_phi, &b_Lepton_pt_eta_phi_mass_phi);
   fChain->SetBranchAddress("Lepton_pt,eta,phi,mass_mass", Lepton_pt_eta_phi_mass_mass, &b_Lepton_pt_eta_phi_mass_mass);
   fChain->SetBranchAddress("nLepton_charge", &nLepton_charge, &b_nLepton_charge);
   fChain->SetBranchAddress("Lepton_charge", Lepton_charge, &b_Lepton_charge);
   fChain->SetBranchAddress("masktophad_m", &masktophad_m, &b_masktophad_m);
   fChain->SetBranchAddress("ntophad_m", &ntophad_m, &b_ntophad_m);
   fChain->SetBranchAddress("tophad_m", tophad_m, &b_tophad_m);
   fChain->SetBranchAddress("masktophad_pt", &masktophad_pt, &b_masktophad_pt);
   fChain->SetBranchAddress("ntophad_pt", &ntophad_pt, &b_ntophad_pt);
   fChain->SetBranchAddress("tophad_pt", tophad_pt, &b_tophad_pt);
   fChain->SetBranchAddress("masktophad_eta", &masktophad_eta, &b_masktophad_eta);
   fChain->SetBranchAddress("ntophad_eta", &ntophad_eta, &b_ntophad_eta);
   fChain->SetBranchAddress("tophad_eta", tophad_eta, &b_tophad_eta);
   fChain->SetBranchAddress("masktophad_phi", &masktophad_phi, &b_masktophad_phi);
   fChain->SetBranchAddress("ntophad_phi", &ntophad_phi, &b_ntophad_phi);
   fChain->SetBranchAddress("tophad_phi", tophad_phi, &b_tophad_phi);
   fChain->SetBranchAddress("masktoplep_m", &masktoplep_m, &b_masktoplep_m);
   fChain->SetBranchAddress("ntoplep_m", &ntoplep_m, &b_ntoplep_m);
   fChain->SetBranchAddress("toplep_m", toplep_m, &b_toplep_m);
   fChain->SetBranchAddress("masktoplep_pt", &masktoplep_pt, &b_masktoplep_pt);
   fChain->SetBranchAddress("ntoplep_pt", &ntoplep_pt, &b_ntoplep_pt);
   fChain->SetBranchAddress("toplep_pt", toplep_pt, &b_toplep_pt);
   fChain->SetBranchAddress("masktoplep_eta", &masktoplep_eta, &b_masktoplep_eta);
   fChain->SetBranchAddress("ntoplep_eta", &ntoplep_eta, &b_ntoplep_eta);
   fChain->SetBranchAddress("toplep_eta", toplep_eta, &b_toplep_eta);
   fChain->SetBranchAddress("masktoplep_phi", &masktoplep_phi, &b_masktoplep_phi);
   fChain->SetBranchAddress("ntoplep_phi", &ntoplep_phi, &b_ntoplep_phi);
   fChain->SetBranchAddress("toplep_phi", toplep_phi, &b_toplep_phi);
   fChain->SetBranchAddress("chi2", &chi2, &b_chi2);
   fChain->SetBranchAddress("nJet_pt,eta,phi,mass", &nJet_pt_eta_phi_mass, &b_nJet_pt_eta_phi_mass);
   fChain->SetBranchAddress("Jet_pt,eta,phi,mass_pt", Jet_pt_eta_phi_mass_pt, &b_Jet_pt_eta_phi_mass_pt);
   fChain->SetBranchAddress("Jet_pt,eta,phi,mass_eta", Jet_pt_eta_phi_mass_eta, &b_Jet_pt_eta_phi_mass_eta);
   fChain->SetBranchAddress("Jet_pt,eta,phi,mass_phi", Jet_pt_eta_phi_mass_phi, &b_Jet_pt_eta_phi_mass_phi);
   fChain->SetBranchAddress("Jet_pt,eta,phi,mass_mass", Jet_pt_eta_phi_mass_mass, &b_Jet_pt_eta_phi_mass_mass);
   fChain->SetBranchAddress("nJet_btagged", &nJet_btagged, &b_nJet_btagged);
   fChain->SetBranchAddress("Jet_btagged", Jet_btagged, &b_Jet_btagged);
   fChain->SetBranchAddress("MET_pt,phi_pt", &MET_pt_phi_pt, &b_MET_pt_phi_pt);
   fChain->SetBranchAddress("MET_pt,phi_phi", &MET_pt_phi_phi, &b_MET_pt_phi_phi);
   fChain->SetBranchAddress("nGenTopPos_pt,phi,eta,mass", &nGenTopPos_pt_phi_eta_mass, &b_nGenTopPos_pt_phi_eta_mass);
   fChain->SetBranchAddress("GenTopPos_pt,phi,eta,mass_pt", GenTopPos_pt_phi_eta_mass_pt, &b_GenTopPos_pt_phi_eta_mass_pt);
   fChain->SetBranchAddress("GenTopPos_pt,phi,eta,mass_phi", GenTopPos_pt_phi_eta_mass_phi, &b_GenTopPos_pt_phi_eta_mass_phi);
   fChain->SetBranchAddress("GenTopPos_pt,phi,eta,mass_eta", GenTopPos_pt_phi_eta_mass_eta, &b_GenTopPos_pt_phi_eta_mass_eta);
   fChain->SetBranchAddress("GenTopPos_pt,phi,eta,mass_mass", GenTopPos_pt_phi_eta_mass_mass, &b_GenTopPos_pt_phi_eta_mass_mass);
   fChain->SetBranchAddress("nGenTopNeg_pt,phi,eta,mass", &nGenTopNeg_pt_phi_eta_mass, &b_nGenTopNeg_pt_phi_eta_mass);
   fChain->SetBranchAddress("GenTopNeg_pt,phi,eta,mass_pt", GenTopNeg_pt_phi_eta_mass_pt, &b_GenTopNeg_pt_phi_eta_mass_pt);
   fChain->SetBranchAddress("GenTopNeg_pt,phi,eta,mass_phi", GenTopNeg_pt_phi_eta_mass_phi, &b_GenTopNeg_pt_phi_eta_mass_phi);
   fChain->SetBranchAddress("GenTopNeg_pt,phi,eta,mass_eta", GenTopNeg_pt_phi_eta_mass_eta, &b_GenTopNeg_pt_phi_eta_mass_eta);
   fChain->SetBranchAddress("GenTopNeg_pt,phi,eta,mass_mass", GenTopNeg_pt_phi_eta_mass_mass, &b_GenTopNeg_pt_phi_eta_mass_mass);
   fChain->SetBranchAddress("nNeutrino1_pt,phi,eta,mass", &nNeutrino1_pt_phi_eta_mass, &b_nNeutrino1_pt_phi_eta_mass);
   fChain->SetBranchAddress("Neutrino1_pt,phi,eta,mass_pt", Neutrino1_pt_phi_eta_mass_pt, &b_Neutrino1_pt_phi_eta_mass_pt);
   fChain->SetBranchAddress("Neutrino1_pt,phi,eta,mass_phi", Neutrino1_pt_phi_eta_mass_phi, &b_Neutrino1_pt_phi_eta_mass_phi);
   fChain->SetBranchAddress("Neutrino1_pt,phi,eta,mass_eta", Neutrino1_pt_phi_eta_mass_eta, &b_Neutrino1_pt_phi_eta_mass_eta);
   fChain->SetBranchAddress("Neutrino1_pt,phi,eta,mass_mass", Neutrino1_pt_phi_eta_mass_mass, &b_Neutrino1_pt_phi_eta_mass_mass);
   fChain->SetBranchAddress("weight", &weight, &b_weight);
   fChain->SetBranchAddress("systematics_weight", &systematics_weight, &b_systematics_weight);
   Notify();
}

Bool_t Events::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Events::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Events::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Events_cxx
