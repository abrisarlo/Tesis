/*
This script reads a post Delphes root file,
-lepton veto
-exactly 4 b-jets
-calculates thrust
-saves minimum info

run as:
root -l -q 'SkimAndCompute_bbbb.cxx+("../Samples/GammaGammaHHESpreadAllSiD2024XCC.root")'


*/
	
R__ADD_INCLUDE_PATH(/home/boc/delphes)
R__ADD_INCLUDE_PATH(/home/boc/delphes/external)
R__ADD_INCLUDE_PATH(/home/boc/delphes/external/ExRootAnalysis)
R__LOAD_LIBRARY(/home/boc/delphes/libDelphes.so)

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TSystem.h>
#include <TMath.h>

#include <vector>
#include <algorithm>

#include "classes/DelphesClasses.h"

float Magnitude(TLorentzVector v) { return v.P(); }

float findThrust(std::vector<TLorentzVector>& momenta,
                 TLorentzVector& thrustAxis)
{
    float maxThrust = 0.0;

    for (float theta = 0; theta < TMath::Pi(); theta += 0.1) {
        for (float phi = 0; phi < 2*TMath::Pi(); phi += 0.1) {

            TLorentzVector axis(
                TMath::Sin(theta)*TMath::Cos(phi),
                TMath::Sin(theta)*TMath::Sin(phi),
                TMath::Cos(theta),
                0.0
            );

            float sumP = 0.0;
            float sumParallel = 0.0;

            for (auto& p : momenta) {
                sumP += p.P();
                sumParallel += TMath::Abs(p.Dot(axis) / Magnitude(axis));
            }

            float thrust = sumParallel / sumP;
            if (thrust > maxThrust) {
                maxThrust = thrust;
                thrustAxis = axis;
            }
        }
    }
    return maxThrust;
}

void SkimAndCompute_bbbb(const char* inputFile)
{
    gROOT->SetBatch(kTRUE);
    gSystem->Load("libDelphes.so");
    TString inFile(inputFile);
    TString baseName = gSystem->BaseName(inFile);
    TString outDir   = "./SkimedSamples";

    gSystem->mkdir(outDir, kTRUE);

    TString outFile = outDir + "/" + baseName;

    TFile *fin = TFile::Open(inFile);
    if (!fin || fin->IsZombie()) {
        Error("SkimAndCompute_bbbb","Cannot open input file");
        return;
    }

    TFile *fout = new TFile(outFile, "RECREATE");

    TTree *tin = (TTree*)fin->Get("Delphes");

    TClonesArray *jets = nullptr;
    TClonesArray *electrons = nullptr;
    TClonesArray *muons = nullptr;
    TClonesArray *tracks = nullptr;
    TClonesArray *photons = nullptr;
    TClonesArray *neutrals = nullptr;

    tin->SetBranchAddress("JetAntiKt", &jets);
    tin->SetBranchAddress("Electron", &electrons);
    tin->SetBranchAddress("Muon", &muons);
    tin->SetBranchAddress("EFlowTrack", &tracks);
    tin->SetBranchAddress("EFlowPhoton", &photons);
    tin->SetBranchAddress("EFlowNeutralHadron", &neutrals);

    TTree *tout = new TTree("Events", "bbbb skim");

    Long64_t eventNumber;
    Float_t thrust;
    Float_t thrustAxis_x, thrustAxis_y, thrustAxis_z;
    Float_t jet_pt[4], jet_eta[4], jet_phi[4];

    tout->Branch("eventNumber", &eventNumber, "eventNumber/L");
    tout->Branch("thrust", &thrust, "thrust/F");

    tout->Branch("thrustAxis_x", &thrustAxis_x, "thrustAxis_x/F");
    tout->Branch("thrustAxis_y", &thrustAxis_y, "thrustAxis_y/F");
    tout->Branch("thrustAxis_z", &thrustAxis_z, "thrustAxis_z/F");

    tout->Branch("jet_pt", jet_pt, "jet_pt[4]/F");
    tout->Branch("jet_eta", jet_eta, "jet_eta[4]/F");
    tout->Branch("jet_phi", jet_phi, "jet_phi[4]/F");

    Long64_t n = tin->GetEntries();

    for (Long64_t i = 0; i < n; ++i) {
        tin->GetEntry(i);
	eventNumber = i;

        // lepton veto
        if (electrons->GetEntries() + muons->GetEntries() != 0) continue;

        // collect b-jets
        std::vector<Jet*> bjets;
        for (int j = 0; j < jets->GetEntries(); ++j) {
            Jet *jet = (Jet*)jets->At(j);
            if (jet->BTag) bjets.push_back(jet);
        }
        if (bjets.size() != 4) continue;

        // sort by pT
        std::sort(bjets.begin(), bjets.end(),
                  [](Jet* a, Jet* b){ return a->PT > b->PT; });

        for (int k = 0; k < 4; ++k) {
            jet_pt[k]  = bjets[k]->PT;
            jet_eta[k] = bjets[k]->Eta;
            jet_phi[k] = bjets[k]->Phi;
        }

        // thrust momenta
        std::vector<TLorentzVector> momenta;
        for (int j=0;j<tracks->GetEntries();j++)
            momenta.push_back(((Track*)tracks->At(j))->P4());
        for (int j=0;j<photons->GetEntries();j++)
            momenta.push_back(((Tower*)photons->At(j))->P4());
        for (int j=0;j<neutrals->GetEntries();j++)
            momenta.push_back(((Tower*)neutrals->At(j))->P4());

        TLorentzVector axis;
        thrust = findThrust(momenta, axis);
	thrustAxis_x = axis.X();
        thrustAxis_y = axis.Y();
        thrustAxis_z = axis.Z();

        tout->Fill();
    }

    fout->Write();
    fout->Close();
}

