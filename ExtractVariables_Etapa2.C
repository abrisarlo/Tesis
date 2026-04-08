// =============================================================
//  ExtractVariables_Etapa2.C
//
//  LEE: output_test.root  (arbol de Delphes)
//  ESCRIBE: variables_etapa2.root
//
//  SIRVE PARA SIGNAL Y BACKGROUNDS hasta donde entiendo
//  Solo cambiar la variable inputFile para cada proceso.
//
//  COMO CORRERLO (desde la carpeta ~/Documentos/Facu/Tesis):
//      conda activate env-root (En el enviroment en el que tengas root)
//      Anda a la carpeta donde tengas delphes
//      root -l /Documentos/Facu/Tesis/ExtractVariables_Etapa2.C
// =============================================================

#ifdef __CLING__
R__LOAD_LIBRARY(/home/abri/delphes/libDelphes.so)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "TVectorD.h"
#include "TVector3.h"
#include <iostream>
#include <vector>
#include <limits>
#endif

using namespace std;

// ---------------------------------------------------------------
// cos(theta) exacto a partir de eta: cos(theta) = cos(2*atan(e^{-eta}))
// ---------------------------------------------------------------
float calcCosTheta(float eta) {
    float eEta = pow(TMath::E(), -eta);
    float theta = 2 * TMath::ATan(eEta);
    return (float) TMath::Cos(theta);
}

// ---------------------------------------------------------------
// Forma del evento: sphericity o aplanarity
// Usa el tensor de momentos normalizado con los 4 jets
// ---------------------------------------------------------------
float findEventShape(TLorentzVector j1, TLorentzVector j2,
                     TLorentzVector j3, TLorentzVector j4, string shape) {
    vector<TVector3> vecs = {
        TVector3(j1.Px(), j1.Py(), j1.Pz()),
        TVector3(j2.Px(), j2.Py(), j2.Pz()),
        TVector3(j3.Px(), j3.Py(), j3.Pz()),
        TVector3(j4.Px(), j4.Py(), j4.Pz())
    };
    TMatrixDSym M(3);
    for (auto& p : vecs)
        for (int k = 0; k < 3; k++)
            for (int m = 0; m <= k; m++)
                M[k][m] += p[k] * p[m];
    M *= 1.0 / (M[0][0] + M[1][1] + M[2][2]);
    TMatrixDSymEigen eigen(M);
    TVectorD ev = eigen.GetEigenValues();
    // eigenvalues ordenados de mayor a menor: ev[0] >= ev[1] >= ev[2]
    float l1 = ev[0], l2 = ev[1], l3 = ev[2];
    float sum = l1 + l2 + l3;
    if (shape == "sphericity") return 1.5f * (l2 + l3) / sum;
    if (shape == "aplanarity") return 1.5f * l3 / sum;
    return 0.f;
}

// ---------------------------------------------------------------
// Thrust: maximiza sum|p.n| / sum|p| sobre todos los vectores unitarios n
// Se hace por barrido angular (igual que en el original)
// ---------------------------------------------------------------
float findThrust(vector<TLorentzVector>& momenta, TLorentzVector& thrustAxis) {
    float maxThrust = 0.0;
    for (float theta = 0.0; theta < TMath::Pi(); theta += 0.1) {
        for (float phi = 0.0; phi < 2.0 * TMath::Pi(); phi += 0.1) {
            TLorentzVector axis(TMath::Sin(theta) * TMath::Cos(phi),
                                TMath::Sin(theta) * TMath::Sin(phi),
                                TMath::Cos(theta), 0.0);
            float sumP = 0.0, sumP_par = 0.0;
            for (auto& p : momenta) {
                float mag = TMath::Sqrt(axis.Px()*axis.Px() +
                                       axis.Py()*axis.Py() +
                                       axis.Pz()*axis.Pz());
                sumP += p.P();
                sumP_par += TMath::Abs(p.Dot(axis) / mag);
            }
            float t = (sumP > 0) ? sumP_par / sumP : 0;
            if (t > maxThrust) { maxThrust = t; thrustAxis = axis; }
        }
    }
    return maxThrust;
}

// ---------------------------------------------------------------
// Pairing de jets que minimiza chi2 respecto a masa del H (125 GeV)
// Devuelve los dos TLorentzVector de los pares optimos
// ---------------------------------------------------------------
void findJetPairsHH(TLorentzVector j1, TLorentzVector j2,
                    TLorentzVector j3, TLorentzVector j4,
                    TLorentzVector& pair1, TLorentzVector& pair2) {
    TLorentzVector p12 = j1+j2, p34 = j3+j4;
    TLorentzVector p13 = j1+j3, p24 = j2+j4;
    TLorentzVector p14 = j1+j4, p23 = j2+j3;
    double chi12 = pow(p12.M()-125,2) + pow(p34.M()-125,2);
    double chi13 = pow(p13.M()-125,2) + pow(p24.M()-125,2);
    double chi14 = pow(p14.M()-125,2) + pow(p23.M()-125,2);
    double minChi = TMath::Min(TMath::Min(chi12, chi13), chi14);
    if      (minChi == chi12) { pair1 = p12; pair2 = p34; }
    else if (minChi == chi13) { pair1 = p13; pair2 = p24; }
    else                      { pair1 = p14; pair2 = p23; }
}

// ---------------------------------------------------------------
// Pairing de jets que minimiza chi2 respecto a masa del Z (90 GeV)
// ---------------------------------------------------------------
void findJetPairsZZ(TLorentzVector j1, TLorentzVector j2,
                    TLorentzVector j3, TLorentzVector j4,
                    TLorentzVector& pair1, TLorentzVector& pair2,
                    float& distZ1, float& distZ2) {
    TLorentzVector p12 = j1+j2, p34 = j3+j4;
    TLorentzVector p13 = j1+j3, p24 = j2+j4;
    TLorentzVector p14 = j1+j4, p23 = j2+j3;
    double chi12 = pow(p12.M()-90,2) + pow(p34.M()-90,2);
    double chi13 = pow(p13.M()-90,2) + pow(p24.M()-90,2);
    double chi14 = pow(p14.M()-90,2) + pow(p23.M()-90,2);
    double minChi = TMath::Min(TMath::Min(chi12, chi13), chi14);
    if      (minChi == chi12) { pair1 = p12; pair2 = p34; }
    else if (minChi == chi13) { pair1 = p13; pair2 = p24; }
    else                      { pair1 = p14; pair2 = p23; }
    distZ1 = (float) fabs(pair1.M() - 90.0);
    distZ2 = (float) fabs(pair2.M() - 90.0);
}

// ---------------------------------------------------------------
// Masa minima entre todos los pares posibles de jets
// ---------------------------------------------------------------
float findMinJetM(TLorentzVector j1, TLorentzVector j2,
                  TLorentzVector j3, TLorentzVector j4) {
    float m = TMath::Min((float)(j1+j2).M(), (float)(j3+j4).M());
    m = TMath::Min(m, (float)(j1+j3).M());
    m = TMath::Min(m, (float)(j2+j4).M());
    m = TMath::Min(m, (float)(j1+j4).M());
    m = TMath::Min(m, (float)(j2+j3).M());
    return m;
}

// ---------------------------------------------------------------
// FUNCION PRINCIPAL
// ---------------------------------------------------------------
void ExtractVariables_Etapa2() {

    // ----------------------------------------------------------
    // Archivos de entrada y salida
    // Cambiar inputFile para correr sobre backgrounds
    // ----------------------------------------------------------
    const char *inputFile  = "/home/abri/Documentos/Facu/Tesis/output_test.root";
    const char *outputFile = "/home/abri/Documentos/Facu/Tesis/variables_etapa2.root";

    // ----------------------------------------------------------
    // Abrir el arbol de Delphes
    // ----------------------------------------------------------
    TChain chain("Delphes");
    chain.Add(inputFile);

    ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
    Long64_t numberOfEntries = treeReader->GetEntries();
    cout << "Eventos en el archivo: " << numberOfEntries << endl;

    // ----------------------------------------------------------
    // Branches de entrada
    // ----------------------------------------------------------
    TClonesArray *branchJet            = treeReader->UseBranch("Jet");
    TClonesArray *branchMissingET      = treeReader->UseBranch("MissingET");
    TClonesArray *branchElectron       = treeReader->UseBranch("Electron");
    TClonesArray *branchMuon           = treeReader->UseBranch("Muon");
    TClonesArray *branchEFlowTrack     = treeReader->UseBranch("EFlowTrack");
    TClonesArray *branchEFlowPhoton    = treeReader->UseBranch("EFlowPhoton");
    TClonesArray *branchEFlowNeutHad   = treeReader->UseBranch("EFlowNeutralHadron");

    // ----------------------------------------------------------
    // Archivo y arbol de salida
    // ----------------------------------------------------------
    TFile *outputTFile = new TFile(outputFile, "recreate");
    TTree *outTree = new TTree("Variables", "Variables fisicas para BDT - Etapa 2");

    // ----------------------------------------------------------
    // === ETAPA 1: variables originales ===
    // ----------------------------------------------------------
    float cosThetaJet1, cosThetaJet2, cosThetaJet3, cosThetaJet4;
    float ptJet1, ptJet2, ptJet3, ptJet4, sumPt;
    float massJet1, massJet2, massJet3, massJet4;
    float missingET;

    outTree->Branch("cosThetaJet1", &cosThetaJet1, "cosThetaJet1/F");
    outTree->Branch("cosThetaJet2", &cosThetaJet2, "cosThetaJet2/F");
    outTree->Branch("cosThetaJet3", &cosThetaJet3, "cosThetaJet3/F");
    outTree->Branch("cosThetaJet4", &cosThetaJet4, "cosThetaJet4/F");
    outTree->Branch("ptJet1",  &ptJet1,  "ptJet1/F");
    outTree->Branch("ptJet2",  &ptJet2,  "ptJet2/F");
    outTree->Branch("ptJet3",  &ptJet3,  "ptJet3/F");
    outTree->Branch("ptJet4",  &ptJet4,  "ptJet4/F");
    outTree->Branch("sumPt",   &sumPt,   "sumPt/F");
    outTree->Branch("massJet1",&massJet1,"massJet1/F");
    outTree->Branch("massJet2",&massJet2,"massJet2/F");
    outTree->Branch("massJet3",&massJet3,"massJet3/F");
    outTree->Branch("massJet4",&massJet4,"massJet4/F");
    outTree->Branch("missingET",&missingET,"missingET/F");

    // ----------------------------------------------------------
    // === ETAPA 2: variables nuevas ===
    // ----------------------------------------------------------

    // -- Forma del evento --
    float aplanarity, sphericity, thrust;
    outTree->Branch("aplanarity", &aplanarity, "aplanarity/F");
    outTree->Branch("sphericity", &sphericity, "sphericity/F");
    outTree->Branch("thrust",     &thrust,     "thrust/F");

    // -- Masa invariante de pares (pairing por chi2 minimo a mH=125) --
    float invMassB1, invMassB2;
    float invMass4Jets;
    float minJetM;
    float jetPairMassDelta;
    outTree->Branch("invMassB1",       &invMassB1,       "invMassB1/F");
    outTree->Branch("invMassB2",       &invMassB2,       "invMassB2/F");
    outTree->Branch("invMass4Jets",    &invMass4Jets,    "invMass4Jets/F");
    outTree->Branch("minJetM",         &minJetM,         "minJetM/F");
    outTree->Branch("jetPairMassDelta",&jetPairMassDelta,"jetPairMassDelta/F");

    // -- Masa invariante de pares recosntruyendo Z (para discriminar ZZ) --
    float invMassZZ1, invMassZZ2;
    float distanceZ1MinChiSquaredZZMass, distanceZ2MinChiSquaredZZMass;
    outTree->Branch("invMassZZ1",    &invMassZZ1,    "invMassZZ1/F");
    outTree->Branch("invMassZZ2",    &invMassZZ2,    "invMassZZ2/F");
    outTree->Branch("distanceZ1MinChiSquaredZZMass", &distanceZ1MinChiSquaredZZMass,
                    "distanceZ1MinChiSquaredZZMass/F");
    outTree->Branch("distanceZ2MinChiSquaredZZMass", &distanceZ2MinChiSquaredZZMass,
                    "distanceZ2MinChiSquaredZZMass/F");

    // -- Boost longitudinal de cada jet y del sistema --
    float boostB1, boostB2, boostB3, boostB4, boostSystem;
    outTree->Branch("boostB1",    &boostB1,    "boostB1/F");
    outTree->Branch("boostB2",    &boostB2,    "boostB2/F");
    outTree->Branch("boostB3",    &boostB3,    "boostB3/F");
    outTree->Branch("boostB4",    &boostB4,    "boostB4/F");
    outTree->Branch("boostSystem",&boostSystem,"boostSystem/F");

    // -- Asimetrias de energia y pT entre los pares de jets --
    float eAssymetry, pTAssymetry;
    outTree->Branch("eAssymetry",  &eAssymetry,  "eAssymetry/F");
    outTree->Branch("pTAssymetry", &pTAssymetry, "pTAssymetry/F");

    // -- Distancia angular DeltaR entre los dos pares de jets --
    float deltaRJetPairs;
    outTree->Branch("deltaRJetPairs",&deltaRJetPairs,"deltaRJetPairs/F");

    // -- Constituyentes de jets --
    float constSizeB1, constSizeB2, constSizeB3, constSizeB4, minConstSize;
    float jetB1NCharged, jetB2NCharged, jetB3NCharged, jetB4NCharged;
    float jetB1NNeutrals, jetB2NNeutrals, jetB3NNeutrals, jetB4NNeutrals;
    float jetNObjects, minJetNObjects;
    float nParticles;
    outTree->Branch("constSizeB1",   &constSizeB1,   "constSizeB1/F");
    outTree->Branch("constSizeB2",   &constSizeB2,   "constSizeB2/F");
    outTree->Branch("constSizeB3",   &constSizeB3,   "constSizeB3/F");
    outTree->Branch("constSizeB4",   &constSizeB4,   "constSizeB4/F");
    outTree->Branch("minConstSize",  &minConstSize,  "minConstSize/F");
    outTree->Branch("jetB1NCharged", &jetB1NCharged, "jetB1NCharged/F");
    outTree->Branch("jetB2NCharged", &jetB2NCharged, "jetB2NCharged/F");
    outTree->Branch("jetB3NCharged", &jetB3NCharged, "jetB3NCharged/F");
    outTree->Branch("jetB4NCharged", &jetB4NCharged, "jetB4NCharged/F");
    outTree->Branch("jetB1NNeutrals",&jetB1NNeutrals,"jetB1NNeutrals/F");
    outTree->Branch("jetB2NNeutrals",&jetB2NNeutrals,"jetB2NNeutrals/F");
    outTree->Branch("jetB3NNeutrals",&jetB3NNeutrals,"jetB3NNeutrals/F");
    outTree->Branch("jetB4NNeutrals",&jetB4NNeutrals,"jetB4NNeutrals/F");
    outTree->Branch("jetNObjects",   &jetNObjects,   "jetNObjects/F");
    outTree->Branch("minJetNObjects",&minJetNObjects,"minJetNObjects/F");
    outTree->Branch("nParticles",    &nParticles,    "nParticles/F");

    // -- Variables de Durham jet clustering exclusivo --
    // Se leen del primer jet del branch (Delphes las guarda ahi)
    float exclYmerge12, exclYmerge23, exclYmerge34, exclYmerge45, exclYmerge56;
    outTree->Branch("exclYmerge12",&exclYmerge12,"exclYmerge12/F");
    outTree->Branch("exclYmerge23",&exclYmerge23,"exclYmerge23/F");
    outTree->Branch("exclYmerge34",&exclYmerge34,"exclYmerge34/F");
    outTree->Branch("exclYmerge45",&exclYmerge45,"exclYmerge45/F");
    outTree->Branch("exclYmerge56",&exclYmerge56,"exclYmerge56/F");

    // ----------------------------------------------------------
    // Contadores para el resumen final
    // ----------------------------------------------------------
    int contTotal = 0, contPasaron = 0, contVeto = 0;

    // ----------------------------------------------------------
    // LOOP PRINCIPAL
    // ----------------------------------------------------------
    for (Long64_t entry = 0; entry < numberOfEntries; entry++) {

        treeReader->ReadEntry(entry);
        contTotal++;

        // ---- CORTE 1: al menos 4 jets ----
        int nJets = branchJet->GetEntries();
        if (nJets < 4) continue;

        // ---- CORTE 2: veto de leptones ----
        int nElectrons = branchElectron->GetEntries();
        int nMuons     = branchMuon->GetEntries();
        if (nElectrons > 0 || nMuons > 0) { contVeto++; continue; }

        // ---- Leer los 4 jets con mayor pT (Delphes los ordena por pT) ----
        Jet *jet1 = (Jet*) branchJet->At(0);
        Jet *jet2 = (Jet*) branchJet->At(1);
        Jet *jet3 = (Jet*) branchJet->At(2);
        Jet *jet4 = (Jet*) branchJet->At(3);

        TLorentzVector j1 = jet1->P4();
        TLorentzVector j2 = jet2->P4();
        TLorentzVector j3 = jet3->P4();
        TLorentzVector j4 = jet4->P4();

        // =====================================================
        // ETAPA 1: variables originales
        // =====================================================
        cosThetaJet1 = calcCosTheta(jet1->Eta);
        cosThetaJet2 = calcCosTheta(jet2->Eta);
        cosThetaJet3 = calcCosTheta(jet3->Eta);
        cosThetaJet4 = calcCosTheta(jet4->Eta);

        ptJet1 = jet1->PT;
        ptJet2 = jet2->PT;
        ptJet3 = jet3->PT;
        ptJet4 = jet4->PT;
        sumPt  = ptJet1 + ptJet2 + ptJet3 + ptJet4;

        massJet1 = jet1->Mass;
        massJet2 = jet2->Mass;
        massJet3 = jet3->Mass;
        massJet4 = jet4->Mass;

        MissingET *met = (MissingET*) branchMissingET->At(0);
        missingET = met->MET;

        // =====================================================
        // ETAPA 2: variables nuevas
        // =====================================================

        // ---- Forma del evento ----
        aplanarity = findEventShape(j1, j2, j3, j4, "aplanarity");
        sphericity = findEventShape(j1, j2, j3, j4, "sphericity");

        // ---- Thrust (sobre todos los objetos EFlow del evento) ----
        vector<TLorentzVector> momenta;
        for (int i = 0; i < branchEFlowTrack->GetEntries(); i++) {
            Track *t = (Track*) branchEFlowTrack->At(i);
            momenta.push_back(t->P4());
        }
        for (int i = 0; i < branchEFlowPhoton->GetEntries(); i++) {
            Tower *p = (Tower*) branchEFlowPhoton->At(i);
            momenta.push_back(p->P4());
        }
        for (int i = 0; i < branchEFlowNeutHad->GetEntries(); i++) {
            Tower *n = (Tower*) branchEFlowNeutHad->At(i);
            momenta.push_back(n->P4());
        }
        TLorentzVector thrustAxis;
        thrust = findThrust(momenta, thrustAxis);

        // ---- Masa invariante de pares (chi2 minimo a mH=125) ----
        TLorentzVector pairHH1, pairHH2;
        findJetPairsHH(j1, j2, j3, j4, pairHH1, pairHH2);
        invMassB1 = (float) pairHH1.M();
        invMassB2 = (float) pairHH2.M();
        invMass4Jets = (float) (j1 + j2 + j3 + j4).M();
        minJetM = findMinJetM(j1, j2, j3, j4);
        jetPairMassDelta = pow(invMassB1 - 125.2f, 2) + pow(invMassB2 - 125.2f, 2);

        // ---- Masa invariante de pares reconstruyendo Z ----
        TLorentzVector pairZZ1, pairZZ2;
        findJetPairsZZ(j1, j2, j3, j4, pairZZ1, pairZZ2,
                       distanceZ1MinChiSquaredZZMass,
                       distanceZ2MinChiSquaredZZMass);
        invMassZZ1 = (float) pairZZ1.M();
        invMassZZ2 = (float) pairZZ2.M();

        // ---- Boost longitudinal: beta_z = pz/E ----
        TLorentzVector totalSys = j1 + j2 + j3 + j4;
        boostB1    = (j1.E() > 0) ? (float)(j1.Pz()       / j1.E())        : 0.f;
        boostB2    = (j2.E() > 0) ? (float)(j2.Pz()       / j2.E())        : 0.f;
        boostB3    = (j3.E() > 0) ? (float)(j3.Pz()       / j3.E())        : 0.f;
        boostB4    = (j4.E() > 0) ? (float)(j4.Pz()       / j4.E())        : 0.f;
        boostSystem = (totalSys.E() > 0) ? (float)(totalSys.Pz() / totalSys.E()) : 0.f;

        // ---- Asimetrias entre pares HH ----
        float E1 = pairHH1.E(), E2 = pairHH2.E();
        eAssymetry  = ((E1 + E2) > 0)  ? fabs(E1 - E2) / (E1 + E2)                : 0.f;
        float PT1 = pairHH1.Pt(), PT2 = pairHH2.Pt();
        pTAssymetry = ((PT1 + PT2) > 0) ? fabs(PT1 - PT2) / (PT1 + PT2)           : 0.f;

        // ---- DeltaR entre los dos pares de jets ----
        deltaRJetPairs = (float) pairHH1.DeltaR(pairHH2);

        // ---- Constituyentes de cada jet ----
        TRefArray* c1 = &(jet1->Constituents);
        TRefArray* c2 = &(jet2->Constituents);
        TRefArray* c3 = &(jet3->Constituents);
        TRefArray* c4 = &(jet4->Constituents);
        constSizeB1 = (float) c1->GetSize();
        constSizeB2 = (float) c2->GetSize();
        constSizeB3 = (float) c3->GetSize();
        constSizeB4 = (float) c4->GetSize();
        minConstSize = TMath::Min(TMath::Min(constSizeB1, constSizeB2),
                                 TMath::Min(constSizeB3, constSizeB4));

        jetB1NCharged  = (float) jet1->NCharged;
        jetB2NCharged  = (float) jet2->NCharged;
        jetB3NCharged  = (float) jet3->NCharged;
        jetB4NCharged  = (float) jet4->NCharged;
        jetB1NNeutrals = (float) jet1->NNeutrals;
        jetB2NNeutrals = (float) jet2->NNeutrals;
        jetB3NNeutrals = (float) jet3->NNeutrals;
        jetB4NNeutrals = (float) jet4->NNeutrals;

        float nObj1 = jetB1NCharged + jetB1NNeutrals;
        float nObj2 = jetB2NCharged + jetB2NNeutrals;
        float nObj3 = jetB3NCharged + jetB3NNeutrals;
        float nObj4 = jetB4NCharged + jetB4NNeutrals;
        jetNObjects    = nObj1 + nObj2 + nObj3 + nObj4;
        minJetNObjects = TMath::Min(TMath::Min(nObj1, nObj2),
                                   TMath::Min(nObj3, nObj4));

        // nParticles: numero de objetos EFlow totales del evento
        nParticles = (float)(branchEFlowTrack->GetEntries()
                           + branchEFlowPhoton->GetEntries()
                           + branchEFlowNeutHad->GetEntries());

        // ---- Variables de Durham jet clustering exclusivo ----
        // Delphes guarda los ExclYmergeXX en el primer jet del branch
        exclYmerge12 = jet1->ExclYmerge12;
        exclYmerge23 = jet1->ExclYmerge23;
        exclYmerge34 = jet1->ExclYmerge34;
        exclYmerge45 = jet1->ExclYmerge45;
        exclYmerge56 = jet1->ExclYmerge56;

        // ---- Guardar evento ----
        outTree->Fill();
        contPasaron++;

        if (entry % 1000 == 0)
            cout << "Evento " << entry << " / " << numberOfEntries
                 << "  |  guardados: " << contPasaron << endl;
    }

    // ----------------------------------------------------------
    // Guardar y cerrar
    // ----------------------------------------------------------
    outputTFile->cd();
    outTree->Write();
    outputTFile->Close();

    cout << endl;
    cout << "========================================" << endl;
    cout << "Archivo guardado: " << outputFile          << endl;
    cout << "Eventos leidos:   " << contTotal            << endl;
    cout << "Rechazados veto leptones: " << contVeto     << endl;
    cout << "Eventos guardados: " << contPasaron          << endl;
    cout << "========================================" << endl;
    cout << endl;
    cout << "Para verificar:" << endl;
    cout << "  root -l variables_etapa2.root"           << endl;
    cout << "  Variables->Print()"                      << endl;
    cout << "  Variables->Draw(\"invMassB1\")"          << endl;
    cout << "  Variables->Draw(\"aplanarity\")"         << endl;

    delete treeReader;
}
