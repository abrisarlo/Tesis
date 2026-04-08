// =============================================================
//  ExtractVariables_Etapa1.C
//
//  LEE: output_test.root  (arbol de Delphes)
//  ESCRIBE: variables_etapa1.root
//
//  Variables que defino: cos(theta), pT, masa de los 4 jets y MissingET.
//  Vamos a ir sumando variables en etapas siguientes.
//
//  COMO CORRERLO (desde la carpeta ~/Documentos/Facu/Tesis):
//      conda activate env-root
//      root -l ExtractVariables_Etapa1.C
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
#include <iostream>
#endif

using namespace std;

// ---------------------------------------------------------------
// cos(theta) a partir de eta: formula exacta.
// Eta Sale de los Jets
// ---------------------------------------------------------------
float calcCosTheta(float eta) {
    return (float) tanh((double)eta);
}

// ---------------------------------------------------------------
// FUNCION PRINCIPAL
// ---------------------------------------------------------------
void ExtractVariables_Etapa1() {

    // ----------------------------------------------------------
    // Archivos de entrada y salida
    // ----------------------------------------------------------
    const char *inputFile  = "/home/abri/Documentos/Facu/Tesis/output_test.root";
    const char *outputFile = "/home/abri/Documentos/Facu/Tesis/variables_etapa1.root";

    // ----------------------------------------------------------
    // Abrir el arbol de Delphes.
    // Te crea una cadena donde pone tu archivo de salida. Ademas pone un lector de columnas y un contador de eventos.
    // ----------------------------------------------------------
    TChain chain("Delphes");
    chain.Add(inputFile);

    ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
    Long64_t numberOfEntries = treeReader->GetEntries();
    cout << "Eventos en el archivo: " << numberOfEntries << endl;

    // ----------------------------------------------------------
    // Declarar los branches que vamos a Leer del archivo. 
    // Osea los que vamos a usar para definir las nuevas variables.
    // ----------------------------------------------------------
    TClonesArray *branchJet       = treeReader->UseBranch("Jet");
    TClonesArray *branchMissingET = treeReader->UseBranch("MissingET");
    TClonesArray *branchElectron  = treeReader->UseBranch("Electron");
    TClonesArray *branchMuon      = treeReader->UseBranch("Muon");

    // ----------------------------------------------------------
    // Crear el archivo y el arbol de salida
    // ----------------------------------------------------------
    TFile *outputTFile = new TFile(outputFile, "recreate");
    TTree *outTree = new TTree("Variables", "Variables fisicas para BDT");

    // ----------------------------------------------------------
    // Variables que vamos a GUARDAR - Etapa 1
    // ----------------------------------------------------------

    // cos(theta) de cada jet
    float cosThetaJet1, cosThetaJet2, cosThetaJet3, cosThetaJet4;

    // pT de cada jet (momentum transverso, en GeV)
    float ptJet1, ptJet2, ptJet3, ptJet4;

    // Suma de pT de los 4 jets
    float sumPt;

    // Masa de cada jet (en GeV)
    float massJet1, massJet2, massJet3, massJet4;

    // Missing ET (energia transversa faltante)
    float missingET;

    // ----------------------------------------------------------
    // Conectar las variables al arbol de salida
    // Sintaxis: Branch("nombre_en_root", &variable, "nombre/tipo")
    // F = float (Esto es muy util :))
    // ----------------------------------------------------------
    outTree->Branch("cosThetaJet1", &cosThetaJet1, "cosThetaJet1/F");
    outTree->Branch("cosThetaJet2", &cosThetaJet2, "cosThetaJet2/F");
    outTree->Branch("cosThetaJet3", &cosThetaJet3, "cosThetaJet3/F");
    outTree->Branch("cosThetaJet4", &cosThetaJet4, "cosThetaJet4/F");

    outTree->Branch("ptJet1", &ptJet1, "ptJet1/F");
    outTree->Branch("ptJet2", &ptJet2, "ptJet2/F");
    outTree->Branch("ptJet3", &ptJet3, "ptJet3/F");
    outTree->Branch("ptJet4", &ptJet4, "ptJet4/F");
    outTree->Branch("sumPt",  &sumPt,  "sumPt/F");

    outTree->Branch("massJet1", &massJet1, "massJet1/F");
    outTree->Branch("massJet2", &massJet2, "massJet2/F");
    outTree->Branch("massJet3", &massJet3, "massJet3/F");
    outTree->Branch("massJet4", &massJet4, "massJet4/F");

    outTree->Branch("missingET", &missingET, "missingET/F");

    // ----------------------------------------------------------
    // Contadores para el resumen final
    // ----------------------------------------------------------
    int contTotal   = 0;
    int contPasaron = 0;
    int contVeto    = 0;

    // ----------------------------------------------------------
    // LOOP PRINCIPAL: recorrer cada evento uno por uno
    // ----------------------------------------------------------
    for (Long64_t entry = 0; entry < numberOfEntries; entry++) {

        // Cargar los datos del evento numero "entry"
        treeReader->ReadEntry(entry);
        contTotal++;

        // ---- CORTE 1: necesitamos al menos 4 jets ----
        int nJets = branchJet->GetEntries();
        if (nJets < 4) continue;  // saltar al siguiente evento

        // ---- CORTE 2: veto de leptones ----
        // En el analisis HH->bbbb no deben haber electrones ni muones
        int nElectrons = branchElectron->GetEntries();
        int nMuons     = branchMuon->GetEntries();
        if (nElectrons > 0 || nMuons > 0) {
            contVeto++;
            continue;
        }

        // ---- Leer los 4 jets con mayor pT ----
        // Delphes los ordena automaticamente por pT de mayor a menor
	// Te quedas con los de mayor PT que tienen mas chances de ser de HH
        // At(0) = jet con mayor pT, At(1) = segundo mayor, etc.
        Jet *jet1 = (Jet*) branchJet->At(0);
        Jet *jet2 = (Jet*) branchJet->At(1);
        Jet *jet3 = (Jet*) branchJet->At(2);
        Jet *jet4 = (Jet*) branchJet->At(3);

        // ---- Calcular cos(theta) ----
        // Eta es la pseudorapidez; cos(theta) = tanh(eta) exactamente
        cosThetaJet1 = calcCosTheta(jet1->Eta);
        cosThetaJet2 = calcCosTheta(jet2->Eta);
        cosThetaJet3 = calcCosTheta(jet3->Eta);
        cosThetaJet4 = calcCosTheta(jet4->Eta);

        // ---- pT de cada jet (ya viene calculado por Delphes) ----
        ptJet1 = jet1->PT;
        ptJet2 = jet2->PT;
        ptJet3 = jet3->PT;
        ptJet4 = jet4->PT;
        sumPt  = ptJet1 + ptJet2 + ptJet3 + ptJet4;

        // ---- Masa de cada jet ----
        massJet1 = jet1->Mass;
        massJet2 = jet2->Mass;
        massJet3 = jet3->Mass;
        massJet4 = jet4->Mass;

        // ---- Missing ET ----
        // branchMissingET->At(0) porque siempre hay exactamente 1 objeto MET
        MissingET *met = (MissingET*) branchMissingET->At(0);
        missingET = met->MET;

        // ---- Guardar este evento en el arbol de salida ----
        // Fill() agrega una fila nueva con todos los valores actuales
        outTree->Fill();
        contPasaron++;

        // Mostrar progreso cada 1000 eventos
        if (entry % 1000 == 0) {
            cout << "Evento " << entry << " / " << numberOfEntries
                 << "  |  guardados: " << contPasaron << endl;
        }
    }

    // ----------------------------------------------------------
    // Guardar el arbol en el archivo y cerrar
    // Sin este paso el archivo queda vacio!
    // ----------------------------------------------------------
    outputTFile->cd();
    outTree->Write();
    outputTFile->Close();

    // Resumen final
    cout << endl;
    cout << "========================================" << endl;
    cout << "Archivo guardado: " << outputFile          << endl;
    cout << "Eventos leidos:   " << contTotal            << endl;
    cout << "Rechazados por veto de leptones: " << contVeto << endl;
    cout << "Eventos guardados: " << contPasaron          << endl;
    cout << "========================================" << endl;
    cout << endl;
    cout << "Para verificar, abri ROOT y ejecuta:"           << endl;
    cout << "  root -l variables_etapa1.root"                << endl;
    cout << "  Variables->Print()         <- lista columnas" << endl;
    cout << "  Variables->Draw(\"ptJet1\")  <- grafica pT"    << endl;

    delete treeReader;
}
