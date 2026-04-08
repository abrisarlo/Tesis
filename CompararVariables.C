// =============================================================
//  CompararVariables.C
//
//  Compara las variables del output del FSRGammaGamma (TreeS)
//  con las del output de ExtractVariables_Etapa2 (Variables)
//
//  COMO CORRERLO (desde ~/delphes):
//      root -l ~/Documentos/Facu/Tesis/CompararVariables.C
// =============================================================

#ifdef __CLING__
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TMath.h"
#include <iostream>
#endif

using namespace std;

void dibujar(TTree *tFSR, TTree *tMio,
             const char *varFSR, const char *varMio,
             const char *etiqueta,
             int nBins, float xMin, float xMax) {

    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    gPad->SetTopMargin(0.08);

    TH1F *hFSR = new TH1F(Form("hFSR_%s", varFSR), "", nBins, xMin, xMax);
    TH1F *hMio = new TH1F(Form("hMio_%s", varMio), "", nBins, xMin, xMax);

    tFSR->Draw(Form("%s>>hFSR_%s", varFSR, varFSR), "", "goff");
    tMio->Draw(Form("%s>>hMio_%s", varMio, varMio), "", "goff");

    if (hFSR->Integral() > 0) hFSR->Scale(1.0 / hFSR->Integral());
    if (hMio->Integral() > 0) hMio->Scale(1.0 / hMio->Integral());

    hFSR->SetLineColor(kAzure+2);
    hFSR->SetLineWidth(2);
    hFSR->SetFillColorAlpha(kAzure+2, 0.2);
    hFSR->SetFillStyle(1001);

    hMio->SetLineColor(kRed+1);
    hMio->SetLineWidth(2);
    hMio->SetLineStyle(2);

    float yMax = TMath::Max(hFSR->GetMaximum(), hMio->GetMaximum()) * 1.35;
    hFSR->SetMaximum(yMax);
    hFSR->SetMinimum(0);

    hFSR->GetXaxis()->SetTitle(etiqueta);
    hFSR->GetYaxis()->SetTitle("Fraccion de eventos");
    hFSR->GetXaxis()->SetTitleSize(0.055);
    hFSR->GetYaxis()->SetTitleSize(0.050);
    hFSR->GetXaxis()->SetLabelSize(0.045);
    hFSR->GetYaxis()->SetLabelSize(0.045);
    hFSR->GetXaxis()->SetTitleOffset(1.1);
    hFSR->GetYaxis()->SetTitleOffset(1.3);
    hFSR->SetTitle("");

    hFSR->Draw("HIST");
    hMio->Draw("HIST SAME");

    TLegend *leg = new TLegend(0.50, 0.75, 0.93, 0.93);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.048);
    leg->AddEntry(hFSR, "Codigo Santiago", "lf");
    leg->AddEntry(hMio, "Codigo Abri", "l");
    leg->Draw();
}

void CompararVariables() {

    const char *archivoFSR = "/home/abri/Documentos/Facu/Tesis/outputTreeSHHbbbbESpreadDurham1034BSplitSampleN.root";
    const char *archivoMio = "/home/abri/Documentos/Facu/Tesis/variables_etapa2.root";

    TFile *fFSR = TFile::Open(archivoFSR);
    if (!fFSR || fFSR->IsZombie()) { cout << "ERROR abriendo FSR" << endl; return; }
    TTree *tFSR = (TTree*) fFSR->Get("TreeS");
    if (!tFSR) { cout << "ERROR: no se encontro TreeS" << endl; return; }

    TFile *fMio = TFile::Open(archivoMio);
    if (!fMio || fMio->IsZombie()) { cout << "ERROR abriendo mi archivo" << endl; return; }
    TTree *tMio = (TTree*) fMio->Get("Variables");
    if (!tMio) { cout << "ERROR: no se encontro Variables" << endl; return; }

    cout << "Eventos FSR:  " << tFSR->GetEntries() << endl;
    cout << "Eventos mios: " << tMio->GetEntries()  << endl;

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gROOT->ForceStyle();

    // ==========================================================
    // CANVAS 1: cos(theta) y pT de los 4 jets
    // ==========================================================
    TCanvas *c1 = new TCanvas("c1", "cosTheta y pT", 1600, 800);
    c1->Divide(4, 2);

    c1->cd(1); dibujar(tFSR, tMio, "cosThetaB1", "cosThetaJet1", "cos(#theta) jet1", 50,  -1.0,   1.0);
    c1->cd(2); dibujar(tFSR, tMio, "cosThetaB2", "cosThetaJet2", "cos(#theta) jet2", 50,  -1.0,   1.0);
    c1->cd(3); dibujar(tFSR, tMio, "cosThetaB3", "cosThetaJet3", "cos(#theta) jet3", 50,  -1.0,   1.0);
    c1->cd(4); dibujar(tFSR, tMio, "cosThetaB4", "cosThetaJet4", "cos(#theta) jet4", 50,  -1.0,   1.0);
    c1->cd(5); dibujar(tFSR, tMio, "jetB1Pt",    "ptJet1",       "pT jet1 [GeV]",    50,  13.0, 185.0);
    c1->cd(6); dibujar(tFSR, tMio, "jetB2Pt",    "ptJet2",       "pT jet2 [GeV]",    50,   8.0, 160.0);
    c1->cd(7); dibujar(tFSR, tMio, "jetB3Pt",    "ptJet3",       "pT jet3 [GeV]",    50,   0.0,  90.0);
    c1->cd(8); dibujar(tFSR, tMio, "jetB4Pt",    "ptJet4",       "pT jet4 [GeV]",    50,   0.0,  90.0);
    c1->SaveAs("/home/abri/Documentos/Facu/Tesis/comp_costheta_pt.pdf");

    // ==========================================================
    // CANVAS 2: Masa de jets, sumPt, missingET
    // ==========================================================
    TCanvas *c2 = new TCanvas("c2", "Masas jets", 1600, 800);
    c2->Divide(4, 2);

    c2->cd(1); dibujar(tFSR, tMio, "jetB1M",      "massJet1",     "Masa jet1 [GeV]",        50,   0.0, 115.0);
    c2->cd(2); dibujar(tFSR, tMio, "jetB2M",      "massJet2",     "Masa jet2 [GeV]",        50,   0.0, 112.0);
    c2->cd(3); dibujar(tFSR, tMio, "jetB3M",      "massJet3",     "Masa jet3 [GeV]",        50,   0.0,  78.0);
    c2->cd(4); dibujar(tFSR, tMio, "jetB4M",      "massJet4",     "Masa jet4 [GeV]",        50,   0.0,  76.0);
    c2->cd(5); dibujar(tFSR, tMio, "sumPt",       "sumPt",        "#SigmapT [GeV]",         50,  40.0, 385.0);
    c2->cd(6); dibujar(tFSR, tMio, "missingET",   "missingET",    "Missing ET [GeV]",       50,   0.0, 125.0);
    c2->cd(7); dibujar(tFSR, tMio, "minJetM",     "minJetM",      "Masa inv. min [GeV]",    50,  12.0, 130.0);
    c2->cd(8); dibujar(tFSR, tMio, "invMass4Jets","invMass4Jets",  "Masa inv. 4 jets [GeV]", 50,  65.0, 400.0);
    c2->SaveAs("/home/abri/Documentos/Facu/Tesis/comp_masas_jets.pdf");

    // ==========================================================
    // CANVAS 3: Masa invariante de pares H y Z
    // ==========================================================
    TCanvas *c3 = new TCanvas("c3", "Masas invariantes pares", 1600, 800);
    c3->Divide(3, 2);

    c3->cd(1); dibujar(tFSR, tMio, "invMassB1",  "invMassB1",  "Masa inv. par H1 [GeV]",          60,  22.0, 232.0);
    c3->cd(2); dibujar(tFSR, tMio, "invMassB2",  "invMassB2",  "Masa inv. par H2 [GeV]",          60,  12.0, 215.0);
    c3->cd(3); dibujar(tFSR, tMio, "jetPairMassDelta","jetPairMassDelta","#Delta masa pares [GeV^{2}]", 50,   0.0, 19000.0);
    c3->cd(4); dibujar(tFSR, tMio, "invMassZZ1", "invMassZZ1", "Masa inv. par ZZ1 [GeV]",         60,  15.0, 265.0);
    c3->cd(5); dibujar(tFSR, tMio, "invMassZZ2", "invMassZZ2", "Masa inv. par ZZ2 [GeV]",         60,  13.0, 272.0);
    c3->cd(6); dibujar(tFSR, tMio, "distanceZ1MinChiSquaredZZMass","distanceZ1MinChiSquaredZZMass","Dist. a mZ par1 [GeV]", 50, 0.0, 175.0);
    c3->SaveAs("/home/abri/Documentos/Facu/Tesis/comp_masas_pares.pdf");

    // ==========================================================
    // CANVAS 4: Forma del evento
    // ==========================================================
    TCanvas *c4 = new TCanvas("c4", "Forma del evento", 1200, 400);
    c4->Divide(3, 1);

    c4->cd(1); dibujar(tFSR, tMio, "aplanarity", "aplanarity", "Aplanarity", 50,  0.0,  0.5);
    c4->cd(2); dibujar(tFSR, tMio, "sphericity", "sphericity", "Sphericity", 50,  0.0,  1.0);
    c4->cd(3); dibujar(tFSR, tMio, "thrust",     "thrust",     "Thrust",     50,  0.5,  1.0);
    c4->SaveAs("/home/abri/Documentos/Facu/Tesis/comp_eventshape.pdf");

    // ==========================================================
    // CANVAS 5: Boost y asimetrias
    // ==========================================================
    TCanvas *c5 = new TCanvas("c5", "Boost y asimetrias", 1600, 800);
    c5->Divide(4, 2);

    c5->cd(1); dibujar(tFSR, tMio, "boostB1",       "boostB1",       "Boost jet1 (pz/E)",    50, -1.0,  1.0);
    c5->cd(2); dibujar(tFSR, tMio, "boostB2",       "boostB2",       "Boost jet2 (pz/E)",    50, -1.0,  1.0);
    c5->cd(3); dibujar(tFSR, tMio, "boostB3",       "boostB3",       "Boost jet3 (pz/E)",    50, -1.0,  1.0);
    c5->cd(4); dibujar(tFSR, tMio, "boostB4",       "boostB4",       "Boost jet4 (pz/E)",    50, -1.0,  1.0);
    c5->cd(5); dibujar(tFSR, tMio, "boostSystem",   "boostSystem",   "Boost sistema (pz/E)", 50, -0.8,  0.8);
    c5->cd(6); dibujar(tFSR, tMio, "eAssymetry",    "eAssymetry",    "Asimetria energia",    50,  0.0,  0.8);
    c5->cd(7); dibujar(tFSR, tMio, "pTAssymetry",   "pTAssymetry",   "Asimetria pT",         50,  0.0,  1.0);
    c5->cd(8); dibujar(tFSR, tMio, "deltaRJetPairs","deltaRJetPairs","#DeltaR pares",         50,  0.0, 10.0);
    c5->SaveAs("/home/abri/Documentos/Facu/Tesis/comp_boost_asimetrias.pdf");

    // ==========================================================
    // CANVAS 6: Constituyentes de jets
    // ==========================================================
    TCanvas *c6 = new TCanvas("c6", "Constituyentes", 1600, 800);
    c6->Divide(4, 2);

    c6->cd(1); dibujar(tFSR, tMio, "constSizeB1",   "constSizeB1",   "Constituyentes jet1",   50, 19.0,  72.0);
    c6->cd(2); dibujar(tFSR, tMio, "constSizeB2",   "constSizeB2",   "Constituyentes jet2",   50, 24.0,  67.0);
    c6->cd(3); dibujar(tFSR, tMio, "constSizeB3",   "constSizeB3",   "Constituyentes jet3",   50, 13.0,  59.0);
    c6->cd(4); dibujar(tFSR, tMio, "constSizeB4",   "constSizeB4",   "Constituyentes jet4",   50, 10.0,  59.0);
    c6->cd(5); dibujar(tFSR, tMio, "minConstSize",  "minConstSize",  "Min. constituyentes",   50, 10.0,  59.0);
    c6->cd(6); dibujar(tFSR, tMio, "jetNObjects",   "jetNObjects",   "Objetos EFlow totales", 50,  7.0, 145.0);
    c6->cd(7); dibujar(tFSR, tMio, "minJetNObjects","minJetNObjects","Min. objetos EFlow",    50,  1.0,  32.0);
    c6->cd(8); dibujar(tFSR, tMio, "nParticles",    "nParticles",    "N particulas evento",   60,  7.0, 435.0);
    c6->SaveAs("/home/abri/Documentos/Facu/Tesis/comp_constituyentes.pdf");

    cout << endl;
    cout << "========================================" << endl;
    cout << "PDFs guardados en ~/Documentos/Facu/Tesis/" << endl;
    cout << "  comp_costheta_pt.pdf       <- cos(theta) y pT jets" << endl;
    cout << "  comp_masas_jets.pdf        <- masa jets, sumPt, MET" << endl;
    cout << "  comp_masas_pares.pdf       <- masas invariantes pares" << endl;
    cout << "  comp_eventshape.pdf        <- aplanarity, sphericity, thrust" << endl;
    cout << "  comp_boost_asimetrias.pdf  <- boost, asimetrias, deltaR" << endl;
    cout << "  comp_constituyentes.pdf    <- constituyentes de jets" << endl;
    cout << "NOTA: exclYmerge no se grafica porque tu archivo de prueba" << endl;
    cout << "      no tiene esos valores. Funcionara con el archivo real." << endl;
    cout << "========================================" << endl;

    fFSR->Close();
    fMio->Close();
}
