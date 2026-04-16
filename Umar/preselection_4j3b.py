#!/usr/bin/env python3
"""
Preselection script for Higgs to bb analysis:
- Require 4 jets, 3 of them b-tagged
- Extract jet kinematics, dijet combinations, and Higgs candidates
"""

import os
import glob
import argparse
import numpy as np
import uproot
import awkward as ak
from collections import OrderedDict
import sys
from pathlib import Path


class Tee:
    def __init__(self, filename, mode="w"):
        self.file = open(filename, mode, buffering=1)
    
    def write(self, data):
        sys.__stdout__.write(data)
        self.file.write(data)
    
    def flush(self):
        sys.__stdout__.flush()
        self.file.flush()


# File configuration - total generated events (before preselection)
# FILE_CONFIG = {
#     "eGammaqqHXAllSiD2024XCC1.root": 3282.34,
#     "eGammaqqHXAllSiD2024XCC2.root": 3282.34,
#     "GammaGammaZHESpreadAllSiD2024XCC.root": 8202.0,
#     "GammaGammaZZESpreadAllSiD2024XCC.root": 1378000.0,
#     "eGammaqqqqXAllSiD2024XCC.root": 7681.04,
#     "eGammaqqXAllSiD2024XCC.root": 41195.17,
#     "pebbESpreadAllSiD2024XCC.root": 753_615.07,
#     "GammaGammabbbbqqESpreadAllSiD2024XCC.root": 7681.04,
#     "pebbqqESpreadAllSiD2024XCC.root": 3282.34,
#     "GammaGammaHHESpreadAllSiD2024XCC.root": 3317.54386,
#     "peqqHESpreadAllSiD2024XCC.root": 753615.07,
#     "GammaGammattAllSiD2024XCC.root": 2866000.0,
#     "pettESpreadAllSiD2024XCC.root": 57001.14,
#     "GammaGammaWWESpreadAllSiD2024XCC.root": 3813000.0,
# }

# # 10-year yield mapping (cross-sections or total yields for 10 years of running)56
# TEN_YR_YIELDS = {
#     "eGammaqqHXAllSiD2024XCC1.root": 3282.34,
#     "eGammaqqHXAllSiD2024XCC2.root": 3282.34,
#     "GammaGammaZHESpreadAllSiD2024XCC.root": 8202.0,
#     "GammaGammaZZESpreadAllSiD2024XCC.root": 1378000.0,
#     "eGammaqqqqXAllSiD2024XCC.root": 7681.04,
#     "eGammaqqXAllSiD2024XCC.root": 41195.17,
#     "pebbESpreadAllSiD2024XCC.root": 753_615.07,
#     "GammaGammabbbbqqESpreadAllSiD2024XCC.root": 7681.04,
#     "pebbqqESpreadAllSiD2024XCC.root": 3282.34,
#     "GammaGammaHHESpreadAllSiD2024XCC.root": 3317.54386,
#     "peqqHESpreadAllSiD2024XCC.root": 753615.07,
#     "GammaGammattAllSiD2024XCC.root": 2866000.0,
#     "pettESpreadAllSiD2024XCC.root": 57001.14,
#     "GammaGammaWWESpreadAllSiD2024XCC.root": 3813000.0,
# }

FILE_CONFIG = {
    "GammaGammaHHESpreadAllSiD2024XCC.root": 3317.54386,
}

TEN_YR_YIELDS = {
    "GammaGammaHHESpreadAllSiD2024XCC.root": 3317.54386,
}

def delta_phi(phi1, phi2):
    """Compute delta phi with proper wrapping"""
    dphi = phi1 - phi2
    dphi = (dphi + np.pi) % (2 * np.pi) - np.pi
    return np.abs(dphi)


def delta_eta(eta1, eta2):
    """Compute delta eta"""
    return eta1 - eta2


def delta_r(eta1, eta2, phi1, phi2):
    """Compute delta R"""
    deta = delta_eta(eta1, eta2)
    dphi = delta_phi(phi1, phi2)
    return np.sqrt(deta**2 + dphi**2)


def invariant_mass(pt1, eta1, phi1, m1, pt2, eta2, phi2, m2):
    """Compute invariant mass of two particles"""
    px1 = pt1 * np.cos(phi1)
    py1 = pt1 * np.sin(phi1)
    pz1 = pt1 * np.sinh(eta1)
    
    px2 = pt2 * np.cos(phi2)
    py2 = pt2 * np.sin(phi2)
    pz2 = pt2 * np.sinh(eta2)
    
    # Energy
    p1_sq = px1**2 + py1**2 + pz1**2
    E1 = np.sqrt(np.maximum(p1_sq + m1**2, 0.0))
    
    p2_sq = px2**2 + py2**2 + pz2**2
    E2 = np.sqrt(np.maximum(p2_sq + m2**2, 0.0))
    
    # Combined 4-vector
    Px = px1 + px2
    Py = py1 + py2
    Pz = pz1 + pz2
    E = E1 + E2
    
    m_inv_sq = E**2 - (Px**2 + Py**2 + Pz**2)
    return np.sqrt(np.maximum(m_inv_sq, 0.0))


def durham_distance(pt1, eta1, phi1, pt2, eta2, phi2):
    """Compute Durham distance (jet clustering metric)"""
    dR = delta_r(eta1, eta2, phi1, phi2)
    return np.minimum(pt1, pt2)**2 * dR**2


def print_cutflow_table(total, counts, title=None):
    """Print a formatted cutflow table"""
    names = list(counts.keys())
    print("-" * 80)
    if title:
        print(title)
    print(f"{'Cut':<30} {'Pass':>12}  {'Eff(cum)':>10}  {'Eff(step)':>10}")
    print("-" * 80)
    cum_prev = total
    for name in names:
        cum_pass = counts[name]
        eff_cum = cum_pass / total if total > 0 else 0.0
        eff_step = cum_pass / cum_prev if cum_prev > 0 else 0.0
        print(f"{name:<30} {cum_pass:12,d}  {eff_cum:10.4f}  {eff_step:10.4f}")
        cum_prev = cum_pass
    print("-" * 80)


def compute_jet_features(pt, eta, phi, m, btag, met_pt, met_eta, met_phi):
    """
    Compute all requested features for a 4-jet event with 3 b-tags
    
    Returns a dictionary of features (all as numpy arrays)
    """
    features = OrderedDict()
    
    # Individual jet kinematics
    for i in range(4):
        features[f"jet{i+1}_pt"] = pt[:, i]
        features[f"jet{i+1}_eta"] = eta[:, i]
        features[f"jet{i+1}_phi"] = phi[:, i]
        features[f"jet{i+1}_m"] = m[:, i]
    
    # All pairwise combinations (6 total)
    pairs = [(0,1), (0,2), (0,3), (1,2), (1,3), (2,3)]
    
    # Delta R for all pairs
    for i, j in pairs:
        dr = delta_r(eta[:, i], eta[:, j], phi[:, i], phi[:, j])
        features[f"dR_j{i+1}j{j+1}"] = dr
    
    # Dijet invariant mass for all pairs
    for i, j in pairs:
        m_ij = invariant_mass(pt[:, i], eta[:, i], phi[:, i], m[:, i],
                             pt[:, j], eta[:, j], phi[:, j], m[:, j])
        features[f"m_j{i+1}j{j+1}"] = m_ij
    
    # Durham distance for all pairs
    for i, j in pairs:
        d_ij = durham_distance(pt[:, i], eta[:, i], phi[:, i],
                              pt[:, j], eta[:, j], phi[:, j])
        features[f"durham_j{i+1}j{j+1}"] = d_ij
    
    # Triple combinations (4 total)
    triples = [(0,1,2), (0,1,3), (0,2,3), (1,2,3)]
    
    for i, j, k in triples:
        m_ijk = invariant_mass(
            pt[:, i] + pt[:, j], 
            eta[:, i], phi[:, i],  # rough approximation for combined eta/phi
            m[:, i] + m[:, j],
            pt[:, k], eta[:, k], phi[:, k], m[:, k]
        )
        features[f"m_j{i+1}j{j+1}j{k+1}"] = m_ijk
    
    # All 4 jets mass
    m_4j = invariant_mass(
        pt[:, 0] + pt[:, 1],
        eta[:, 0], phi[:, 0],
        m[:, 0] + m[:, 1],
        pt[:, 2] + pt[:, 3],
        eta[:, 2], phi[:, 2],
        m[:, 2] + m[:, 3]
    )                               #Esta bien esto que hace????
    features["m_4j"] = m_4j
    
    # Higgs candidates: minimize chi2 = (m_ij - m_H)^2/sigma_m^2 + (m_kl - m_H)^2/sigma_m^2
    M_H = 125.0
    SIGMA_M = 10.0
    
    # Three possible pairings
    m12 = features["m_j1j2"]
    m34 = features["m_j3j4"]
    m13 = features["m_j1j3"]
    m24 = features["m_j2j4"]
    m14 = features["m_j1j4"]
    m23 = features["m_j2j3"]
    
    chi2_12_34 = ((m12 - M_H)**2 + (m34 - M_H)**2) / (SIGMA_M**2)
    chi2_13_24 = ((m13 - M_H)**2 + (m24 - M_H)**2) / (SIGMA_M**2)
    chi2_14_23 = ((m14 - M_H)**2 + (m23 - M_H)**2) / (SIGMA_M**2)
    
    # Find minimum chi2
    chi2_min = np.minimum(np.minimum(chi2_12_34, chi2_13_24), chi2_14_23)
    
    # Assign Higgs masses based on minimum chi2
    h1_mass = np.where(chi2_12_34 == chi2_min, m12,
                       np.where(chi2_13_24 == chi2_min, m13, m14))
    h2_mass = np.where(chi2_12_34 == chi2_min, m34,
                       np.where(chi2_13_24 == chi2_min, m24, m23))
    
    # Determine which jets form the Higgs candidates
    h1_j1 = np.where(chi2_12_34 == chi2_min, 0, np.where(chi2_13_24 == chi2_min, 0, 0))
    h1_j2 = np.where(chi2_12_34 == chi2_min, 1, np.where(chi2_13_24 == chi2_min, 2, 3))
    
    # Higgs 2 kinematics (combined from remaining jets)
    h2_j1 = np.where(chi2_12_34 == chi2_min, 2, np.where(chi2_13_24 == chi2_min, 1, 1))
    h2_j2 = np.where(chi2_12_34 == chi2_min, 3, np.where(chi2_13_24 == chi2_min, 3, 2))
    
    # Proper 2D indexing with row indices
    n_evt = len(pt)
    row_idx = np.arange(n_evt)
    
    # Higgs 1 kinematics (combined)
    h1_pt = pt[row_idx, h1_j1] + pt[row_idx, h1_j2]
    h1_eta = (eta[row_idx, h1_j1] + eta[row_idx, h1_j2]) / 2.0
    h1_phi = (phi[row_idx, h1_j1] + phi[row_idx, h1_j2]) / 2.0
    
    # Higgs 2 kinematics (combined)
    h2_pt = pt[row_idx, h2_j1] + pt[row_idx, h2_j2]
    h2_eta = (eta[row_idx, h2_j1] + eta[row_idx, h2_j2]) / 2.0
    h2_phi = (phi[row_idx, h2_j1] + phi[row_idx, h2_j2]) / 2.0
    
    features["h1_mass"] = h1_mass
    features["h1_pt"] = h1_pt
    features["h1_eta"] = h1_eta
    features["h1_phi"] = h1_phi
    
    features["h2_mass"] = h2_mass
    features["h2_pt"] = h2_pt
    features["h2_eta"] = h2_eta
    features["h2_phi"] = h2_phi
    
    features["chi2_min"] = chi2_min
    
    # MET kinematics
    features["met_pt"] = met_pt
    features["met_eta"] = met_eta
    features["met_phi"] = met_phi
    
    return features


def main():
    ap = argparse.ArgumentParser(description="Preselection for Higgs->bb: 4 jets, 3 b-tagged")
    ap.add_argument("indir", help="Input directory with ROOT files")
    ap.add_argument("outdir", help="Output directory")
    ap.add_argument("--tree", default="Delphes", help="ROOT tree name")
    ap.add_argument("--step", type=int, default=50000, help="Chunk size")
    
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    
    # Setup logging
    log_path = os.path.join(args.outdir, "preselection_log.txt")
    tee = Tee(log_path, "w")
    sys.stdout = tee
    sys.stderr = tee
    
    print("=" * 80)
    print("PRESELECTION: 4 JETS, 3 B-TAGGED")
    print("=" * 80)
    print(f"Input directory: {args.indir}")
    print(f"Output directory: {args.outdir}")
    print()
    
    # Find ROOT files
    root_files = sorted(glob.glob(os.path.join(args.indir, "*.root")))
    print(f"Found {len(root_files)} ROOT files")
    print()
    
    # Branches to read
    branches = [
        "Jet10.PT", "Jet10.Eta", "Jet10.Phi", "Jet10.Mass", "Jet10.BTag",
        "MissingET.MET", "MissingET.Eta", "MissingET.Phi", "Electron.PT", "Muon.PT"
    ]
    
    grand_total = 0
    grand_pass = 0
    grand_4jets = 0
    grand_3btags = 0
    
    # Dictionary to track per-file statistics
    file_stats = OrderedDict()
    
    # Process each file
    for root_file in root_files:
        filename = os.path.basename(root_file)
        
        if filename not in FILE_CONFIG:
            print(f"[SKIP] {filename} - not in FILE_CONFIG")
            continue
        
        print(f"\n{'='*80}")
        print(f"Processing: {filename}")
        print(f"Expected events: {FILE_CONFIG[filename]:,.0f}")
        print(f"{'='*80}")
        
        file_total = 0
        file_4jets = 0
        file_pass = 0
        all_features = OrderedDict()
        
        try:
            # Iterate through chunks
            for chunk_idx, arr in enumerate(uproot.iterate(
                f"{root_file}:{args.tree}",
                branches,
                library="ak",
                step_size=args.step
            )):

                
                 #Esto lo agregue yo (Abri)
                ##################################
        
                print("\n====================")
                print("CHUNK", chunk_idx)
                print("Eventos en chunk:", len(arr["Jet10.PT"]))
                
                # debug evento 0 crudo
                print("\n--- EVENTO 0 RAW ---")
                print("Jets pt:", arr["Jet10.PT"][0])
                print("Electrones:", arr["Electron.PT"][0])
                print("Muones:", arr["Muon.PT"][0])
                ##################################
                
                n_events = len(arr["Jet10.PT"])
                file_total += n_events
                grand_total += n_events
                
                # Get jet info
                j_pt = arr["Jet10.PT"]
                j_eta = arr["Jet10.Eta"]
                j_phi = arr["Jet10.Phi"]
                j_m = arr["Jet10.Mass"]
                j_btag = arr["Jet10.BTag"]
                
                # Get MET info
                met_pt = ak.to_numpy(arr["MissingET.MET"][:, 0])
                met_eta = ak.to_numpy(arr["MissingET.Eta"][:, 0])
                met_phi = ak.to_numpy(arr["MissingET.Phi"][:, 0])
                
                # Preselection cuts
                n_jets = ak.num(j_pt, axis=1)
                n_btags = ak.sum(j_btag > 0.5, axis=1)  # BTag > 0.5 is typically b-tagged
                
                #Agregado
                ###################################3
                print("\n--- DEBUG CUTS ---")
                print("Ejemplo n_jets:", n_jets[:5])
                print("Ejemplo n_btags:", n_btags[:5])
                #############################
                
                # Cut 1: >= 4 jets
                has_4jets = (n_jets >= 4)
                n_4jets_chunk = ak.sum(has_4jets)
                file_4jets += n_4jets_chunk
                grand_4jets += n_4jets_chunk
                
                # Cut 2: >= 3 b-tagged
                #presel_mask = has_4jets & (n_btags >= 3)
                presel_mask = has_4jets                        #HICE EL CAMBIO ACA
                
                n_pass_chunk = ak.sum(presel_mask)
                print("Eventos que pasan todo (chunk):", n_pass_chunk)
                file_pass += n_pass_chunk
                grand_pass += n_pass_chunk
                
                # Keep only passing events and select leading 4 jets
                if n_pass_chunk > 0:
                    j_pt_sel = j_pt[presel_mask][:, :4]
                    j_eta_sel = j_eta[presel_mask][:, :4]
                    j_phi_sel = j_phi[presel_mask][:, :4]
                    j_m_sel = j_m[presel_mask][:, :4]
                    
                    met_pt_sel = met_pt[ak.to_numpy(presel_mask)]
                    met_eta_sel = met_eta[ak.to_numpy(presel_mask)]
                    met_phi_sel = met_phi[ak.to_numpy(presel_mask)]

                    #Agregdo
                    ##########################
                    print("\n--- DEBUG MET ---")
                    for i in range(min(3, len(met_pt_sel))):
                        print("EVENTO", i, "MET:", met_pt_sel[i])
                    ###############################################3
                    
                    # Convert to numpy
                    j_pt_np = ak.to_numpy(j_pt_sel)
                    j_eta_np = ak.to_numpy(j_eta_sel)
                    j_phi_np = ak.to_numpy(j_phi_sel)
                    j_m_np = ak.to_numpy(j_m_sel)

                    #Agregado
                    ######################################
                    print("\n--- DEBUG JETS SELECCIONADOS ---")
                    
                    for i in range(min(3, len(j_pt_np))):
                        print("\nEVENTO", i)
                        print("pt:", j_pt_np[i])
                        print("eta:", j_eta_np[i])
                        print("phi:", j_phi_np[i])
                        print("m:", j_m_np[i])
                    #######################################
                    
                    # Compute features
                    chunk_features = compute_jet_features(
                        j_pt_np, j_eta_np, j_phi_np, j_m_np,
                        None,  # btag array (not used in feature computation)
                        met_pt_sel, met_eta_sel, met_phi_sel
                    )
                    #Agregado
                    #########################
                    print("\n--- DEBUG MASAS ---")

                    for i in range(min(3, len(j_pt_np))):
                        print("\nEVENTO", i)
                        print("m_j1j2:", chunk_features["m_j1j2"][i])
                        print("m_j1j3:", chunk_features["m_j1j3"][i])
                        print("m_j1j4:", chunk_features["m_j1j4"][i])
                        print("m_j2j3:", chunk_features["m_j2j3"][i])
                        print("m_j2j4:", chunk_features["m_j2j4"][i])
                        print("m_j3j4:", chunk_features["m_j3j4"][i])
                    ###########################################################
                    
                    
                    # Store features
                    for name, values in chunk_features.items():
                        if name not in all_features:
                            all_features[name] = []
                        all_features[name].append(values)
                
                if (chunk_idx + 1) % 5 == 0 or chunk_idx == 0:
                    print(f"  Chunk {chunk_idx+1}: Processed {n_events:,} events, "
                          f"passed {n_pass_chunk:,}")
        
        except Exception as e:
            print(f"[ERROR] Failed to process {filename}: {e}")
            continue
        
        # Store statistics
        file_stats[filename] = {
            "total": file_total,
            "n_4jets": file_4jets,
            "n_pass": file_pass,
            "expected_events": FILE_CONFIG.get(filename, 0)
        }
        
        # Concatenate all chunks for this file
        if all_features:
            for name in all_features:
                all_features[name] = np.concatenate(all_features[name], axis=0)
            
            # Save to NPZ
            safe_name = filename.replace(".root", "")
            out_npz = os.path.join(args.outdir, f"{safe_name}_presel_4j3b.npz")
            
            # Prepare dict for saving
            save_dict = dict(all_features)
            save_dict["feature_names"] = np.array(list(all_features.keys()))
            
            np.savez_compressed(out_npz, **save_dict)
            print(f"  Saved: {out_npz}")
        
        frac = (file_pass / file_total) if file_total > 0 else 0.0
        print(f"\n[SUMMARY] {filename}")
        print(f"  Total events:  {file_total:,}")
        print(f"  Pass presel:   {file_pass:,}")
        print(f"  Efficiency:    {frac:.4f}")
    
    # Grand summary with detailed cutflow
    print(f"\n{'='*80}")
    print(f"GRAND SUMMARY - CUTFLOW TABLE")
    print(f"{'='*80}")
    
    # Create cutflow dictionary
    cutflow = OrderedDict()
    cutflow["All events"] = grand_total
    cutflow["N_jets >= 4"] = grand_4jets
    cutflow["N_btags >= 3"] = grand_pass
    
    print_cutflow_table(grand_total, cutflow, title="PRESELECTION CUTFLOW (ALL FILES)")
    
    # Per-file summary table with 10-year yields
    print(f"\n{'='*80}")
    print("PER-FILE SUMMARY WITH 10-YEAR SCALED YIELDS")
    print(f"{'='*80}")
    print(f"{'Filename':<35} {'Total':>10}  {'Pass':>10}  {'Eff':>8}  {'10-Yr Yield':>14}")
    print("-" * 90)
    
    total_10yr_unscaled = 0.0
    total_10yr_scaled = 0.0
    
    for fname, stats in file_stats.items():
        total = stats["total"]
        n_pass = stats["n_pass"]
        eff = (n_pass / total) if total > 0 else 0.0
        
        # Get 10-year yield
        yr10_yield = TEN_YR_YIELDS.get(fname, 0.0)
        # Scale by efficiency
        yr10_scaled = yr10_yield * eff
        
        total_10yr_unscaled += yr10_yield
        total_10yr_scaled += yr10_scaled
        
        print(f"{fname:<35} {total:10,d}  {n_pass:10,d}  {eff:8.4f}  {yr10_scaled:14,.1f}")
    
    print("-" * 90)
    grand_eff = (grand_pass / grand_total) if grand_total > 0 else 0.0
    print(f"{'TOTAL':<35} {grand_total:10,d}  {grand_pass:10,d}  {grand_eff:8.4f}  {total_10yr_scaled:14,.1f}")
    print(f"{'='*90}")
    
    # Summary statistics
    print(f"\nSUMMARY STATISTICS:")
    print(f"  Total events processed:        {grand_total:,}")
    print(f"  Events with >= 4 jets:        {grand_4jets:,}  ({grand_4jets/grand_total*100:.2f}%)")
    print(f"  Events passing full presel:   {grand_pass:,}  ({grand_pass/grand_total*100:.2f}%)")
    print(f"  ")
    print(f"  10-Year yields (unscaled):    {total_10yr_unscaled:,.1f}")
    print(f"  10-Year yields (post-presel): {total_10yr_scaled:,.1f}")
    print(f"  Preselection efficiency:      {grand_eff:.4f}")
    print(f"{'='*90}")


if __name__ == "__main__":
    main()

