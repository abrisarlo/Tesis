"""
comparar_outputs.py
Compara el output del codigo de Umar (.npz) con el de Abri (.root)
Genera plots lado a lado de las variables en comun.

Uso:
    python3 comparar_outputs.py

Requiere: uproot, numpy, matplotlib, awkward
"""

import numpy as np
import uproot
import matplotlib.pyplot as plt
import os

# ─── Rutas (cambialas si hace falta) ────────────────────────────────────────
NPZ_FILE  = "GammaGammaHHESpreadAllSiD2024XCC_presel_4j3b.npz"
ROOT_FILE = "variables_HHbbbb_380GeV.root"
OUT_DIR   = "plots_comparacion"
os.makedirs(OUT_DIR, exist_ok=True)

# ─── Cargar datos ────────────────────────────────────────────────────────────
print("Cargando NPZ (Umar)...")
npz = np.load(NPZ_FILE, allow_pickle=True)

print("Cargando ROOT (Abri)...")
root_tree = uproot.open(ROOT_FILE)["Variables"]

print(f"  Eventos Umar (.npz) : {len(npz['jet1_pt']):,}")
print(f"  Eventos Abri    (.root) : {root_tree.num_entries:,}")
print()

# ─── Mapeo de variables en comun ─────────────────────────────────────────────
# Formato: (nombre_Umar_npz, nombre_Abri_root, label_eje_x, n_bins, (xmin, xmax))
# CAMBIOS respecto a la version anterior:
#   - "m_j1j2" -> "h1_mass" para comparar con invMassB11Best (ambas son el mejor par por chi2)
#   - "h1_mass" e "h2_mass" comparadas con invMassB1/invMassB2 (equivalentes)
VARIABLES_COMUNES = [
    ("jet1_pt",  "ptJet1",        "pT jet1 [GeV]",             50, (0,  180)),
    ("jet2_pt",  "ptJet2",        "pT jet2 [GeV]",             50, (0,  160)),
    ("jet3_pt",  "ptJet3",        "pT jet3 [GeV]",             50, (0,  100)),
    ("jet4_pt",  "ptJet4",        "pT jet4 [GeV]",             50, (0,   90)),
    ("jet1_m",   "massJet1",      "Masa jet1 [GeV]",           50, (0,   35)),
    ("jet2_m",   "massJet2",      "Masa jet2 [GeV]",           50, (0,   30)),
    ("jet3_m",   "massJet3",      "Masa jet3 [GeV]",           50, (0,   25)),
    ("jet4_m",   "massJet4",      "Masa jet4 [GeV]",           50, (0,   20)),
    ("h1_mass",  "invMassB1",     "Masa inv. H1 [GeV]",        50, (0,  250)),
    ("h2_mass",  "invMassB2",     "Masa inv. H2 [GeV]",        50, (0,  250)),
    ("h1_mass",  "invMassB11Best","Masa par mejor chi2 [GeV]", 50, (0,  300)),
    ("m_4j",     "invMass4Jets",  "Masa inv. 4 jets [GeV]",    50, (0,  500)),
    ("met_pt",   "missingET",     "Missing ET [GeV]",          50, (0,  160)),
    ("jet1_eta", "cosThetaJet1",  None, None, None),  # no comparables directamente
]

# Filtrar solo las que son directamente comparables (tienen label)
comparables = [(s, a, lbl, nb, rng) for s, a, lbl, nb, rng in VARIABLES_COMUNES
               if lbl is not None]

# ─── Funcion para hacer un plot de comparacion ───────────────────────────────
def plot_comparison(var_umar, var_abri, xlabel, nbins, xrange, filename):
    data_s = npz[var_umar].astype(float)
    data_a = root_tree[var_abri].array(library="np").astype(float)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    fig.suptitle(f"Comparacion: {var_umar} vs {var_abri}", fontsize=13, fontweight="bold")

    bins = np.linspace(xrange[0], xrange[1], nbins + 1)

    # ── Histogramas normalizados superpuestos ──
    ax = axes[0]
    hs, _ = np.histogram(data_s, bins=bins)
    ha, _ = np.histogram(data_a, bins=bins)
    hs = hs / hs.sum() if hs.sum() > 0 else hs
    ha = ha / ha.sum() if ha.sum() > 0 else ha

    ax.step(bins[:-1], hs, where="post", color="steelblue", label=f"Umar  (n={len(data_s):,})", linewidth=1.5)
    ax.step(bins[:-1], ha, where="post", color="crimson",   label=f"Abri      (n={len(data_a):,})", linewidth=1.5, linestyle="--")
    ax.set_xlabel(xlabel, fontsize=11)
    ax.set_ylabel("Fraccion de eventos", fontsize=11)
    ax.legend(fontsize=10)
    ax.set_title("Superpuestos (normalizados)")

    # ── Ratio Umar / Abri ──
    ax2 = axes[1]
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(ha > 0, hs / ha, np.nan)
    bin_centers = 0.5 * (bins[:-1] + bins[1:])
    ax2.plot(bin_centers, ratio, "o-", color="darkorange", markersize=3, linewidth=1.2)
    ax2.axhline(1.0, color="gray", linestyle="--", linewidth=1)
    ax2.set_xlabel(xlabel, fontsize=11)
    ax2.set_ylabel("Ratio Umar / Abri", fontsize=11)
    ax2.set_ylim(0, 3)
    ax2.set_title("Ratio bin a bin")

    plt.tight_layout()
    plt.savefig(filename, dpi=120, bbox_inches="tight")
    plt.close()
    print(f"  Guardado: {filename}")


# ─── Generar todos los plots ─────────────────────────────────────────────────
print("Generando plots...")
for var_s, var_a, lbl, nb, rng in comparables:
    fname = os.path.join(OUT_DIR, f"comp_{var_s}_vs_{var_a}.png")
    plot_comparison(var_s, var_a, lbl, nb, rng, fname)

# ─── Plot resumen: todos los pT en una figura ────────────────────────────────
print("\nGenerando resumen de pT jets...")
fig, axes = plt.subplots(2, 2, figsize=(12, 9))
fig.suptitle("pT de los 4 jets: Umar vs Abri", fontsize=14, fontweight="bold")

for idx, (var_s, var_a, lbl, nb, rng) in enumerate([
    ("jet1_pt","ptJet1","pT jet1 [GeV]",50,(0,180)),
    ("jet2_pt","ptJet2","pT jet2 [GeV]",50,(0,160)),
    ("jet3_pt","ptJet3","pT jet3 [GeV]",50,(0,100)),
    ("jet4_pt","ptJet4","pT jet4 [GeV]",50,(0, 90)),
]):
    ax = axes[idx // 2][idx % 2]
    bins = np.linspace(rng[0], rng[1], nb + 1)
    ds = npz[var_s].astype(float)
    da = root_tree[var_a].array(library="np").astype(float)
    hs, _ = np.histogram(ds, bins=bins)
    ha, _ = np.histogram(da, bins=bins)
    hs = hs / hs.sum(); ha = ha / ha.sum()
    ax.step(bins[:-1], hs, where="post", color="steelblue", label="Umar", lw=1.5)
    ax.step(bins[:-1], ha, where="post", color="crimson",   label="Abri",     lw=1.5, ls="--")
    ax.set_xlabel(lbl); ax.set_ylabel("Fraccion"); ax.legend()

plt.tight_layout()
resumen_path = os.path.join(OUT_DIR, "resumen_pT_jets.png")
plt.savefig(resumen_path, dpi=130, bbox_inches="tight")
plt.close()
print(f"  Guardado: {resumen_path}")

# ─── Estadisticas descriptivas ───────────────────────────────────────────────
print("\n" + "="*70)
print(f"{'Variable':<25} {'Media_Santi':>12} {'Media_Abri':>12} {'Std_Santi':>12} {'Std_Abri':>12}")
print("="*70)
for var_s, var_a, lbl, nb, rng in comparables:
    ds = npz[var_s].astype(float)
    da = root_tree[var_a].array(library="np").astype(float)
    print(f"{var_s+' vs '+var_a:<25} {ds.mean():12.3f} {da.mean():12.3f} {ds.std():12.3f} {da.std():12.3f}")
print("="*70)
print("\nListo! Plots guardados en:", OUT_DIR)
