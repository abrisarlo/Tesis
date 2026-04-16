"""
Comparacion3.py
Compara distribuciones de variables entre tres outputs:
  - Abri    : variables_HHbbbb_380GeV.root         (TTree: Variables)
  - Umar    : GammaGammaHHESpreadAllSiD2024XCC_presel_4j3b.npz
  - Santiago: outputTreeSHHbbbbESpreadDurham1034BSplitSampleN.root (TTree: TreeS)

Umar usa nombres de variables distintos. El mapeo a los nombres canonicos
del header_names se define en UMAR_MAP abajo.

Solo se grafican las variables que están presentes en AL MENOS DOS de los tres archivos.
Si una variable falta en alguno de los tres, ese dataset simplemente no se incluye en ese plot.
"""

import uproot
import numpy as np
import matplotlib
matplotlib.use("Agg")          # sin pantalla; cambiar a "TkAgg" o quitar si se corre interactivo
import matplotlib.pyplot as plt
import os

# ─────────────────────────────────────────────
# RUTAS DE LOS ARCHIVOS  ← ajustar si es necesario
# ─────────────────────────────────────────────
FILE_ABRI     = "variables_HHbbbb_380GeV.root"
FILE_UMAR     = "GammaGammaHHESpreadAllSiD2024XCC_presel_4j3b.npz"
FILE_SANTIAGO = "outputTreeSHHbbbbESpreadDurham1034BSplitSampleN.root"

OUTPUT_DIR = "plots_Comparacion3"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ─────────────────────────────────────────────
# MAPEO DE VARIABLES DE UMAR
# Formato: nombre_canonico (header_names) -> nombre_en_npz_de_Umar
# Solo se listan las que tienen equivalente directo confirmado por estadísticas.
# ─────────────────────────────────────────────
UMAR_MAP = {
    "invMassB1"   : "h1_mass",   # mean ~116.5 vs 116.4  ✓
    "invMassB2"   : "h2_mass",   # mean ~103.0 vs 102.9  ✓
    "missingET"   : "met_pt",    # mean ~17.3  vs  17.2  ✓
    "invMass4Jets": "m_4j",      # mismo concepto, distinto proceso (se grafica igual)
}

# ─────────────────────────────────────────────
# VARIABLES A COMPARAR
# ─────────────────────────────────────────────
header_names = [
    "aplanarity",
    "invMassB1",
    "invMassB2",
    "minJetM",
    "sphericity",
    "cosThetaB1",
    "cosThetaB2",
    "cosThetaB3",
    "cosThetaB4",
    "sumPt",
    "jetB1Pt",
    "jetB2Pt",
    "jetB3Pt",
    "jetB4Pt",
    "jetB1M",
    "jetB2M",
    "jetB3M",
    "jetB4M",
    "jetNObjects",
    "minJetNObjects",
    "invMassB1AntiKt",
    "invMassB2AntiKt",
    "nJetsAntiKt",
    "invMassB11Best",
    "invMassB21Best",
    "invMassB12Best",
    "invMassB22Best",
    "invMassB13Best",
    "invMassB23Best",
    "invMassB14Best",
    "invMassB24Best",
    "invMassB15Best",
    "invMassB25Best",
    "invMassB16Best",
    "invMassB26Best",
    "invMassB17Best",
    "invMassB27Best",
    "invMassB18Best",
    "invMassB28Best",
    "distanceZ1MinChiSquaredZZMass",
    "distanceZ2MinChiSquaredZZMass",
    "exclYmerge12",
    "exclYmerge23",
    "exclYmerge34",
    "exclYmerge45",
    "exclYmerge56",
    "invMassZZ1",
    "invMassZZ2",
    "thrust",
    "boostB1",
    "boostB2",
    "boostB3",
    "boostB4",
    "boostSystem",
    "missingET",
    "invMass4Jets",
    "deltaRJetPairs",
    "pTAssymetry",
    "jetPairMassDelta",
    "invMassB1FitBest",
    "invMassB2FitBest",
    "chi2ndfBest",
    "jetPairMassDeltaFit",
]

# ─────────────────────────────────────────────
# CARGA DE DATOS
# ─────────────────────────────────────────────

# ── Abri (ROOT) ──────────────────────────────
print("Cargando Abri...")
with uproot.open(FILE_ABRI) as f_abri:
    tree_abri = f_abri["Variables"]
    abri_data = {k: tree_abri[k].array(library="np") for k in tree_abri.keys()}

# ── Umar (NPZ) ───────────────────────────────
print("Cargando Umar...")
umar_raw  = np.load(FILE_UMAR, allow_pickle=True)
# Renombrar al nombre canónico usando UMAR_MAP
umar_data = {}
for canonical, umar_key in UMAR_MAP.items():
    if umar_key in umar_raw.files:
        umar_data[canonical] = umar_raw[umar_key]

# ── Santiago (ROOT) ──────────────────────────
print("Cargando Santiago...")
with uproot.open(FILE_SANTIAGO) as f_santi:
    tree_santi = f_santi["TreeS"]
    santi_data = {k: tree_santi[k].array(library="np") for k in tree_santi.keys()}

print(f"  Eventos Abri    : {len(list(abri_data.values())[0]):,}")
print(f"  Eventos Umar    : {len(list(umar_data.values())[0]):,}")
print(f"  Eventos Santiago: {len(list(santi_data.values())[0]):,}")
print()

datasets = {
    "Abri"    : abri_data,
    "Umar"    : umar_data,
    "Santiago": santi_data,
}

colors = {
    "Abri"    : "#1f77b4",   # azul
    "Umar"    : "#ff7f0e",   # naranja
    "Santiago": "#2ca02c",   # verde
}

# ─────────────────────────────────────────────
# GRAFICADO
# ─────────────────────────────────────────────

skipped = []
plotted = []

for var in header_names:
    available = {name: ds[var] for name, ds in datasets.items() if var in ds}

    if len(available) < 2:
        skipped.append((var, list(available.keys())))
        continue

    # ------ figura ------
    fig, ax = plt.subplots(figsize=(7, 5))

    all_values = np.concatenate([v for v in available.values()])
    all_values = all_values[np.isfinite(all_values)]
    if len(all_values) == 0:
        skipped.append((var, list(available.keys())))
        plt.close(fig)
        continue

    lo, hi = np.percentile(all_values, 0.5), np.percentile(all_values, 99.5)
    if lo == hi:
        lo -= 0.5
        hi += 0.5
    bins = np.linspace(lo, hi, 60)

    for name, values in available.items():
        vals = values[np.isfinite(values)]
        vals = vals[(vals >= lo) & (vals <= hi)]
        counts, edges = np.histogram(vals, bins=bins)
        # normalizar a densidad
        norm = counts / (counts.sum() * np.diff(edges))
        centers = 0.5 * (edges[:-1] + edges[1:])
        ax.step(centers, norm, where="mid",
                color=colors[name], label=f"{name}  (N={len(vals):,})",
                linewidth=1.5)

    # Nota si la variable de Umar es un renombre
    if var in UMAR_MAP and "Umar" in available:
        ax.set_title(f"{var}\n(Umar: {UMAR_MAP[var]})", fontsize=11, fontweight="bold")
    else:
        ax.set_title(var, fontsize=12, fontweight="bold")

    ax.set_xlabel(var, fontsize=11)
    ax.set_ylabel("Densidad normalizada", fontsize=11)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()

    out_path = os.path.join(OUTPUT_DIR, f"{var}.png")
    fig.savefig(out_path, dpi=120)
    plt.close(fig)
    plotted.append(var)
    print(f"  ✓ {var}" + (f"  [Umar={UMAR_MAP[var]}]" if var in UMAR_MAP and "Umar" in available else ""))

# ─────────────────────────────────────────────
# RESUMEN
# ─────────────────────────────────────────────
print(f"\n{'='*60}")
print(f"  Plots generados : {len(plotted)}  → carpeta '{OUTPUT_DIR}/'")
print(f"  Variables saltadas (< 2 datasets): {len(skipped)}")
for v, who in skipped:
    print(f"    · {v}  (solo en: {who if who else 'ninguno'})")
print("="*60)
