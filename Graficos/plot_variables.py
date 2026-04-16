"""
plot_variables.py
Lee variables_etapa2.root y guarda todas las distribuciones en un solo PDF.

Uso:
    python3 plot_variables.py
"""

import uproot
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import os

# ─── Rutas ───────────────────────────────────────────────────────────────────
ROOT_FILE  = "variables_etapa2.root"
OUTPUT_PDF = "distribuciones_etapa2.pdf"

# ─── Cargar datos ─────────────────────────────────────────────────────────────
print(f"Cargando {ROOT_FILE}...")
tree = uproot.open(ROOT_FILE)["Variables"]
branches = tree.keys()
N = tree.num_entries
print(f"  Variables: {len(branches)}   Eventos: {N:,}")

# ─── Configuracion de plots por variable (xmin, xmax o None para auto) ───────
# None = dejar que matplotlib elija el rango automaticamente
RANGOS = {
    "cosThetaJet1"                   : (-1,    1),
    "cosThetaJet2"                   : (-1,    1),
    "cosThetaJet3"                   : (-1,    1),
    "cosThetaJet4"                   : (-1,    1),
    "ptJet1"                         : (0,   200),
    "ptJet2"                         : (0,   180),
    "ptJet3"                         : (0,   120),
    "ptJet4"                         : (0,   100),
    "sumPt"                          : (0,   500),
    "massJet1"                       : (0,    40),
    "massJet2"                       : (0,    35),
    "massJet3"                       : (0,    30),
    "massJet4"                       : (0,    25),
    "missingET"                      : (0,   100),
    "aplanarity"                     : (0,   0.5),
    "sphericity"                     : (0,     1),
    "thrust"                         : (0.5,   1),
    "invMassB1"                      : (0,   300),
    "invMassB2"                      : (0,   300),
    "invMass4Jets"                   : (0,   500),
    "minJetM"                        : (0,    50),
    "jetPairMassDelta"               : (0, 10000),
    "invMassZZ1"                     : (0,   300),
    "invMassZZ2"                     : (0,   300),
    "distanceZ1MinChiSquaredZZMass"  : (0,   150),
    "distanceZ2MinChiSquaredZZMass"  : (0,   150),
    "boostB1"                        : (-1,    1),
    "boostB2"                        : (-1,    1),
    "boostB3"                        : (-1,    1),
    "boostB4"                        : (-1,    1),
    "boostSystem"                    : (-0.5, 0.5),
    "nEFlowObjects"                  : (0,   300),
    "minJetNObjects"                 : (0,    60),
    "eAssymetry"                     : (0,     1),
    "pTAssymetry"                    : (0,     1),
    "deltaRJetPairs"                 : (0,     7),
    "invMassB11Best"                 : (0,   300),
    "invMassB21Best"                 : (0,   300),
    "invMassB12Best"                 : (0,   300),
    "invMassB22Best"                 : (0,   300),
    "invMassB13Best"                 : (0,   300),
    "invMassB23Best"                 : (0,   300),
    "exclYmerge12"                   : (0,     1),
    "exclYmerge23"                   : (0,   0.3),
    "exclYmerge34"                   : (0,   0.1),
    "exclYmerge45"                   : (0,  0.03),
    "exclYmerge56"                   : (0,  0.01),
}

NBINS = 60
PLOTS_POR_PAGINA = 6   # 2 columnas x 3 filas

# ─── Generar PDF ──────────────────────────────────────────────────────────────
print(f"Generando {OUTPUT_PDF}...")

with PdfPages(OUTPUT_PDF) as pdf:

    # Pagina de portada
    fig_cover = plt.figure(figsize=(11, 8.5))
    fig_cover.text(0.5, 0.6, "Distribuciones - variables_etapa2.root",
                   ha="center", va="center", fontsize=22, fontweight="bold")
    fig_cover.text(0.5, 0.48, f"Variables: {len(branches)}     Eventos: {N:,}",
                   ha="center", va="center", fontsize=14, color="gray")
    fig_cover.text(0.5, 0.38, f"Archivo: {ROOT_FILE}",
                   ha="center", va="center", fontsize=11, color="gray", style="italic")
    pdf.savefig(fig_cover, bbox_inches="tight")
    plt.close(fig_cover)

    # Una pagina por cada grupo de PLOTS_POR_PAGINA variables
    branch_list = list(branches)
    n_cols = 2
    n_rows = PLOTS_POR_PAGINA // n_cols

    for page_start in range(0, len(branch_list), PLOTS_POR_PAGINA):
        chunk = branch_list[page_start : page_start + PLOTS_POR_PAGINA]

        fig, axes = plt.subplots(n_rows, n_cols, figsize=(14, 10))
        axes = axes.flatten()

        for ax_idx, var in enumerate(chunk):
            ax = axes[ax_idx]
            vals = tree[var].array(library="np").astype(float)
            vals = vals[np.isfinite(vals)]

            # Rango
            rng = RANGOS.get(var, None)
            if rng is not None:
                lo, hi = rng
            else:
                lo = np.percentile(vals, 0.5)
                hi = np.percentile(vals, 99.5)
            if lo == hi:
                lo -= 0.5; hi += 0.5

            vals_clipped = vals[(vals >= lo) & (vals <= hi)]
            bins = np.linspace(lo, hi, NBINS + 1)

            counts, edges = np.histogram(vals_clipped, bins=bins)
            norm = counts / (counts.sum() * np.diff(edges)) if counts.sum() > 0 else counts
            centers = 0.5 * (edges[:-1] + edges[1:])

            ax.step(centers, norm, where="mid", color="#1f77b4", linewidth=1.5)
            ax.fill_between(centers, norm, step="mid", alpha=0.25, color="#1f77b4")

            ax.set_title(var, fontsize=10, fontweight="bold")
            ax.set_xlabel(var, fontsize=8)
            ax.set_ylabel("Densidad", fontsize=8)
            ax.tick_params(labelsize=7)
            ax.grid(True, alpha=0.3)

            # Stats en el plot
            stats_txt = f"N={len(vals):,}\nMedia={vals.mean():.2g}\nStd={vals.std():.2g}"
            ax.text(0.97, 0.97, stats_txt, transform=ax.transAxes,
                    fontsize=6.5, va="top", ha="right",
                    bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.7))

        # Ocultar axes sobrantes si el chunk es menor que PLOTS_POR_PAGINA
        for ax_idx in range(len(chunk), len(axes)):
            axes[ax_idx].set_visible(False)

        fig.suptitle(f"variables_etapa2.root  —  pag. {page_start // PLOTS_POR_PAGINA + 2}",
                     fontsize=11, color="gray")
        fig.tight_layout(rect=[0, 0, 1, 0.97])
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)
        print(f"  Pagina {page_start // PLOTS_POR_PAGINA + 2}: {[v for v in chunk]}")

    # Metadata del PDF
    d = pdf.infodict()
    d["Title"]   = "Distribuciones variables_etapa2"
    d["Subject"] = "Variables fisicas para BDT - Etapa 2"

print(f"\nListo! PDF guardado en: {OUTPUT_PDF}")
print(f"Paginas totales: {len(branch_list) // PLOTS_POR_PAGINA + 2}")
