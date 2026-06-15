#!/usr/bin/env python3
"""Render a decay-adaptor design: target-site duplex + modular adaptor + metrics/decision."""
import json, sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
NAVY, BLUE, RED, GREY, AMBER, TEAL, GREEN = "#1F3864", "#2E6DB4", "#C0504D", "#BBBBBB", "#E0922E", "#2E9B8F", "#2E8B57"

PAIR = {"A": "U", "U": "A", "G": "C", "C": "G"}


def draw(jpath, out, title="Decay Adaptor Designer — in-silico death-tag design"):
    d = json.load(open(jpath))
    win = d["target_window_5to3"]
    be = d["binding_element_5to3"]            # antisense, 5'->3'
    be_anti = be[::-1]                         # 3'->5' so it aligns under the 5'->3' target
    fig = plt.figure(figsize=(11, 5.4))

    # ---- Top-left: target site / antisense duplex ----
    ax = fig.add_axes([0.04, 0.55, 0.52, 0.34]); ax.axis("off")
    ax.set_xlim(0, max(len(win), 1)); ax.set_ylim(0, 3)
    ax.text(0, 2.78, "Target site (accessible window, 5′→3′)", fontsize=9, color=NAVY)
    for i, b in enumerate(win):
        ax.text(i + 0.5, 2.25, b, ha="center", va="center", fontsize=9.5, family="monospace", color="#1A1A1A")
    for i in range(len(win)):
        comp_ok = i < len(be_anti) and PAIR.get(win[i]) == be_anti[i]
        ax.plot([i + 0.5, i + 0.5], [2.05, 1.55], "-", color=(GREEN if comp_ok else GREY), lw=1.0)
    for i, b in enumerate(be_anti):
        ax.text(i + 0.5, 1.3, b, ha="center", va="center", fontsize=9.5, family="monospace", color=BLUE)
    ax.text(0, 0.85, "Antisense RNA-binding element (3′→5′; 2′-MOE/LNA)", fontsize=9, color=BLUE)

    # ---- Bottom-left: modular adaptor schematic ----
    axm = fig.add_axes([0.04, 0.06, 0.52, 0.38]); axm.axis("off")
    axm.set_xlim(0, 10); axm.set_ylim(0, 3)
    axm.add_patch(FancyBboxPatch((0.2, 1.1), 4.2, 1.0, boxstyle="round,pad=0.05,rounding_size=0.12",
                                 fc="#DCE6F4", ec=BLUE, lw=2))
    axm.text(2.3, 1.6, "RNA-binding element\n(antisense, target-specific)", ha="center", va="center",
             fontsize=8.5, color=NAVY)
    splice = (d["compartment"] == "splice")
    if not splice:
        axm.annotate("", xy=(5.5, 1.6), xytext=(4.5, 1.6), arrowprops=dict(arrowstyle="-", color=GREY, lw=2))
        axm.text(5.0, 1.95, "linker", ha="center", fontsize=7.5, color=GREY, style="italic")
        axm.add_patch(FancyBboxPatch((5.6, 1.1), 4.0, 1.0, boxstyle="round,pad=0.05,rounding_size=0.12",
                                     fc="#FAEBD2", ec=AMBER, lw=2))
        axm.text(7.6, 1.6, "Decay-Triggering Module\n(constant, compartment-matched)", ha="center",
                 va="center", fontsize=8.5, color="#8A5A12")
    else:
        axm.text(7.4, 1.6, "splice-switch ASO\n(pseudo-exon → frameshift PTC → NMD)", ha="center",
                 va="center", fontsize=8.5, color="#8A5A12")
    axm.text(0.2, 0.4, f"Assembled adaptor: {d['adaptor_len']} nt  ·  5′-{d['adaptor_5to3']}-3′",
             fontsize=7.6, family="monospace", color="#444")
    axm.set_title("Modular adaptor", fontsize=10, color=NAVY, loc="left")

    # ---- Right-top: metrics ----
    ax2 = fig.add_axes([0.63, 0.52, 0.33, 0.37])
    labels = ["accessibility\n(0–1)", "|duplex ΔG|\n(kcal/mol)", "Tm unmod.\n(°C)", "GC\n(%)"]
    acc = d["target_site_accessibility"]; dg = abs(d["duplex_dG_kcal"])
    tm = d["binding_element_tm_C"]; gc = d["binding_element_gc_pct"]
    vals = [acc, dg, tm, gc]
    cols = [TEAL, BLUE, AMBER, GREY]
    bars = ax2.bar(labels, vals, color=cols, width=0.62)
    for b, v in zip(bars, vals):
        ax2.text(b.get_x() + b.get_width() / 2, v, f"{v:.1f}" if v < 5 else f"{v:.0f}",
                 ha="center", va="bottom", fontsize=8.5)
    ax2.set_ylim(0, max(vals) * 1.18)
    ax2.set_title("Predicted binding metrics", fontsize=9.5, color=NAVY)
    ax2.tick_params(labelsize=7.5); ax2.set_yticks([])
    for sp in ["top", "right", "left"]:
        ax2.spines[sp].set_visible(False)

    # ---- Right-bottom: decision + QC ----
    ax3 = fig.add_axes([0.63, 0.05, 0.33, 0.36]); ax3.axis("off")
    comp_name = {"nuclear": "Nuclear", "cytoplasmic": "Cytoplasmic",
                 "5utr": "5′UTR-proximal", "splice": "Splice-switch"}[d["compartment"]]
    txt = (f"Compartment → decay route\n"
           f"  {comp_name}\n"
           f"  → {d['dtm_strategy']}\n"
           f"  effector: {d['dtm_effector']}")
    ax3.text(0, 0.98, txt, fontsize=8.3, va="top", color="#1A1A1A")
    ax3.text(0, 0.30, "In-silico target-site + binding design only;\ndecay recruitment is a Phase II wet-lab readout.",
             fontsize=7.5, va="top", color=GREY, style="italic")

    fig.suptitle(title, fontsize=12, color=NAVY, fontweight="bold", y=1.02)
    fig.savefig(out, dpi=150, bbox_inches="tight"); plt.close(fig); print("wrote", out)


if __name__ == "__main__":
    draw(sys.argv[1], sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else
         "Decay Adaptor Designer — in-silico death-tag design")
