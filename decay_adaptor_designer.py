#!/usr/bin/env python3
"""
Decay Adaptor Designer - programmable RNA "death tag" designer (Stage-1, in silico)
===================================================================================
Developed by Ravi Goyal, MD, PhD, Epigenuity LLC.

Input  : a target transcript (gene name or RNA sequence) and its subcellular
         compartment (cytoplasmic / nuclear / 5'UTR-proximal / splice-switch).
Output : a modular decay adaptor = a target-specific RNA-binding element
         (15-25 nt antisense, 2'-MOE/LNA chemistry) joined through a linker to a
         constant Decay-Triggering Module (DTM) matched to the compartment:
           - nuclear      -> nuclear-exosome routing (Strategy D/E; NEXT/PAXT-mimic
                             or structure-disruption of a protective fold)
           - cytoplasmic  -> NMD recruitment (Strategy A/C; UPF1-recruiting / EJC-mimic)
           - 5'UTR        -> uORF cassette creating an NMD-eligible context (Strategy B)
           - splice       -> splice-switch "death tag": pseudo-exon -> frameshift PTC
                             converts the transcript into an endogenous NMD substrate (Strategy F)
         with a ViennaRNA-predicted target-site accessibility, the antisense:target
         duplex free energy, a nearest-neighbor melting temperature, GC content, and
         a compartment->DTM decision.

In-silico DESIGN + TARGET-SITE / BINDING PREDICTION only. The actual decay
recruitment (UPF1/SMD, nuclear exosome via NEXT/PAXT) is a Phase II wet-lab readout,
precedented by artificial tethering of UPF1 (sufficient to trigger decay) and by the
well-characterized NEXT/PAXT-fed nuclear exosome. The DTM motifs here are constant
scaffold stand-ins, not optimized effector sequences. Not a medical device.
"""
import json, argparse, math
try:
    import RNA
except Exception:
    RNA = None

COMP = {"A": "U", "U": "A", "G": "C", "C": "G", "T": "A", "N": "N"}

# Constant Decay-Triggering Module (DTM) stand-ins, one per compartment strategy.
# These are placeholders for the effector motif; the functional sequence is a
# Phase II wet-lab deliverable. The designer's job is target-site + binding design.
DTM = {
    "nuclear":     {"strategy": "D/E - nuclear-exosome routing (NEXT/PAXT-mimic or structure disruption)",
                    "motif": "GAUCGGAAGAGC", "effector": "nuclear RNA exosome (NEXT/PAXT)"},
    "cytoplasmic": {"strategy": "A/C - NMD recruitment (UPF1-recruiting / EJC-mimic)",
                    "motif": "GAAACAGGUGA", "effector": "UPF1 / NMD"},
    "5utr":        {"strategy": "B - uORF cassette (NMD-eligible context)",
                    "motif": "AUGUAAUGAUAA", "effector": "uORF-triggered NMD"},
    "splice":      {"strategy": "F - splice-switch death tag (pseudo-exon -> frameshift PTC)",
                    "motif": "", "effector": "endogenous NMD (splice-switch ASO; no appended motif)"},
}
LINKER = "UUUU"

# SantaLucia (1998) unified nearest-neighbor DNA parameters (dH kcal/mol, dS cal/mol-K).
_NN_H = {"AA": -7.9, "AT": -7.2, "TA": -7.2, "CA": -8.5, "GT": -8.4, "CT": -7.8,
         "GA": -8.2, "CG": -10.6, "GC": -9.8, "GG": -8.0,
         "AC": -8.4, "TG": -8.5, "AG": -7.8, "TC": -8.2, "TT": -7.9, "CC": -8.0}
_NN_S = {"AA": -22.2, "AT": -20.4, "TA": -21.3, "CA": -22.7, "GT": -22.4, "CT": -21.0,
         "GA": -22.2, "CG": -27.2, "GC": -24.4, "GG": -19.9,
         "AC": -22.4, "TG": -22.7, "AG": -21.0, "TC": -22.2, "TT": -22.2, "CC": -19.9}

def revcomp(s):
    s = s.upper().replace("T", "U")
    return "".join(COMP.get(b, "N") for b in reversed(s))

def clean(seq):
    return "".join(c for c in seq.upper().replace("T", "U") if c in "ACGU")

def gc_percent(s):
    s = s.upper().replace("U", "T")
    return 100.0 * sum(c in "GC" for c in s) / max(1, len(s))

def nn_tm(seq, conc_nM=200.0, na_M=1.0):
    """Nearest-neighbor Tm (deg C) for the antisense oligo as a DNA analog (SantaLucia 1998).
    Unmodified estimate; 2'-MOE/LNA chemistry raises the effective Tm well above this."""
    s = seq.upper().replace("U", "T")
    if len(s) < 2:
        return float("nan")
    dH = 0.2 + (2.2 if s[0] in "AT" else 0.0) + (2.2 if s[-1] in "AT" else 0.0)  # init + terminal-AT
    dS = -5.7 + (6.9 if s[0] in "AT" else 0.0) + (6.9 if s[-1] in "AT" else 0.0)
    for i in range(len(s) - 1):
        nn = s[i:i+2]
        dH += _NN_H.get(nn, -8.0); dS += _NN_S.get(nn, -22.0)
    R = 1.987
    ct = conc_nM * 1e-9 / 4.0  # for self-non-complementary duplex
    tm = (dH * 1000.0) / (dS + R * math.log(ct)) - 273.15
    tm += 16.6 * math.log10(max(1e-3, na_M))  # salt correction
    return round(tm, 1)

def accessible_window(seq, L=20, W=120):
    """Pick the most single-stranded (accessible) window of length L by local folding."""
    seq = clean(seq); n = len(seq)
    if n <= L:
        return 0, seq, 1.0
    if RNA is None:
        return 0, seq[:L], float("nan")
    Wl = min(W, n)
    pl = RNA.pfl_fold(seq, Wl, min(Wl, 100), 1e-4)
    paired = [0.0] * n
    for e in pl:
        try: i, j, p = e.i, e.j, e.p
        except AttributeError: i, j, p = e[0], e[1], e[2]
        paired[i-1] += p; paired[j-1] += p
    unp = [max(0.0, 1 - x) for x in paired]
    best_i, best = 0, -1
    for i in range(0, n - L + 1):
        s = sum(unp[i:i+L]) / L
        if s > best:
            best, best_i = s, i
    return best_i, seq[best_i:best_i+L], round(best, 3)

def duplex_dg(binding_element, target_window):
    """ViennaRNA antisense:target hybrid free energy (kcal/mol; more negative = tighter)."""
    if RNA is None:
        return float("nan")
    try:
        d = RNA.duplexfold(binding_element, target_window)
        return round(d.energy, 1)
    except Exception:
        return float("nan")

def design(target, compartment="nuclear", be_len=20):
    target = clean(target)
    compartment = compartment.lower()
    if compartment not in DTM:
        compartment = "nuclear"
    start, window, acc = accessible_window(target, L=be_len)
    binding_element = revcomp(window)               # antisense RNA-binding element (2'-MOE/LNA chimera)
    dtm = DTM[compartment]
    if compartment == "splice":
        adaptor = binding_element                   # splice-switch ASO; no appended effector motif
    else:
        adaptor = binding_element + LINKER + dtm["motif"]
    res = {
        "compartment": compartment,
        "dtm_strategy": dtm["strategy"],
        "dtm_effector": dtm["effector"],
        "target_window_5to3": window,
        "target_window_start": start,
        "binding_element_5to3": binding_element,
        "binding_element_len": len(binding_element),
        "dtm_motif_5to3": dtm["motif"],
        "linker": "" if compartment == "splice" else LINKER,
        "adaptor_5to3": adaptor,
        "adaptor_len": len(adaptor),
        "target_site_accessibility": acc,
        "duplex_dG_kcal": duplex_dg(binding_element, window),
        "binding_element_gc_pct": round(gc_percent(binding_element), 1),
        "binding_element_tm_C": nn_tm(binding_element),
        "design_rules": {
            "accessible_window_min_nt": 15,
            "binding_element_nt": "15-25",
            "chemistry": "2'-MOE/LNA, phosphorothioate",
            "target_affinity": "Kd < 50 nM (predicted)",
            "tm_target_C": ">= 60 (with 2'-MOE/LNA; NN estimate below is for the unmodified oligo)",
            "uniqueness": "BLAST transcriptome-wide: <= 70% off-target identity (run separately)",
        },
        "qc_note": ("In-silico target-site + binding design only. Decay recruitment "
                    "(UPF1/SMD, nuclear exosome via NEXT/PAXT) is a Phase II wet-lab "
                    "readout; the DTM motif is a constant scaffold stand-in."),
    }
    return res

def main():
    ap = argparse.ArgumentParser(description="Decay Adaptor Designer (programmable RNA death tag, in silico)")
    ap.add_argument("--target", required=True, help="target RNA: sequence or .fa/.txt path")
    ap.add_argument("--compartment", default="nuclear",
                    choices=["nuclear", "cytoplasmic", "5utr", "splice"],
                    help="target compartment -> selects the matched DTM strategy")
    ap.add_argument("--be-len", type=int, default=20, help="binding-element length (15-25 nt)")
    ap.add_argument("--out", default="adaptor.json")
    a = ap.parse_args()
    def load(x):
        import os
        if os.path.exists(x):
            t = open(x).read()
            return "".join(l for l in t.splitlines() if not l.startswith(">"))
        return x
    r = design(load(a.target), compartment=a.compartment, be_len=a.be_len)
    json.dump(r, open(a.out, "w"), indent=2)
    print("compartment:", r["compartment"], "->", r["dtm_strategy"])
    print("binding element (5'->3'):", r["binding_element_5to3"])
    print("full adaptor   (5'->3'):", r["adaptor_5to3"])
    print(f"accessibility={r['target_site_accessibility']}  duplex_dG={r['duplex_dG_kcal']} kcal/mol  "
          f"Tm(unmod)={r['binding_element_tm_C']} C  GC={r['binding_element_gc_pct']}%")
    print("wrote", a.out)

if __name__ == "__main__":
    main()
