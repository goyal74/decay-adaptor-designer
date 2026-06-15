# Decay Adaptor Designer — in-silico programmable RNA "death tag" design

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20707857.svg)](https://doi.org/10.5281/zenodo.20707857)

> Run it in your browser (no coding): **https://goyal74-decay-adaptor-designer-app-ysuuno.streamlit.app**
> Developed by **Ravi Goyal, MD, PhD, Epigenuity LLC**. MIT License.

Designs a modular **decay-pathway recruitment adaptor** — a programmable RNA "death tag" — that
redirects the cell's own RNA-surveillance machinery to destroy a disease-driving transcript. Name the
**condition** and the **transcript to destroy**; the tool finds an accessible target site, designs the
antisense **RNA-binding element** (15–25 nt, 2′-MOE/LNA chemistry), and joins it to a constant
**Decay-Triggering Module (DTM)** matched to the target's compartment:

| Compartment | Route (DTM strategy) | Endogenous effector |
|-------------|----------------------|---------------------|
| Nuclear (lncRNA / repeat) — **lead** | D/E — nuclear-exosome routing (NEXT/PAXT-mimic or structure disruption) | nuclear RNA exosome |
| Cytoplasmic mRNA | A/C — NMD recruitment (UPF1-recruiting / EJC-mimic) | UPF1 / NMD |
| 5′UTR-proximal | B — uORF cassette (NMD-eligible context) | uORF-triggered NMD |
| Splice-switch | F — pseudo-exon → frameshift PTC | endogenous NMD |

It reports, using ViennaRNA:
- **Target-site accessibility** (mean unpaired probability of the chosen window, local folding),
- **Antisense:target duplex free energy** (ΔG, kcal/mol),
- **Nearest-neighbor Tm** of the binding element (unmodified estimate) and **GC content**,
- the assembled **adaptor sequence** and the **compartment → decay-route** decision.

## Scope / honesty
This is **in-silico target-site + binding design only** (the computationally tractable part). The actual
**decay recruitment** — UPF1/SMD in the cytoplasm, the NEXT/PAXT-fed **nuclear exosome** in the nucleus —
is a **Phase II wet-lab readout**, precedented by artificial tethering of UPF1 (sufficient to trigger
decay) and by the well-characterized nuclear-exosome adaptor complexes. The DTM motifs here are **constant
scaffold stand-ins**, not optimized effector sequences. Transcriptome-wide BLAST uniqueness is run
separately. **Not a medical device.**

## Install & run
```bash
pip install ViennaRNA numpy matplotlib pandas streamlit requests
streamlit run app.py            # web app
# or command line:
python3 decay_adaptor_designer.py --target MALAT1_segment.fa --compartment nuclear --out adaptor.json
python3 make_adaptor_figure.py adaptor.json adaptor.png
```

## Files
| File | Purpose |
|------|---------|
| `decay_adaptor_designer.py` | core: accessible-window selection, antisense binding-element design, duplex ΔG / Tm / GC, compartment→DTM matching, adaptor assembly |
| `make_adaptor_figure.py` | target-site duplex + modular-adaptor schematic + metrics/decision figure |
| `app.py` | Streamlit web app (condition + target gene; infers compartment, suggests DTM) |

Part of the NIH TRDNT Challenge submission (Programmable RNA "Death Tags" — Decay-Pathway Recruitment Adaptors, Epigenuity LLC).
