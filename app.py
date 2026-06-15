#!/usr/bin/env python3
"""
Decay Adaptor Designer - web app (Epigenuity LLC)
=================================================
Name a disease condition and the disease-driving transcript you want to destroy.
The tool finds an accessible target site, designs the antisense RNA-binding element
(2'-MOE/LNA chemistry), and attaches the Decay-Triggering Module (DTM) matched to the
transcript's compartment - routing nuclear RNAs to the nuclear exosome and cytoplasmic
mRNAs to NMD - then reports predicted accessibility, duplex free energy, Tm, and GC.

Developed by Ravi Goyal, MD, PhD, Epigenuity LLC. In-silico target-site + binding
design only; the decay-recruitment step is a Phase II wet-lab readout. Not a medical device.
"""
import streamlit as st
import matplotlib
matplotlib.use("Agg")
import decay_adaptor_designer as da
import requests

st.set_page_config(page_title="Decay Adaptor Designer · Epigenuity LLC", layout="wide")

st.title("Decay Adaptor Designer — programmable RNA “death tags”")
st.markdown("**Developed by Ravi Goyal, MD, PhD · Epigenuity LLC** &nbsp;·&nbsp; "
            "[github.com/goyal74/decay-adaptor-designer](https://github.com/goyal74/decay-adaptor-designer)")
st.caption("Name the **condition** and the **transcript to destroy**. The tool designs a modular adaptor — an "
           "antisense RNA-binding element joined to a Decay-Triggering Module that redirects the cell's own "
           "RNA-surveillance machinery (nuclear exosome for nuclear RNAs, NMD for cytoplasmic mRNAs) onto that "
           "transcript. In-silico design only; confirm experimentally.")


def fetch_by_name(symbol, species="homo_sapiens", timeout=25):
    symbol = (symbol or "").strip()
    if not symbol:
        return None
    try:
        b = "https://rest.ensembl.org"
        r = requests.get(f"{b}/lookup/symbol/{species}/{symbol}", params={"expand": "1"},
                         headers={"Content-Type": "application/json"}, timeout=timeout)
        if r.ok:
            g = r.json(); tid = (g.get("canonical_transcript") or "").split(".")[0]
            if not tid and g.get("Transcript"):
                tid = g["Transcript"][0]["id"]
            if tid:
                s = requests.get(f"{b}/sequence/id/{tid}", params={"type": "cdna"},
                                 headers={"Content-Type": "text/x-fasta"}, timeout=timeout)
                if s.ok and s.text.lstrip().startswith(">"):
                    return "".join(l for l in s.text.splitlines() if not l.startswith(">"))
    except Exception:
        pass
    return None


# Condition -> (target gene, compartment, one-line rationale). The compartment selects
# the matched decay route; "nuclear" -> exosome (lead), "cytoplasmic" -> NMD.
CONDITIONS = {
    "— choose a condition —": ("", "nuclear", ""),
    "Metastatic cancer — MALAT1 (lead)": ("MALAT1", "nuclear",
        "Nuclear, triple-helix-stabilized lncRNA → nuclear-exosome routing (Strategy D/E)."),
    "Myotonic dystrophy (DM1) — DMPK": ("DMPK", "nuclear",
        "Largely nuclear CUG-repeat transcript → nuclear-exosome routing; frees sequestered MBNL."),
    "Huntington's disease — mutant HTT": ("HTT", "cytoplasmic",
        "Allele-selective NMD via a linked SNP spares wild-type HTT (Strategy A/C)."),
    "KRAS-driven cancer — mutant KRAS": ("KRAS", "cytoplasmic",
        "Codon-selective degradation of an undruggable oncogene mRNA via NMD (Strategy A/C)."),
    "Neuroblastoma — MYCN": ("MYCN", "cytoplasmic",
        "Over-expressed cytoplasmic oncogene mRNA → NMD recruitment."),
    "Lung adenocarcinoma — NEAT1 (nuclear lncRNA)": ("NEAT1", "nuclear",
        "Nuclear paraspeckle lncRNA → nuclear-exosome routing (Strategy D/E)."),
}


def _apply_condition():
    k = st.session_state.get("cond_pick", "")
    gene, comp, _ = CONDITIONS.get(k, ("", "nuclear", ""))
    if gene:
        st.session_state.tgt_name = gene
        st.session_state.tgt_seq = ""
        st.session_state.comp_pick = COMP_LABEL[comp]


COMP_LABEL = {"nuclear": "Nuclear (lncRNA / repeat) → exosome",
              "cytoplasmic": "Cytoplasmic mRNA → NMD",
              "5utr": "5'UTR-proximal → uORF cassette",
              "splice": "Splice-switch death tag (pseudo-exon → PTC)"}
LABEL_COMP = {v: k for k, v in COMP_LABEL.items()}


def _load_example():
    st.session_state.tgt_name = "MALAT1"
    st.session_state.tgt_seq = ""
    st.session_state.comp_pick = COMP_LABEL["nuclear"]


with st.sidebar:
    st.markdown("**Destroy one disease transcript, spare the rest.** Name the *condition* and the "
                "*transcript to destroy*; the device recruits the cell's own machinery to degrade it.")
    st.button("⚡ Load an example (MALAT1)", on_click=_load_example, use_container_width=True)
    st.divider()

    st.subheader("① Condition")
    st.selectbox("Start from a condition (fills the target + route below)", list(CONDITIONS.keys()),
                 key="cond_pick", on_change=_apply_condition)
    k = st.session_state.get("cond_pick", "")
    if CONDITIONS.get(k, ("", "", ""))[2]:
        st.caption(CONDITIONS[k][2])
    st.divider()

    st.subheader("② Transcript to destroy")
    st.caption("The disease-driving RNA to degrade — a nuclear lncRNA, a repeat-expansion transcript, "
               "or a mutant/over-expressed mRNA. Type a gene symbol or paste a sequence under Advanced.")
    tgt_name = st.text_input("Target gene", placeholder="e.g. MALAT1, HTT, KRAS, DMPK", key="tgt_name")
    st.divider()

    st.subheader("③ Decay route")
    st.selectbox("Compartment → matched decay machinery", list(LABEL_COMP.keys()), key="comp_pick")
    st.caption("Nuclear targets route to the nuclear exosome (the near-term lead); cytoplasmic mRNAs "
               "recruit NMD. The route is matched to where the target lives.")

    with st.expander("Advanced — paste sequence / settings"):
        tgt_seq = st.text_area("Target RNA (instead of a gene name)", height=80, key="tgt_seq")
        be_len = st.slider("Binding-element length (nt)", 15, 25, 20)
    go = st.button("Design adaptor", type="primary", use_container_width=True)
    st.caption("© 2026 Epigenuity LLC · MIT License")


def resolve(text, name):
    if name.strip():
        seq = fetch_by_name(name)
        if not seq:
            st.error(f"Could not find '{name}' as a gene. Use a gene symbol (e.g. MALAT1, KRAS), "
                     "pick a condition from the dropdown, or paste a sequence under Advanced."); st.stop()
        st.success(f"Fetched {name} ({len(seq)} nt)")
        return seq
    return text


if go:
    comp = LABEL_COMP.get(st.session_state.get("comp_pick", ""), "nuclear")
    with st.spinner("Resolving target sequence…"):
        target = resolve(st.session_state.get("tgt_seq", ""), st.session_state.get("tgt_name", ""))
    if not target:
        st.error("Name a target gene or paste a sequence under Advanced."); st.stop()
    with st.spinner("Finding accessible site and designing the adaptor…"):
        r = da.design(target, compartment=comp, be_len=be_len)
        import tempfile, json
        jp = tempfile.NamedTemporaryFile("w", suffix=".json", delete=False)
        json.dump(r, jp); jp.close()
        import make_adaptor_figure as mf
        fig_path = jp.name + ".png"
        mf.draw(jp.name, fig_path, "Decay Adaptor Designer — in-silico death-tag design")

    c1, c2, c3 = st.columns(3)
    c1.metric("Target-site accessibility", f"{r['target_site_accessibility']:.2f}")
    c2.metric("Antisense:target duplex ΔG", f"{r['duplex_dG_kcal']:.0f} kcal/mol")
    c3.metric("Binding-element Tm (unmod.)", f"{r['binding_element_tm_C']:.0f} °C")

    st.image(fig_path, use_container_width=True)

    st.subheader(f"Decay route: {r['dtm_strategy']}")
    st.caption(f"Recruited effector: {r['dtm_effector']}")
    st.markdown("**Designed adaptor (5′→3′)**")
    st.code(r["adaptor_5to3"], language="text")
    st.markdown(
        f"- RNA-binding element (antisense, 2′-MOE/LNA): `{r['binding_element_5to3']}` "
        f"(GC {r['binding_element_gc_pct']}%)\n"
        f"- Decay-Triggering Module motif: `{r['dtm_motif_5to3'] or '(splice-switch ASO — no appended motif)'}`\n"
        f"- Target window (accessible site): `{r['target_window_5to3']}` (start {r['target_window_start']})")
    st.download_button("Download design (JSON)", open(jp.name).read(), file_name="decay_adaptor.json")
    st.info("In-silico target-site + binding design only. Decay recruitment (UPF1/SMD; nuclear exosome via "
            "NEXT/PAXT) is a Phase II wet-lab readout, precedented by artificial UPF1 tethering and the "
            "characterized NEXT/PAXT-fed exosome. The DTM motif is a constant scaffold stand-in. Confirm experimentally.")
else:
    st.write("**New here?** Click **⚡ Load an example (MALAT1)** in the sidebar, then **Design adaptor** — "
             "or name your own target:")
    st.markdown(
        "- **① Condition** — the disease whose transcript you want to destroy.\n"
        "- **② Transcript to destroy** — the disease-driving RNA (a nuclear lncRNA, a repeat-expansion "
        "transcript, or a mutant/over-expressed mRNA).\n"
        "- **③ Decay route** — matched to where the target lives: nuclear → exosome (lead), cytoplasmic → NMD.\n\n"
        "The tool returns the designed antisense binding element, the matched Decay-Triggering Module, the "
        "assembled adaptor, and predicted accessibility, duplex free energy, Tm, and GC.")

st.divider()
st.caption("Decay Adaptor Designer · Developed by **Ravi Goyal, MD, PhD, Epigenuity LLC** "
           "(Tucson, AZ) · © 2026 Epigenuity LLC · MIT License")
