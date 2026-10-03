from pathlib import Path

import altair as alt
import pandas as pd
import streamlit as st

import seq_tools

DEMO_FILE = Path(__file__).parent / "sample_data" / "insulin_NM_000207.fasta"

st.set_page_config(page_title="BioSeq Analyzer", page_icon="🧬", layout="wide")


def load_example():
    st.session_state.pasted_text = DEMO_FILE.read_text()


# ---------- Sidebar ----------
with st.sidebar:
    st.header("About")
    st.write(
        "BioSeq Analyzer computes core properties of a DNA sequence. "
        "The analysis functions are written from scratch in Python and validated "
        "against Biopython with automated tests."
    )
    st.markdown("Built by **Kareem Damilare Oreoluwa**")
    st.markdown("[Source code on GitHub](https://github.com/damilare-kareem/bioseq-analyzer)")

# ---------- Input ----------
st.title("🧬 BioSeq Analyzer")
st.caption("Paste a DNA sequence or upload a FASTA file to analyse it.")

input_method = st.radio("Input method", ["Paste sequence", "Upload FASTA file"], horizontal=True)

raw_text = ""
if input_method == "Paste sequence":
    st.text_area(
        "DNA sequence (plain or FASTA format)",
        key="pasted_text",
        height=150,
        placeholder=">my_gene\nATGGCCATTGTAATGGGC...",
    )
    if DEMO_FILE.exists():
        st.button("Load example: human insulin mRNA (NM_000207.3)", on_click=load_example)
    raw_text = st.session_state.get("pasted_text", "")
else:
    uploaded = st.file_uploader("Upload a FASTA file", type=["fasta", "fa", "fna", "txt"])
    if uploaded is not None:
        raw_text = uploaded.getvalue().decode("utf-8-sig", errors="replace")

records = seq_tools.parse_fasta(raw_text)
if not records:
    st.info("Waiting for a sequence...")
    st.stop()

# A FASTA file can hold several sequences: let the user pick one
if len(records) > 1:
    choice = st.selectbox(
        f"{len(records)} sequences found. Choose one to analyse:",
        range(len(records)),
        format_func=lambda i: records[i][0],
    )
    header, sequence = records[choice]
else:
    header, sequence = records[0]

if not sequence:
    st.warning("This record has a header but no sequence.")
    st.stop()

invalid = seq_tools.find_invalid_characters(sequence)
if invalid:
    st.error(
        f"Found characters that aren't A, C, G or T: {', '.join(invalid)}. "
        "Please remove them and try again."
    )
    st.stop()

st.subheader(header)

# ---------- Summary metrics ----------
col1, col2, col3, col4 = st.columns(4)
col1.metric("Length", f"{len(sequence):,} bp")
col2.metric("GC content", f"{seq_tools.gc_content(sequence):.2f}%")
col3.metric("AT content", f"{seq_tools.at_content(sequence):.2f}%")
col4.metric("Complete codons", f"{len(sequence) // 3:,}")

composition_tab, sequences_tab, orf_tab, restriction_tab = st.tabs(
    ["📊 Composition", "🔁 Sequences", "🧬 Open reading frames", "✂️ Restriction sites"]
)

# ---------- Nucleotide composition ----------
with composition_tab:
    counts = seq_tools.nucleotide_counts(sequence)
    composition = pd.DataFrame({"Base": list(counts.keys()), "Count": list(counts.values())})
    composition["Percent"] = (composition["Count"] / len(sequence) * 100).round(2)

    chart_col, table_col = st.columns([2, 1])
    with chart_col:
        composition_chart = (
            alt.Chart(composition)
            .mark_bar(cornerRadiusTopLeft=4, cornerRadiusTopRight=4)
            .encode(
                x=alt.X("Base:N", sort=list("ACGT"), axis=alt.Axis(labelAngle=0)),
                y=alt.Y("Count:Q"),
                tooltip=["Base", "Count", alt.Tooltip("Percent:Q", format=".2f")],
            )
        )
        st.altair_chart(composition_chart)
    with table_col:
        st.dataframe(composition, hide_index=True)

# ---------- Sequence transformations ----------
with sequences_tab:
    st.markdown("**Reverse complement (5'→3')**")
    st.code(seq_tools.reverse_complement(sequence), language=None, wrap_lines=True)

    st.markdown("**mRNA (coding strand with T → U)**")
    st.code(seq_tools.transcribe(sequence), language=None, wrap_lines=True)

    st.markdown("**Straight translation (frame +1)**")
    st.code(seq_tools.translate(sequence), language=None, wrap_lines=True)
    st.caption(
        "* marks a stop codon. This reads straight through from the first base; "
        "see the Open reading frames tab for true start-to-stop coding regions."
    )

# ---------- Open reading frames ----------
with orf_tab:
    min_aa = st.slider("Minimum protein length (amino acids)", 5, 300, 30, step=5)
    orfs = seq_tools.find_orfs(sequence, min_protein_length=min_aa)

    if not orfs:
        st.info(f"No ORFs of at least {min_aa} amino acids found. Try lowering the minimum.")
    else:
        longest = orfs[0]
        st.success(
            f"Found {len(orfs)} ORF(s). Longest: frame {longest['Frame']}, "
            f"position {longest['Start']}–{longest['End']}, {longest['Length (aa)']} amino acids."
        )
        st.markdown("**Protein encoded by the longest ORF**")
        st.code(longest["Protein"], language=None, wrap_lines=True)

        st.markdown("**All ORFs**")
        st.dataframe(pd.DataFrame(orfs), hide_index=True)
    st.caption(
        "An ORF runs from an ATG start codon to the first in-frame stop codon (TAA, TAG or TGA). "
        "All six frames are scanned: +1 to +3 on the given strand, −1 to −3 on the reverse complement. "
        "For minus-strand ORFs, Start is larger than End."
    )

# ---------- Restriction sites ----------
with restriction_tab:
    sites = seq_tools.find_restriction_sites(sequence)
    summary = pd.DataFrame(
        {
            "Enzyme": list(sites.keys()),
            "Recognition site": [seq_tools.RESTRICTION_ENZYMES[e] for e in sites],
            "Number of sites": [len(p) for p in sites.values()],
            "Positions": [", ".join(map(str, p)) if p else "—" for p in sites.values()],
        }
    )

    cutters = summary[summary["Number of sites"] > 0]
    non_cutters = summary[summary["Number of sites"] == 0]["Enzyme"].tolist()

    if cutters.empty:
        st.info("None of the 12 enzymes cut this sequence.")
    else:
        site_rows = pd.DataFrame(
            [{"Enzyme": e, "Position": p} for e, positions in sites.items() for p in positions]
        )
        site_map = (
            alt.Chart(site_rows)
            .mark_tick(thickness=4, size=24)
            .encode(
                x=alt.X(
                    "Position:Q",
                    scale=alt.Scale(domain=[1, len(sequence)]),
                    title="Position in sequence (bp)",
                ),
                y=alt.Y("Enzyme:N", title=None),
                tooltip=["Enzyme", "Position"],
            )
        )
        st.markdown("**Restriction map**")
        st.altair_chart(site_map)
        st.dataframe(cutters, hide_index=True)

    if non_cutters:
        st.markdown(f"**Non-cutters** (no site in this sequence): {', '.join(non_cutters)}")
        st.caption(
            "Non-cutters are useful for cloning: they can open a vector "
            "without cutting the insert."
        )