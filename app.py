import pandas as pd
import streamlit as st

import seq_tools

EXAMPLE = ">example_gene Demo sequence\nATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"

st.set_page_config(page_title="BioSeq Analyzer", page_icon="🧬", layout="wide")


def load_example():
    st.session_state.pasted_text = EXAMPLE


# ---------- Sidebar ----------
with st.sidebar:
    st.header("About")
    st.write(
        "BioSeq Analyzer computes core properties of a DNA sequence. "
        "The analysis functions are written from scratch and validated against Biopython."
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
    st.button("Load example sequence", on_click=load_example)
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

# ---------- Nucleotide composition ----------
st.subheader("Nucleotide composition")
counts = seq_tools.nucleotide_counts(sequence)
composition = pd.DataFrame({"Base": list(counts.keys()), "Count": list(counts.values())})
composition["Percent"] = (composition["Count"] / len(sequence) * 100).round(2)

chart_col, table_col = st.columns([2, 1])
with chart_col:
    st.bar_chart(composition, x="Base", y="Count")
with table_col:
    st.dataframe(composition, hide_index=True)

# ---------- Sequence transformations ----------
st.subheader("Reverse complement (5'→3')")
st.code(seq_tools.reverse_complement(sequence), language=None, wrap_lines=True)

st.subheader("mRNA (coding strand with T → U)")
st.code(seq_tools.transcribe(sequence), language=None, wrap_lines=True)

st.subheader("Protein (reading frame 1)")
st.code(seq_tools.translate(sequence), language=None, wrap_lines=True)
st.caption(
    "* marks a stop codon. This reads straight through from the first base; "
    "the ORF finder (coming next) will locate true start-to-stop reading frames."
)