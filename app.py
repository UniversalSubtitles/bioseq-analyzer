import streamlit as st

# Browser tab title and icon
st.set_page_config(page_title="BioSeq Analyzer", page_icon="🧬")

st.title("🧬 BioSeq Analyzer")
st.write("Paste a DNA sequence below to analyse it.")

# Text box where the user pastes a sequence
raw_sequence = st.text_area(
	"DNA sequence",
	height=150,
	placeholder="e.g. ATGGCCATTGTAATGGGCCGC",
)

if raw_sequence:
	# Clean the input: remove spaces/line breaks, convert to capitals
	sequence = "".join(raw_sequence.split()).upper()
	st.metric("Sequence length", f"{len(sequence)} bp")
