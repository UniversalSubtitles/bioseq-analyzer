"""Core DNA sequence analysis functions for BioSeq Analyzer."""

VALID_BASES = "ACGT"

# Watson-Crick base pairing
COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G"}

# The standard genetic code (NCBI table 1), laid out in T-C-A-G order,
# exactly like the codon table in a biochemistry textbook.
BASES = "TCAG"
CODONS = [a + b + c for a in BASES for b in BASES for c in BASES]
AMINO_ACIDS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE = dict(zip(CODONS, AMINO_ACIDS))


def clean_sequence(raw_text):
	"""Remove FASTA header lines, spaces and line breaks; convert to capitals."""
	lines = raw_text.splitlines()
	sequence_lines = [line for line in lines if not line.strip().startswith(">")]
	return "".join("".join(sequence_lines).split()).upper()


def find_invalid_characters(sequence):
	"""Return any characters that are not A, C, G or T (sorted)."""
	return sorted(set(sequence) - set(VALID_BASES))


def nucleotide_counts(sequence):
	"""Count each of the four bases."""
	return {base: sequence.count(base) for base in VALID_BASES}


def gc_content(sequence):
	"""Percentage of bases that are G or C."""
	if not sequence:
		return 0.0
	gc = sequence.count("G") + sequence.count("C")
	return gc / len(sequence) * 100


def at_content(sequence):
	"""Percentage of bases that are A or T."""
	if not sequence:
		return 0.0
	at = sequence.count("A") + sequence.count("T")
	return at / len(sequence) * 100


def reverse_complement(sequence):
	"""Complement every base, then reverse: the opposite strand read 5'->3'."""
	return "".join(COMPLEMENT[base] for base in reversed(sequence))


def transcribe(sequence):
	"""Coding-strand DNA -> mRNA: same sequence, with T replaced by U."""
	return sequence.replace("T", "U")


def translate(sequence):
	"""Translate codon by codon from the first base.
	Stop codons appear as '*'. Leftover bases that don't form a full codon are ignored."""
	protein = ""
	for i in range(0, len(sequence) - 2, 3):
		codon = sequence[i:i + 3]
		protein += CODON_TABLE[codon]
	return protein


def parse_fasta(raw_text):
    """Split FASTA text into a list of (header, sequence) records.
    Text with no '>' header line is treated as one unnamed sequence."""
    records = []
    header, chunks = None, []
    for line in raw_text.splitlines():
        line = line.strip()
        if line.startswith(">"):
            if header is not None or chunks:
                records.append((header or "Unnamed sequence", "".join(chunks).upper()))
            header, chunks = line[1:].strip() or "Unnamed sequence", []
        elif line:
            chunks.append("".join(line.split()))
    if header is not None or chunks:
        records.append((header or "Unnamed sequence", "".join(chunks).upper()))
    return records

STOP_CODONS = {"TAA", "TAG", "TGA"}


def find_orfs(sequence, min_protein_length=30):
    """Find open reading frames (ATG ... stop codon) in all six reading frames.

    Each frame is scanned from left to right: at the first ATG we read codons
    until a stop codon, record the ORF, then continue scanning after that stop.
    ORFs with no stop codon before the sequence ends are ignored.
    Positions are 1-based on the sequence as given. For minus-strand ORFs,
    Start is larger than End because that strand is read right to left."""
    seq_len = len(sequence)
    orfs = []
    for strand, strand_seq in (("+", sequence), ("-", reverse_complement(sequence))):
        for frame in range(3):
            i = frame
            while i + 3 <= seq_len:
                if strand_seq[i:i + 3] == "ATG":
                    stop_index = None
                    for j in range(i, seq_len - 2, 3):
                        if strand_seq[j:j + 3] in STOP_CODONS:
                            stop_index = j
                            break
                    if stop_index is None:
                        break  # no stop codon downstream in this frame
                    protein = translate(strand_seq[i:stop_index])
                    if len(protein) >= min_protein_length:
                        if strand == "+":
                            start, end = i + 1, stop_index + 3
                        else:
                            start, end = seq_len - i, seq_len - stop_index - 2
                        orfs.append({
                            "Frame": f"{strand}{frame + 1}",
                            "Start": start,
                            "End": end,
                            "Length (nt)": stop_index + 3 - i,
                            "Length (aa)": len(protein),
                            "Protein": protein,
                        })
                    i = stop_index  # resume scanning after this ORF's stop codon
                i += 3
    return sorted(orfs, key=lambda orf: orf["Length (aa)"], reverse=True)


# Common type II restriction enzymes and their recognition sites.
# All of these sites are palindromic: they read the same 5'->3' on both strands,
# so searching one strand finds every site.
RESTRICTION_ENZYMES = {
    "BamHI": "GGATCC",
    "EcoRI": "GAATTC",
    "HindIII": "AAGCTT",
    "KpnI": "GGTACC",
    "NcoI": "CCATGG",
    "NdeI": "CATATG",
    "NotI": "GCGGCCGC",
    "PstI": "CTGCAG",
    "SacI": "GAGCTC",
    "SmaI": "CCCGGG",
    "XbaI": "TCTAGA",
    "XhoI": "CTCGAG",
}


def find_restriction_sites(sequence, enzymes=RESTRICTION_ENZYMES):
    """Return {enzyme: [1-based start positions of its recognition site]}."""
    sites = {}
    for enzyme, site in enzymes.items():
        positions = []
        index = sequence.find(site)
        while index != -1:
            positions.append(index + 1)
            index = sequence.find(site, index + 1)
        sites[enzyme] = positions
    return sites

if __name__ == "__main__":
	example = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
	print(f"Length: {len(example)} bp")
	print(f"Counts: {nucleotide_counts(example)}")
	print(f"GC content: {gc_content(example):.2f}%")
	print(f"AT content: {at_content(example):.2f}%")
	print(f"Reverse complement: {reverse_complement(example)}")
	print(f"mRNA: {transcribe(example)}")
	print(f"Protein: {translate(example)}")

