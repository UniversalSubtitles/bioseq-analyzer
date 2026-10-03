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


if __name__ == "__main__":
	example = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
	print(f"Length: {len(example)} bp")
	print(f"Counts: {nucleotide_counts(example)}")
	print(f"GC content: {gc_content(example):.2f}%")
	print(f"AT content: {at_content(example):.2f}%")
	print(f"Reverse complement: {reverse_complement(example)}")
	print(f"mRNA: {transcribe(example)}")
	print(f"Protein: {translate(example)}")
