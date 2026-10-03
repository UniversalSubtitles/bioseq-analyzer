import random

import pytest
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction

import seq_tools

# 100 random DNA sequences (30-300 bp). Seed 42 makes them identical every run.
random.seed(42)
RANDOM_SEQUENCES = [
	"".join(random.choice("ACGT") for _ in range(random.randint(30, 300)))
	for _ in range(100)
]

EXAMPLE = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"


def test_clean_sequence_removes_header_spaces_and_case():
	raw = ">my_gene some description\natg gcc\nATT gta\n"
	assert seq_tools.clean_sequence(raw) == "ATGGCCATTGTA"


def test_invalid_characters_are_detected():
	assert seq_tools.find_invalid_characters("ATGXCCZ") == ["X", "Z"]
	assert seq_tools.find_invalid_characters("ATGC") == []


def test_counts_add_up_to_length():
	for seq in RANDOM_SEQUENCES:
		assert sum(seq_tools.nucleotide_counts(seq).values()) == len(seq)


def test_gc_content_matches_biopython():
	for seq in RANDOM_SEQUENCES:
		assert seq_tools.gc_content(seq) == pytest.approx(gc_fraction(seq) * 100)


def test_reverse_complement_matches_biopython():
	for seq in RANDOM_SEQUENCES:
		assert seq_tools.reverse_complement(seq) == str(Seq(seq).reverse_complement())


def test_transcription_matches_biopython():
	for seq in RANDOM_SEQUENCES:
		assert seq_tools.transcribe(seq) == str(Seq(seq).transcribe())


def test_translation_matches_biopython():
	for seq in RANDOM_SEQUENCES:
		whole_codons = seq[: len(seq) - len(seq) % 3]
		assert seq_tools.translate(seq) == str(Seq(whole_codons).translate())


def test_known_example():
	assert seq_tools.translate(EXAMPLE) == "MAIVMGR*KGAR*"
	assert round(seq_tools.gc_content(EXAMPLE), 2) == 56.41


def test_parse_fasta_handles_multiple_records_and_plain_text():
    multi = ">gene1 first\nATG\nCCC\n>gene2\nttt gga\n"
    assert seq_tools.parse_fasta(multi) == [("gene1 first", "ATGCCC"), ("gene2", "TTTGGA")]
    assert seq_tools.parse_fasta("atg cc\nGG") == [("Unnamed sequence", "ATGCCGG")]
    assert seq_tools.parse_fasta("   \n") == []

# ---------- ORF finder ----------
# Longer random sequences so that ORFs actually occur
random.seed(7)
LONG_SEQUENCES = [
    "".join(random.choice("ACGT") for _ in range(random.randint(300, 1500)))
    for _ in range(50)
]


def biopython_orfs(seq, min_aa):
    """Independent ORF finder: translate each frame with Biopython, split at stops."""
    found = set()
    for strand, nucleotides in (("+", Seq(seq)), ("-", Seq(seq).reverse_complement())):
        for frame in range(3):
            usable = 3 * ((len(seq) - frame) // 3)
            protein = str(nucleotides[frame:frame + usable].translate())
            for piece in protein.split("*")[:-1]:  # the last piece has no stop codon after it
                if "M" in piece:
                    orf_protein = piece[piece.index("M"):]
                    if len(orf_protein) >= min_aa:
                        found.add((f"{strand}{frame + 1}", orf_protein))
    return found


def test_orfs_match_biopython():
    for seq in LONG_SEQUENCES:
        ours = {(orf["Frame"], orf["Protein"]) for orf in seq_tools.find_orfs(seq, min_protein_length=10)}
        assert ours == biopython_orfs(seq, 10)


def test_orf_coordinates_point_to_start_and_stop_codons():
    for seq in LONG_SEQUENCES:
        for orf in seq_tools.find_orfs(seq, min_protein_length=10):
            if orf["Frame"].startswith("+"):
                region = seq[orf["Start"] - 1:orf["End"]]
            else:
                region = seq_tools.reverse_complement(seq[orf["End"] - 1:orf["Start"]])
            assert region.startswith("ATG")
            assert region[-3:] in seq_tools.STOP_CODONS
            assert len(region) == orf["Length (nt)"] == 3 * orf["Length (aa)"] + 3


def test_simple_orf_on_both_strands():
    forward = "CC" + "ATG" + "GCC" * 10 + "TAA" + "GG"
    orfs = seq_tools.find_orfs(forward, min_protein_length=5)
    assert orfs[0]["Frame"] == "+3"
    assert (orfs[0]["Start"], orfs[0]["End"]) == (3, 38)
    reverse = seq_tools.reverse_complement(forward)
    orfs = seq_tools.find_orfs(reverse, min_protein_length=5)
    assert orfs[0]["Frame"].startswith("-")
    assert (orfs[0]["Start"], orfs[0]["End"]) == (38, 3)


# ---------- Restriction sites ----------
def test_recognition_sites_are_palindromic():
    for site in seq_tools.RESTRICTION_ENZYMES.values():
        assert seq_tools.reverse_complement(site) == site


def test_restriction_sites_match_biopython():
    from Bio import Restriction

    for seq in LONG_SEQUENCES:
        ours = seq_tools.find_restriction_sites(seq)
        for name, positions in ours.items():
            enzyme = getattr(Restriction, name)
            # Biopython reports cut positions; convert our site starts the same way
            expected = sorted(p + enzyme.fst5 for p in positions)
            assert sorted(enzyme.search(Seq(seq))) == expected


def test_overlapping_sites_are_all_found():
    assert seq_tools.find_restriction_sites("GAATTCGAATTC")["EcoRI"] == [1, 7]
    assert seq_tools.find_restriction_sites("CCCGGGGCCCGGG")["SmaI"] == [1, 8]