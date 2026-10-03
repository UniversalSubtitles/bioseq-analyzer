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