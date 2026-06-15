"""Basic tests for ncontextmut."""

# import pytest
from ncontextmut import all_strings_same_length, consensus_sequence, count_triplets, find_mutations


def test_all_strings_same_length():
    assert all_strings_same_length(["ATG", "CCC"]) is True
    assert all_strings_same_length(["ATG", "CC"]) is False


def test_consensus_sequence_simple():
    msa = ["ATCG", "ATCG", "ATCG"]
    assert consensus_sequence(msa, threshold=0.6) == "ATCG"


def test_consensus_sequence_with_gaps():
    msa = ["A-CG", "ATCG", "ATCG"]
    cons = consensus_sequence(msa, threshold=0.6)
    assert "T" in cons  # T wins at position 1


def test_count_triplets():
    counts = count_triplets("ATCG")
    assert counts["ATC"] == 1
    assert counts["TCG"] == 1


def test_find_mutations_no_diff():
    assert find_mutations("ATCG", "ATCG") == []


def test_find_mutations_middle():
    # Position 1 (middle) differs: A[T]CG vs A[G]CG
    result = find_mutations("ATCG", "AGCG")
    assert 1 in result
