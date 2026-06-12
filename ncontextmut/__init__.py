"""
nContextMut: N-context mutation analysis tool for multiple sequence alignments.
"""

from .core import (
    read_fasta,
    consensus_sequence,
    count_triplets,
    find_mutations,
    all_strings_same_length,
    run_analysis,
)

__version__ = "1.0.0"
__author__ = "nContextMut Authors"
__all__ = [
    "read_fasta",
    "consensus_sequence",
    "count_triplets",
    "find_mutations",
    "all_strings_same_length",
    "run_analysis",
]
