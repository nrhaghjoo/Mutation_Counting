"""
nContextMut — Trinucleotide-context mutation analysis for multiple sequence alignments.

Author
------
Niloofar Haghjoo
Email  : nhaghjoo@gmail.com
GitHub : https://github.com/nrhaghjoo

Quick start
-----------
>>> from ncontextmut import run_analysis
>>> results = run_analysis("alignment.fasta", "results/", "my_run")

>>> from ncontextmut import read_fasta, consensus_sequence
>>> headers, seqs = read_fasta("alignment.fasta")
>>> cons = consensus_sequence(seqs, threshold=0.6)

Use help() on any function for full documentation:
>>> help(run_analysis)
>>> help(read_fasta)
"""

from .core import (
    all_strings_same_length,
    consensus_sequence,
    count_triplets,
    find_mutations,
    read_fasta,
    run_analysis,
)

__version__ = "1.0.0"
__author__ = "Niloofar Haghjoo"
__email__ = "nhaghjoo@gmail.com"

__all__ = [
    "read_fasta",
    "consensus_sequence",
    "count_triplets",
    "find_mutations",
    "all_strings_same_length",
    "run_analysis",
]
