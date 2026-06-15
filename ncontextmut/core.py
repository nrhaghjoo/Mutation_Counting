"""
nContextMut — Trinucleotide-context mutation analysis for multiple sequence alignments.

Author
------
Niloofar Haghjoo
Email  : nhaghjoo@gmail.com
GitHub : https://github.com/nrhaghjoo

Description
-----------
This module implements the core analytical pipeline of nContextMut.
It quantifies nucleotide mutations in their trinucleotide context (the
nucleotide immediately upstream and downstream of the mutated position),
producing normalized mutation spectra that can reveal the underlying
mutational processes active in a set of sequences (e.g. APOBEC activity,
UV damage, polymerase errors).

Pipeline overview
-----------------
1. Read and validate a multiple sequence alignment (MSA) in FASTA format.
2. Derive a consensus sequence from the MSA (used as the reference).
3. Count every overlapping trinucleotide in each sequence (opportunity table).
4. Align each sequence to the consensus and record mutations with their
   trinucleotide context as  <ref_triplet>_<query_triplet>  (e.g. ATC_AGC).
5. Normalise raw mutation counts by the trinucleotide frequency of the
   reference context, yielding a comparable mutation spectrum per sequence.

Public API
----------
read_fasta(file_path)
    Read and validate a FASTA MSA file.

consensus_sequence(msa, threshold=0.6)
    Derive a consensus sequence from an MSA.

count_triplets(s)
    Count overlapping trinucleotides in a sequence string.

find_mutations(seq1, seq2)
    Return positions where two aligned sequences differ.

run_analysis(fasta_path, output_folder, output_prefix, consensus_threshold=0.6)
    Execute the full pipeline and write all output files.

Examples
--------
>>> from ncontextmut import run_analysis
>>> results = run_analysis("alignment.fasta", "results/", "my_run")
>>> print(results["consensus"])
>>> print(results["normalized_counts"].head())

>>> from ncontextmut import read_fasta, consensus_sequence
>>> headers, sequences = read_fasta("alignment.fasta")
>>> cons = consensus_sequence(sequences, threshold=0.7)
"""

import os
import re
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd
from Bio.Align import PairwiseAligner

__author__ = "Niloofar Haghjoo"
__email__ = "nhaghjoo@gmail.com"
__version__ = "1.0.0"


# ---------------------------------------------------------------------------
# Helper utilities
# ---------------------------------------------------------------------------

def all_strings_same_length(strings):
    """
    Check that every string in a collection has the same length.

    Parameters
    ----------
    strings : iterable of str
        The strings to compare.

    Returns
    -------
    bool
        True if all strings share the same length, False otherwise.

    Examples
    --------
    >>> all_strings_same_length(["ATG", "CCC"])
    True
    >>> all_strings_same_length(["ATG", "CC"])
    False
    """
    return len(set(map(len, strings))) == 1


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------

def read_fasta(file_path):
    """
    Read a FASTA file and return headers and aligned sequences.

    The function validates the file thoroughly before returning any data.
    All nucleotide characters are uppercased; U is converted to T (RNA→DNA);
    any character that is not A, T, C, G, or U is replaced with a gap '-'.

    Parameters
    ----------
    file_path : str or pathlib.Path
        Path to the FASTA file. The file must contain at least one sequence,
        and all sequences must have the same length (i.e. it must be a
        multiple sequence alignment).

    Returns
    -------
    headers : list of str
        Sequence identifiers, one per sequence, with the leading '>' stripped.
        Example: ['Sample_1', 'Sample_2', 'Sample_3']

    sequences : list of str
        Aligned nucleotide sequences in DNA alphabet (A/T/C/G/-).
        Each sequence has the same length as every other sequence.
        Example: ['ATCG--AT', 'ATCGCCAT', 'ATC--CAT']

    Raises
    ------
    FileNotFoundError
        If ``file_path`` does not exist on disk.
    ValueError
        If the file is empty, contains no FASTA headers, has a mismatch
        between the number of headers and sequences, or if the sequences
        are not all the same length (not a valid MSA).

    Examples
    --------
    >>> headers, sequences = read_fasta("alignment.fasta")
    >>> print(headers[0])
    'Sample_1'
    >>> print(len(sequences[0]))
    1200
    """
    file_path = str(file_path)
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Error: The file '{file_path}' does not exist.")

    headers = []
    sequences = []

    with open(file_path, "r") as file:
        lines = file.readlines()

    if not lines:
        raise ValueError(f"Error: The file '{file_path}' is empty.")

    sequence_data = ""
    for line in lines:
        line = line.strip()
        if line.startswith(">"):
            if sequence_data:
                sequence_data = re.sub(r"[^TCGAU]", "-", sequence_data.upper()).replace("U", "T")
                sequences.append(sequence_data)
                sequence_data = ""
            headers.append(line[1:])
        elif line:
            sequence_data += line

    if sequence_data:
        sequence_data = re.sub(r"[^TCGAU]", "-", sequence_data.upper()).replace("U", "T")
        sequences.append(sequence_data)

    if not headers:
        raise ValueError(f"Error: The file '{file_path}' does not contain any headers.")

    if len(headers) != len(sequences):
        raise ValueError(
            f"Error: Mismatch between headers ({len(headers)}) and sequences ({len(sequences)})."
        )

    if not all_strings_same_length(sequences):
        raise ValueError(
            f"Error: The file '{file_path}' is not a multiple sequence alignment."
        )

    return headers, sequences


# ---------------------------------------------------------------------------
# Consensus
# ---------------------------------------------------------------------------

def consensus_sequence(msa, threshold=0.6):
    """
    Derive a consensus sequence from a multiple sequence alignment.

    At each alignment column the most frequent non-gap nucleotide is chosen.
    A position is included in the consensus only when that nucleotide's
    frequency among non-gap characters meets or exceeds ``threshold``.
    Columns that are entirely gaps are silently skipped.

    Parameters
    ----------
    msa : list of str
        Multiple sequence alignment. All strings must have the same length.
        Gap characters ('-') are ignored when computing frequencies.

    threshold : float, optional
        Minimum relative frequency (0.0–1.0) required for a nucleotide to
        be accepted as the consensus at a given position. Default is 0.6
        (i.e. the winning nucleotide must appear in at least 60 % of
        non-gap sequences at that column).

    Returns
    -------
    str
        Consensus sequence in DNA alphabet (A/T/C/G). Gap columns and
        positions that do not meet the threshold are omitted, so the
        returned string may be shorter than the alignment width.

    Examples
    --------
    >>> msa = ["ATCG", "ATCG", "ATCG"]
    >>> consensus_sequence(msa)
    'ATCG'

    >>> msa = ["ATCG", "AGCG", "AGCG"]
    >>> consensus_sequence(msa, threshold=0.6)
    'AGCG'   # G wins at position 1 with 2/3 ≈ 67 %

    >>> consensus_sequence(msa, threshold=0.8)
    'ACG'    # position 1 excluded: no nucleotide reaches 80 %
    """
    consensus = []
    for column in zip(*msa):
        counter = Counter(column)
        del counter["-"]
        seq_number = sum(counter.values())
        if counter.most_common(1):
            ratio = counter.most_common(1)[0][1] / seq_number
            if ratio >= threshold:
                consensus.append(counter.most_common(1)[0][0])
    return "".join(consensus)


# ---------------------------------------------------------------------------
# Triplet counting
# ---------------------------------------------------------------------------

def count_triplets(s):
    """
    Count all overlapping trinucleotides (triplets) in a sequence string.

    Every consecutive window of three characters is counted, including
    overlapping windows. The sequence should contain no gap characters;
    strip gaps before calling this function.

    Parameters
    ----------
    s : str
        A nucleotide sequence string with no gap characters ('-').
        Example: 'ATCGATCG'

    Returns
    -------
    collections.defaultdict
        A defaultdict(int) mapping each observed trinucleotide to the
        number of times it appears in ``s``.
        Example: {'ATC': 2, 'TCG': 2, 'CGA': 1, 'GAT': 1}

    Notes
    -----
    A string of length n contains n-2 overlapping trinucleotides.
    Strings shorter than 3 characters return an empty defaultdict.

    Examples
    --------
    >>> count_triplets("ATCG")
    defaultdict(<class 'int'>, {'ATC': 1, 'TCG': 1})

    >>> count_triplets("AAAA")
    defaultdict(<class 'int'>, {'AAA': 2})
    """
    triplet_counts = defaultdict(int)
    for i in range(len(s) - 2):
        triplet = s[i : i + 3]
        triplet_counts[triplet] += 1
    return triplet_counts


# ---------------------------------------------------------------------------
# Mutation detection
# ---------------------------------------------------------------------------

def find_mutations(seq1, seq2):
    """
    Find positions where two aligned sequences differ.

    Only interior positions are considered: the first and last positions
    of the alignment are excluded because they cannot form a complete
    trinucleotide context on both sides.

    Parameters
    ----------
    seq1 : str
        Reference (consensus) aligned sequence.
    seq2 : str
        Query aligned sequence. Must be the same length as ``seq1``.

    Returns
    -------
    list of int
        Zero-based indices of positions where ``seq1[i] != seq2[i]``,
        restricted to the range 1 ≤ i ≤ len(seq1) - 2.
        Returns an empty list if the sequences are identical or differ
        only at the terminal positions.

    Notes
    -----
    This function intentionally excludes the first (index 0) and last
    (index len-1) positions because the surrounding trinucleotide context
    (positions i-1 and i+1) would be incomplete at the sequence ends.

    Examples
    --------
    >>> find_mutations("ATCG", "ATCG")
    []

    >>> find_mutations("ATCG", "AGCG")
    [1]   # position 1: T → G

    >>> find_mutations("ATCG", "AGCA")
    [1]   # position 3 excluded (last position); only position 1 returned
    """
    return [
        i
        for i in range(len(seq1))
        if seq1[i] != seq2[i] and 1 <= i < len(seq1) - 1
    ]


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _build_aligner():
    """
    Create and return a configured Bio.Align.PairwiseAligner.

    The aligner uses global alignment (Needleman-Wunsch) with the
    following scoring scheme:
      match     : +2
      mismatch  : -1
      gap open  : -5
      gap extend: -1

    Returns
    -------
    Bio.Align.PairwiseAligner
    """
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -1
    return aligner

def trim_reference_by_query_gaps(ref_aln: str, qry_aln: str) -> str:
    """
    ref_aln: aligned reference sequence (string)
    qry_aln: aligned query sequence (string)

    Returns: trimmed reference sequence based on leading/trailing gaps in query.
    """

    left_gaps = 0
    for c in qry_aln:
        if c == '-':
            left_gaps += 1
        else:
            break

    right_gaps = 0
    for c in reversed(qry_aln):
        if c == '-':
            right_gaps += 1
        else:
            break

    if right_gaps == 0:
        trimmed = ref_aln[left_gaps:]
    else:
        trimmed = ref_aln[left_gaps: len(ref_aln) - right_gaps]

    return trimmed


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def run_analysis(fasta_path, output_folder, output_prefix, consensus_threshold=0.6):
    """
    Execute the full nContextMut pipeline on a FASTA MSA file.

    This is the main entry point for the analysis. It chains all steps —
    consensus building, trinucleotide counting, mutation detection, and
    normalisation — and writes four output files to ``output_folder``.

    Parameters
    ----------
    fasta_path : str or pathlib.Path
        Path to the input FASTA file. Must be a multiple sequence alignment
        (all sequences the same length). RNA sequences (U) are accepted and
        automatically converted to DNA (T).

    output_folder : str or pathlib.Path
        Directory where all output files will be written. Created
        automatically (including any missing parent directories) if it
        does not already exist.

    output_prefix : str
        String prepended to every output filename.
        Example: 'my_run' → 'my_run_consensus_sequence.txt', etc.

    consensus_threshold : float, optional
        Minimum relative frequency for a nucleotide to be included in the
        consensus sequence. Must be between 0.0 and 1.0. Default is 0.6.
        See ``consensus_sequence()`` for details.

    Returns
    -------
    dict
        A dictionary with four keys:

        'consensus' : str
            The consensus sequence derived from the MSA.

        'triplet_counts' : pandas.DataFrame
            Shape (n_sequences, n_triplets). Each row is a sequence; each
            column is a trinucleotide (e.g. 'ATC'). Values are raw counts
            of how often that trinucleotide appears in the sequence.

        'mutation_counts' : pandas.DataFrame
            Shape (n_sequences, n_mutation_types). Each row is a sequence;
            each column is a mutation context label of the form
            '<ref_triplet>_<query_triplet>'  (e.g. 'ATC_AGC' means the
            middle base changed from T to G in an A[]C context). Values
            are raw counts.

        'normalized_counts' : pandas.DataFrame
            Same shape as ``mutation_counts``. Each mutation count is
            divided by the frequency of the reference trinucleotide context
            in that sequence (from ``triplet_counts``), correcting for
            differences in sequence composition and length.

    Output files
    ------------
    <prefix>_consensus_sequence.txt
        Plain text file containing the consensus sequence.

    <prefix>_triplet_counts.csv
        CSV of raw trinucleotide frequencies (rows = sequences,
        columns = trinucleotides).

    <prefix>_mutation_counts.csv
        CSV of raw mutation-context counts (rows = sequences,
        columns = '<ref>_<query>' labels).

    <prefix>_Normalized_counts.csv
        CSV of normalised mutation frequencies (rows = sequences,
        columns = '<ref>_<query>' labels).

    Raises
    ------
    FileNotFoundError
        If ``fasta_path`` does not exist.
    ValueError
        If the FASTA file is empty, malformed, or not a valid MSA.

    Examples
    --------
    >>> from ncontextmut import run_analysis
    >>> results = run_analysis(
    ...     fasta_path="sequences.fasta",
    ...     output_folder="results/",
    ...     output_prefix="experiment_1",
    ...     consensus_threshold=0.6,
    ... )
    >>> print(results["consensus"])
    'ATCGATCG...'
    >>> results["normalized_counts"].head()
       ATC_AGC  ATC_ACC  ...
    Sample_1   0.002    0.000  ...
    Sample_2   0.000    0.001  ...
    """
    fasta_path = Path(fasta_path).resolve()
    output_folder = Path(output_folder).resolve()
    output_folder.mkdir(parents=True, exist_ok=True)

    viruses, msa = read_fasta(fasta_path)
    consensus = consensus_sequence(msa, consensus_threshold)

    # Write consensus
    consensus_path = output_folder / f"{output_prefix}_consensus_sequence.txt"
    consensus_path.write_text(consensus)


    # --- Mutation counts ---
    aligner = _build_aligner()
    # --- Triplet counts ---
    triplet_counts: dict[str, dict] = defaultdict(lambda: defaultdict(int))
    mutation_counts: dict[str, dict] = defaultdict(lambda: defaultdict(int))

    for virus, seq_aligned in zip(viruses, msa):

        query = seq_aligned.replace("-", "")
        alignments = aligner.align(consensus, query)
        best = alignments[0]
        aligned_ref, aligned_query = best[0], best[1]

        for pos in find_mutations(aligned_ref, aligned_query):
            ref_ctx = aligned_ref[pos - 1 : pos + 2]
            qry_ctx = aligned_query[pos - 1 : pos + 2]
            if "-" not in ref_ctx and "-" not in qry_ctx:
                mutation_counts[virus][f"{ref_ctx}_{qry_ctx}"] += 1



    # --- Triplet counts ---
        trimmed_consensus = trim_reference_by_query_gaps(aligned_ref, aligned_query).replace("-", "")

        for NT in range(len(trimmed_consensus) - 2):  # Ensuring we have enough characters for a triplet
            triplet = trimmed_consensus[NT:NT + 3]  # Extract 3 consecutive characters
            triplet_counts[viruses[virus]][triplet] += 1

    triplet_df = pd.DataFrame.from_dict(triplet_counts, orient="index").fillna(0)
    triplet_df.to_csv(output_folder / f"{output_prefix}_triplet_counts.csv")

    mutation_df = pd.DataFrame.from_dict(mutation_counts, orient="index").fillna(0)
    mutation_df.to_csv(output_folder / f"{output_prefix}_mutation_counts.csv")

    # --- Normalized counts ---
    common = mutation_df.index.intersection(triplet_df.index)
    mut_common = mutation_df.loc[common]
    tri_common = triplet_df.loc[common]

    normalized = mut_common.copy()
    for col in mut_common.columns:
        prefix = col.split("_")[0]
        if prefix in tri_common.columns:
            normalized[col] = mut_common[col] / tri_common[prefix]

    normalized.to_csv(output_folder / f"{output_prefix}_Normalized_counts.csv")

    return {
        "consensus": consensus,
        "triplet_counts": triplet_df,
        "mutation_counts": mutation_df,
        "normalized_counts": normalized,
    }
