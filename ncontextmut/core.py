"""Core functions for nContextMut analysis."""

import os
import re
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd
from Bio.Align import PairwiseAligner


def all_strings_same_length(strings):
    """Check that all strings in a list have the same length."""
    return len(set(map(len, strings))) == 1


def read_fasta(file_path):
    """
    Read a FASTA file and validate it before processing.

    Parameters
    ----------
    file_path : str or Path
        Path to the FASTA file.

    Returns
    -------
    headers : list of str
        Sequence headers (without '>').
    sequences : list of str
        Aligned sequences (gaps represented as '-', U replaced by T).

    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    ValueError
        If the file is empty, malformed, or not a valid MSA.
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


def consensus_sequence(msa, threshold=0.6):
    """
    Compute a consensus sequence from a list of aligned sequences.

    Parameters
    ----------
    msa : list of str
        Multiple sequence alignment (all sequences must be the same length).
    threshold : float, optional
        Minimum frequency for a nucleotide to be included in the consensus.
        Default is 0.6.

    Returns
    -------
    str
        Consensus sequence (gaps excluded).
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


def count_triplets(s):
    """
    Count all overlapping triplets in a nucleotide string.

    Parameters
    ----------
    s : str
        Nucleotide sequence (no gaps).

    Returns
    -------
    defaultdict
        Mapping of triplet -> count.
    """
    triplet_counts = defaultdict(int)
    for i in range(len(s) - 2):
        triplet = s[i : i + 3]
        triplet_counts[triplet] += 1
    return triplet_counts


def find_mutations(seq1, seq2):
    """
    Find positions (1-indexed interior) where seq1 and seq2 differ.

    Parameters
    ----------
    seq1, seq2 : str
        Two aligned sequences of equal length.

    Returns
    -------
    list of int
        Positions where the sequences differ, excluding the first and last positions.
    """
    return [
        i
        for i in range(len(seq1))
        if seq1[i] != seq2[i] and 1 <= i < len(seq1) - 1
    ]


def _build_aligner():
    """Return a configured PairwiseAligner."""
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -1
    return aligner


def run_analysis(fasta_path, output_folder, output_prefix, consensus_threshold=0.6):
    """
    Run the full nContextMut pipeline.

    Parameters
    ----------
    fasta_path : str or Path
        Path to the input FASTA MSA file.
    output_folder : str or Path
        Directory where output files will be written (created if needed).
    output_prefix : str
        Prefix for all output file names.
    consensus_threshold : float, optional
        Frequency threshold for consensus calling. Default is 0.6.

    Returns
    -------
    dict with keys:
        'consensus'         – consensus sequence string
        'triplet_counts'    – pd.DataFrame of raw triplet counts
        'mutation_counts'   – pd.DataFrame of raw mutation counts
        'normalized_counts' – pd.DataFrame of normalized mutation frequencies

    Raises
    ------
    FileNotFoundError / ValueError
        Propagated from read_fasta() on bad input.
    """
    fasta_path = Path(fasta_path).resolve()
    output_folder = Path(output_folder).resolve()
    output_folder.mkdir(parents=True, exist_ok=True)

    viruses, msa = read_fasta(fasta_path)
    consensus = consensus_sequence(msa, consensus_threshold)

    # Write consensus
    consensus_path = output_folder / f"{output_prefix}_consensus_sequence.txt"
    consensus_path.write_text(consensus)

    # --- Triplet counts ---
    triplet_counts: dict[str, dict] = defaultdict(lambda: defaultdict(int))
    for virus, seq_aligned in zip(viruses, msa):
        seq = seq_aligned.replace("-", "")
        for j in range(len(seq) - 2):
            triplet = seq[j : j + 3]
            triplet_counts[virus][triplet] += 1

    triplet_df = pd.DataFrame.from_dict(triplet_counts, orient="index").fillna(0)
    triplet_df.to_csv(output_folder / f"{output_prefix}_triplet_counts.csv")

    # --- Mutation counts ---
    aligner = _build_aligner()
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
