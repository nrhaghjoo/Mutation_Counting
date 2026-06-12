#!/usr/bin/env python3

from collections import Counter
from collections import defaultdict
from Bio.Align import PairwiseAligner
import re
import pandas as pd
import os
import argparse
from pathlib import Path



def all_strings_same_length(strings):
    return len(set(map(len, strings))) == 1

def read_fasta(file_path):
    """Reads a FASTA file and validates all possible errors before processing."""
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Error: The file '{file_path}' does not exist.")

    headers = []
    sequences = []

    with open(file_path, 'r') as file:
        lines = file.readlines()

    if not lines:
        raise ValueError(f"Error: The file '{file_path}' is empty.")

    sequence_data = ''
    for line in lines:
        line = line.strip()

        if line.startswith(">"):
            if sequence_data:
                sequence_data = re.sub(r"[^TCGAU]", "-", sequence_data.upper()).replace("U", "T")
                sequences.append(sequence_data)
                sequence_data = ''

            headers.append(line[1:])  # Store header (without '>')

        elif line:
            sequence_data += line


    if sequence_data:
        sequence_data = re.sub(r"[^TCGAU]", "-", sequence_data.upper()).replace("U", "T")
        sequences.append(sequence_data)

    if not headers:
        raise ValueError(f"Error: The file '{file_path}' does not contain any headers.")

    if len(headers) != len(sequences):
        raise ValueError(f"Error: Mismatch between headers ({len(headers)}) and sequences ({len(sequences)}).")

    if not all_strings_same_length(sequences):
        raise ValueError(f"Error: The file '{file_path}' is not a multiple sequence alignment.")

    return headers, sequences

def consensus_sequence(msa, threshold):
    """Compute consensus sequence from a list of aligned sequences."""
    consensus = []
    for column in zip(*msa):  # Transpose the alignment matrix
        counter = Counter(column)
        del counter['-']
        seq_number = sum(counter.values())
        if len(counter.most_common(1)) != 0:
            ratio = counter.most_common(1)[0][1]/seq_number
            if ratio >= threshold:
                most_common = counter.most_common(1)[0][0]  # Get most frequent character
                consensus.append(most_common)

    return "".join(consensus)

def count_triplets(s):
    triplet_counts = defaultdict(int)
    for i in range(len(s) - 2):  # Ensuring we have enough characters for a triplet
        triplet = s[i:i + 3]  # Extract 3 consecutive characters
        triplet_counts[triplet] += 1

    return triplet_counts

def find_mutations(seq1, seq2):
    mutations = list(filter(lambda i: seq1[i] != seq2[i] and 1 <= i < len(seq1) - 1, range(len(seq1))))
    return (mutations)

# Create argument parser
parser = argparse.ArgumentParser(description="Process some input values.")

# Define command-line arguments
parser.add_argument("-i","--input", required=True, help="Path to the input FASTA file")
parser.add_argument("-o","--output_folder", required=True, help="Path to output folder (absolute or relative)")
parser.add_argument("-prefix","--output_filename_prefix", required=True, help="Prefix for output files")

# Parse arguments
args = parser.parse_args()

fasta_path = os.path.abspath(args.input)


output_folder = args.output_folder
output_folder = os.path.abspath(args.output_folder)
output_filename_prefix = args.output_filename_prefix
# Ensure output folder exists
Path(output_folder).mkdir(parents=True, exist_ok=True)


try:
    viruses , msa = read_fasta(fasta_path)
    consensus = consensus_sequence(msa, 0.6)
    with open(output_folder + "/" + output_filename_prefix + "_" +"consensus_sequence.txt", "w") as file:
        file.write(consensus)


    triplet_counts = defaultdict(lambda: defaultdict(int))
    for i in range(len(viruses)):
        seq = msa[i].replace("-", "")
        for j in range(len(seq) - 2):  # Ensuring we have enough characters for a triplet
            triplet = seq[j:j + 3]  # Extract 3 consecutive characters
            triplet_counts[viruses[i]][triplet] += 1

    # Convert mutation data to a Pandas DataFrame
    triplet_counts = pd.DataFrame.from_dict(triplet_counts, orient="index").fillna(0)
    triplet_counts.to_csv(output_folder + "/" + output_filename_prefix + "_" + "triplet_counts.csv")


    mutation_counts = defaultdict(lambda: defaultdict(int))
    for i in range(len(viruses)):
        # Initialize the aligner
        aligner = PairwiseAligner()
        # Choose alignment mode (global vs local)
        aligner.mode = 'global'  # Use 'local' for local alignment

        # Set scoring parameters
        aligner.match_score = 2  # Score for a match
        aligner.mismatch_score = -1  # Penalty for a mismatch
        aligner.open_gap_score = -5  # Penalty for opening a gap
        aligner.extend_gap_score = -1

        # Define two sequences
        # Perform alignment
        query = msa[i].replace("-", "")
        alignments = aligner.align(consensus, query)

        # Convert alignment to a zipped list format
        best_alignment = alignments[0]  # Get the top alignment
        aligned_seq1, aligned_seq2 = best_alignment[0], best_alignment[1]
        mutation_list = find_mutations(aligned_seq1, aligned_seq2)
        for mutation_point in mutation_list:
            if "-" not in aligned_seq1[(mutation_point-1):(mutation_point+2)] and \
                "-" not in aligned_seq2[(mutation_point-1):(mutation_point+2)]:
                mutation_counts[viruses[i]][aligned_seq1[(mutation_point-1):(mutation_point+2)] +
                "_"+ aligned_seq2[(mutation_point-1):(mutation_point+2)]] += 1

    # Convert mutation data to a Pandas DataFrame
    mutation_counts = pd.DataFrame.from_dict(mutation_counts, orient="index").fillna(0)

    mutation_counts.to_csv(output_folder + "/" + output_filename_prefix + "_" +"mutation_counts.csv")
    print("Process completed successfully.")

    common_samples = mutation_counts.index.intersection(triplet_counts.index)
    mutation_counts = mutation_counts.loc[common_samples]
    triplet_counts = triplet_counts.loc[common_samples]

    result = mutation_counts.copy()

    for col in mutation_counts.columns:
        prefix = col.split("_")[0]  # ATC از ATC_CTG
        result[col] = mutation_counts[col] / triplet_counts[prefix]

    result.to_csv(output_folder + "/" + output_filename_prefix + "_" + "Normalized_counts.csv", index=True)



except Exception as e:
    print(str(e))











