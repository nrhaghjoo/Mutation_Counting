# nContextMut

Trinucleotide-context mutation analysis tool for multiple sequence alignments (MSA).

---

## What does this pipeline do?

Many mutational processes — such as APOBEC enzyme activity, UV damage, or replication errors — do not mutate nucleotides randomly. Instead, they preferentially target specific nucleotides **depending on the surrounding sequence context**. The most widely studied context is the **trinucleotide context**: the identity of the nucleotide immediately 5' (before) and 3' (after) the mutated position.

For example, APOBEC enzymes preferentially cause `TCA → TTA` and `TCT → TTT` mutations. By counting mutations in their trinucleotide context, you can identify which mutational processes shaped a set of sequences.

**nContextMut** takes a multiple sequence alignment (MSA) of nucleotide sequences and, for each sequence, quantifies how often each trinucleotide context mutates and into what. The pipeline runs in four steps:

---

### Step 1 — Build a consensus sequence

A position-wise consensus is computed from the full alignment. At each column, the most frequent nucleotide is chosen as the consensus base, but only if it appears in at least `threshold` (default 60%) of sequences at that position. This consensus acts as the **reference** against which every individual sequence is compared.

### Step 2 — Count trinucleotides per sequence

Each sequence is globally aligned to the consensus. The consensus is then **trimmed** to the region actually covered by that sequence — removing any leading or trailing positions where the sequence has no coverage (i.e. where the alignment introduces terminal gaps in the query). Triplets are counted on this trimmed, gap-stripped consensus region.

This produces a **trinucleotide frequency table** — the "opportunity" each context had to be mutated — that accurately reflects only the genomic region each sequence spans, avoiding inflation from uncovered regions.

### Step 3 — Identify and count mutations in context

Each sequence is globally aligned to the consensus using pairwise alignment. Positions where the sequence differs from the consensus are identified as mutations. For each mutation, the tool records the **trinucleotide context** in both the reference (consensus) and the query (sequence), producing labels of the form:

```
ATC_AGC   →  middle base mutated from T to G, within context A[T→G]C
```

Mutations adjacent to alignment gaps are excluded to avoid artefacts from indels.

### Step 4 — Normalize mutation counts

Raw mutation counts are divided by the trinucleotide frequency of the reference context from Step 2. This corrects for the fact that some trinucleotide contexts are simply more common in the genome — a context that appears 1000 times and mutates 10 times is very different from one that appears 10 times and mutates 10 times. The result is a **normalized mutation spectrum** comparable across sequences of different lengths and compositions.

---

### Summary of outputs

| Output file | What it contains |
|---|---|
| `_consensus_sequence.txt` | The reference consensus derived from the alignment |
| `_triplet_counts.csv` | Raw trinucleotide frequencies per sequence (the opportunity table) |
| `_mutation_counts.csv` | Raw counts of each reference→mutant trinucleotide context per sequence |
| `_Normalized_counts.csv` | Mutation counts divided by trinucleotide frequency (the mutation spectrum) |

The normalized output is the primary result and can be used directly to compare mutational signatures across samples, time points, or experimental conditions.

---

## Installation

### From source (local)
```bash
pip install .
```

### From a built wheel / PyPI
```bash
pip install ncontextmut
```

## Requirements

- Python ≥ 3.9
- biopython ≥ 1.81
- pandas ≥ 1.5

Dependencies are installed automatically by pip.

---

## Command-line usage

```bash
ncontextmut -i alignment.fasta -o results/ -prefix my_run
```

| Flag | Description |
|---|---|
| `-i / --input` | Path to the input FASTA MSA file |
| `-o / --output_folder` | Directory for output files (created if absent) |
| `-prefix / --output_filename_prefix` | Prefix for every output filename |
| `--threshold` | Consensus frequency threshold (default `0.6`) |

### Output files

| File | Contents |
|---|---|
| `<prefix>_consensus_sequence.txt` | Consensus sequence |
| `<prefix>_triplet_counts.csv` | Raw per-sequence trinucleotide counts |
| `<prefix>_mutation_counts.csv` | Raw per-sequence mutation-context counts |
| `<prefix>_Normalized_counts.csv` | Mutation counts normalised by trinucleotide frequency |

---

## Python API

```python
from ncontextmut import run_analysis, read_fasta, consensus_sequence

# Full pipeline
results = run_analysis(
    fasta_path="alignment.fasta",
    output_folder="results/",
    output_prefix="my_run",
    consensus_threshold=0.6,   # optional
)

print(results["consensus"])
print(results["normalized_counts"].head())

# Lower-level helpers
headers, sequences = read_fasta("alignment.fasta")
cons = consensus_sequence(sequences, threshold=0.7)
```

---

## Input format

A standard FASTA file where **all sequences are already aligned** (same length, gaps represented as `-`).  
RNA sequences (U) are automatically converted to DNA (T).  
Any character that is not `T/C/G/A/U` is replaced with `-`.

---

## License

MIT
