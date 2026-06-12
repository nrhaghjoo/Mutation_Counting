# nContextMut

Trinucleotide-context mutation analysis tool for multiple sequence alignments (MSA).

## Installation

### From source (local)
```bash
pip install .
```

### From a built wheel / PyPI (once published)
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
