# archaea-genome-prep

[![CI](https://github.com/CameronPiepkorn/archaea-genome-prep/actions/workflows/ci.yml/badge.svg)](https://github.com/CameronPiepkorn/archaea-genome-prep/actions)
[![PyPI](https://img.shields.io/pypi/v/archaea-genome-prep.svg)](https://pypi.org/project/archaea-genome-prep/)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

A Python library for validating, cleaning, and preparing archaeal genome FASTA files for downstream annotation and bioinformatic analysis.

Designed as a companion to archaeal genome analysis tools that require annotated GenBank files as input, such as att site scanners and prophage detection pipelines.

---

## What This Does

1. **Validates** your genome FASTA — checks for duplicate headers, invalid bases, assembly size, GC content, and contig length
2. **Cleans** the FASTA — removes short contigs, standardises headers for Prokka compatibility
3. **Annotates** with Prokka — runs Prokka with archaeal-appropriate settings (`--kingdom Archaea`) if it is installed, or provides the exact command to run manually

---

## Why This Exists

Most archaeal genome analysis tools require annotated GenBank (.gbk) files, not raw FASTA files. Getting from a downloaded FASTA to an annotated GBK requires:

- Validating the FASTA is correctly formatted
- Removing short contigs that confuse annotators
- Running Prokka with the `--kingdom Archaea` flag (the default bacterial setting produces poor results for archaea)

This library automates those steps with sensible defaults tuned for methanogenic archaea.

---

## Installation

```bash
pip install archaea-genome-prep
```

Prokka is a soft dependency — the library works without it (validation and cleaning still run), but annotation requires Prokka to be installed:

```bash
conda install -c bioconda prokka
```

---

## Quick Start

```python
from archaea_genome_prep import prepare

# Full pipeline: validate -> clean -> annotate
result = prepare(
    "my_genome.fasta",
    genus="Methanosarcina",
    species="acetivorans",
    strain="C2A",
)
print(result.summary())

# Path to annotated GenBank file (if Prokka ran successfully)
print(result.gbk_path)
```

---

## Individual Steps

### Validate only

```python
from archaea_genome_prep import validate

result = validate("my_genome.fasta")
print(result.summary())

# Check specific properties
print(f"Contigs: {result.contig_count}")
print(f"Total length: {result.total_length:,} bp")
print(f"N50: {result.n50:,} bp")
print(f"GC content: {result.gc_overall:.1%}")
print(f"Valid: {result.valid}")
print(f"Errors: {result.errors}")
print(f"Warnings: {result.warnings}")
```

### Clean only

```python
from archaea_genome_prep import clean

result = clean(
    "my_genome.fasta",
    output_path="my_genome_cleaned.fasta",
    min_contig_len=200,
    rename_headers=True,
)
print(result.summary())
```

### Annotate only (assumes cleaned FASTA)

```python
from archaea_genome_prep import annotate, prokka_available

if prokka_available():
    result = annotate(
        "my_genome_cleaned.fasta",
        output_dir="annotation/",
        prefix="my_genome",
        genus="Methanosarcina",
        species="acetivorans",
        cpus=8,
    )
    print(result.summary())
    print(result.gbk_path)
else:
    print("Prokka not found. Install with: conda install -c bioconda prokka")
```

### Check the Prokka command without running it

```python
result = annotate("my_genome_cleaned.fasta", run=False)
print(result.command_run)  # prints the exact Prokka command
```

---

## Validation Checks

| Check | Description |
|---|---|
| File exists and is non-empty | Basic file sanity |
| Valid FASTA format | Sequences follow > header / sequence structure |
| IUPAC nucleotide characters only | Catches binary files or protein sequences |
| No duplicate sequence IDs | Required for all downstream tools |
| Contig length >= minimum | Short contigs cause annotation errors |
| Assembly size 1-10 Mb | Expected range for methanogenic archaea |
| GC content 28-68% | Expected range for methanogens |

---

## Note on Pyrrolysine-Encoding Methanogens

*Methanosarcina acetivorans*, *M. mazei*, and *M. barkeri* encode pyrrolysine at UAG codons. Prokka does not handle this ambiguity — genes containing pyrrolysine UAG codons will be annotated as truncated.

After annotation, use [amber-codon-scanner](https://github.com/CameronPiepkorn/amber-codon-scanner) to identify pyrrolysine-encoding UAG codons in the annotated genome.

---

## Reference

Prokka annotation:

> Seemann T (2014) Prokka: rapid prokaryotic genome annotation.
> Bioinformatics 30(14):2068-2069. doi:10.1093/bioinformatics/btu153

---

## License

MIT. See [LICENSE](LICENSE).
