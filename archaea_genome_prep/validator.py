"""
archaea_genome_prep.validator
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Validate and clean archaeal genome FASTA files before annotation.

Checks performed
----------------
1. File is valid FASTA format
2. All sequences contain only IUPAC nucleotide characters
3. No duplicate sequence headers
4. Contigs meet minimum length threshold
5. Total assembly size is within expected range for archaea
6. GC content is within expected range for methanogens (30-65%)
7. No extremely short contigs that would confuse annotators

References
----------
Methanogen genome sizes and GC content ranges from NCBI RefSeq:
  Methanosarcina acetivorans C2A: 5.75 Mb, 42.7% GC
  Methanosarcina mazei Go1: 4.10 Mb, 41.5% GC
  Methanosarcina barkeri: 4.84 Mb, 38.3% GC
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator


# IUPAC nucleotide characters (includes ambiguity codes)
_VALID_BASES = set("ACGTacgtNnRrYySsWwKkMmBbDdHhVv")

# Minimum contig length to retain
DEFAULT_MIN_CONTIG_LEN = 200

# Expected archaeal genome size range (bp)
_GENOME_SIZE_MIN = 1_000_000
_GENOME_SIZE_MAX = 10_000_000

# Expected GC range for methanogens (fraction)
_GC_MIN = 0.28
_GC_MAX = 0.68


@dataclass
class ContigInfo:
    """Information about a single contig in the assembly."""
    header: str
    seq_id: str
    length: int
    gc_fraction: float
    has_invalid_bases: bool
    invalid_base_count: int = 0


@dataclass
class ValidationResult:
    """Result of validating a genome FASTA file."""
    path: str
    valid: bool
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    contigs: list[ContigInfo] = field(default_factory=list)

    total_length: int = 0
    contig_count: int = 0
    gc_overall: float = 0.0
    short_contigs_removed: int = 0

    @property
    def n50(self) -> int:
        """Compute N50 of the assembly."""
        lengths = sorted([c.length for c in self.contigs], reverse=True)
        total = sum(lengths)
        cumsum = 0
        for ln in lengths:
            cumsum += ln
            if cumsum >= total / 2:
                return ln
        return 0

    def summary(self) -> str:
        status = "VALID" if self.valid else "INVALID"
        lines = [
            f"Validation result: {status}",
            f"  File: {self.path}",
            f"  Contigs: {self.contig_count}",
            f"  Total length: {self.total_length:,} bp",
            f"  N50: {self.n50:,} bp",
            f"  GC content: {self.gc_overall:.1%}",
        ]
        if self.short_contigs_removed:
            lines.append(
                f"  Short contigs removed: {self.short_contigs_removed}"
            )
        if self.errors:
            lines.append("  Errors:")
            for e in self.errors:
                lines.append(f"    ERROR: {e}")
        if self.warnings:
            lines.append("  Warnings:")
            for w in self.warnings:
                lines.append(f"    WARNING: {w}")
        return "\n".join(lines)


def _parse_fasta_raw(path: Path) -> Iterator[tuple[str, str]]:
    """Yield (header, sequence) pairs from a FASTA file."""
    current_header: str | None = None
    buffer: list[str] = []

    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_header is not None:
                    yield current_header, "".join(buffer)
                current_header = line[1:]
                buffer = []
            else:
                buffer.append(line)

    if current_header is not None:
        yield current_header, "".join(buffer)


def _gc_fraction(seq: str) -> float:
    seq = seq.upper()
    if not seq:
        return 0.0
    gc = seq.count("G") + seq.count("C")
    return gc / len(seq)


def _seq_id(header: str) -> str:
    return header.split()[0]


def validate(
    fasta_path: str | Path,
    min_contig_len: int = DEFAULT_MIN_CONTIG_LEN,
) -> ValidationResult:
    """
    Validate an archaeal genome FASTA file.

    Parameters
    ----------
    fasta_path : str or Path
        Path to the FASTA file to validate.
    min_contig_len : int
        Minimum contig length to retain (default 200 bp).
        Shorter contigs are counted but flagged as warnings.

    Returns
    -------
    ValidationResult with validation status, errors, warnings,
    and per-contig statistics.
    """
    path = Path(fasta_path)
    result = ValidationResult(path=str(path), valid=True)

    # Check file exists
    if not path.exists():
        result.valid = False
        result.errors.append(f"File not found: {path}")
        return result

    # Check file is not empty
    if path.stat().st_size == 0:
        result.valid = False
        result.errors.append("File is empty")
        return result

    # Parse and validate
    seen_ids: set[str] = set()
    total_gc = 0
    total_len = 0
    short_count = 0
    has_sequences = False

    try:
        for header, seq in _parse_fasta_raw(path):
            has_sequences = True
            seq_id = _seq_id(header)

            # Duplicate header check
            if seq_id in seen_ids:
                result.valid = False
                result.errors.append(f"Duplicate sequence ID: '{seq_id}'")
            seen_ids.add(seq_id)

            # Empty sequence check
            if not seq:
                result.valid = False
                result.errors.append(f"Empty sequence for '{seq_id}'")
                continue

            # Invalid base check
            invalid = set(seq) - _VALID_BASES
            has_invalid = bool(invalid)
            invalid_count = sum(1 for b in seq if b not in _VALID_BASES)
            if has_invalid:
                result.valid = False
                result.errors.append(
                    f"Invalid characters in '{seq_id}': "
                    f"{invalid} ({invalid_count} bases)"
                )

            gc = _gc_fraction(seq)
            total_gc += seq.upper().count("G") + seq.upper().count("C")
            total_len += len(seq)

            if len(seq) < min_contig_len:
                short_count += 1

            result.contigs.append(ContigInfo(
                header=header,
                seq_id=seq_id,
                length=len(seq),
                gc_fraction=gc,
                has_invalid_bases=has_invalid,
                invalid_base_count=invalid_count,
            ))

    except UnicodeDecodeError:
        result.valid = False
        result.errors.append("File is not valid UTF-8 text")
        return result

    if not has_sequences:
        result.valid = False
        result.errors.append("No FASTA sequences found in file")
        return result

    # Assembly-level stats
    result.contig_count = len(result.contigs)
    result.total_length = total_len
    result.short_contigs_removed = short_count
    result.gc_overall = total_gc / total_len if total_len > 0 else 0.0

    # Assembly-level warnings
    if total_len < _GENOME_SIZE_MIN:
        result.warnings.append(
            f"Total assembly size ({total_len:,} bp) is very small for an "
            f"archaeal genome. Expected >= {_GENOME_SIZE_MIN:,} bp. "
            f"Check that this is a complete genome, not a partial assembly."
        )
    elif total_len > _GENOME_SIZE_MAX:
        result.warnings.append(
            f"Total assembly size ({total_len:,} bp) is very large. "
            f"Expected <= {_GENOME_SIZE_MAX:,} bp for most methanogens."
        )

    if result.gc_overall < _GC_MIN:
        result.warnings.append(
            f"Overall GC content ({result.gc_overall:.1%}) is lower than "
            f"expected for methanogens ({_GC_MIN:.0%}-{_GC_MAX:.0%}). "
            f"Verify this is the correct organism."
        )
    elif result.gc_overall > _GC_MAX:
        result.warnings.append(
            f"Overall GC content ({result.gc_overall:.1%}) is higher than "
            f"expected for methanogens ({_GC_MIN:.0%}-{_GC_MAX:.0%}). "
            f"Verify this is the correct organism."
        )

    if short_count > 0:
        result.warnings.append(
            f"{short_count} contig(s) shorter than {min_contig_len} bp. "
            f"These will be excluded from annotation by Prokka. "
            f"Use clean() to remove them."
        )

    n_contigs = len(result.contigs)
    if n_contigs > 500:
        result.warnings.append(
            f"Assembly is highly fragmented ({n_contigs} contigs). "
            f"Annotation quality may be reduced. Consider scaffolding first."
        )

    return result
