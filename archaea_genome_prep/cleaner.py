"""
archaea_genome_prep.cleaner
~~~~~~~~~~~~~~~~~~~~~~~~~~~
Clean and reformat archaeal genome FASTA files for downstream annotation.

Operations
----------
- Remove contigs below minimum length threshold
- Standardise sequence headers to simple numeric IDs
- Replace ambiguous bases with N
- Write cleaned FASTA to output path
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .validator import _parse_fasta_raw, _seq_id, DEFAULT_MIN_CONTIG_LEN


@dataclass
class CleanResult:
    """Result of cleaning a genome FASTA file."""
    input_path: str
    output_path: str
    contigs_input: int
    contigs_output: int
    contigs_removed: int
    total_length_input: int
    total_length_output: int
    headers_renamed: bool

    def summary(self) -> str:
        lines = [
            f"Cleaning complete:",
            f"  Input:   {self.input_path}",
            f"  Output:  {self.output_path}",
            f"  Contigs: {self.contigs_input} -> {self.contigs_output} "
            f"({self.contigs_removed} removed)",
            f"  Length:  {self.total_length_input:,} -> "
            f"{self.total_length_output:,} bp",
        ]
        if self.headers_renamed:
            lines.append("  Headers renamed to contig_1, contig_2, ...")
        return "\n".join(lines)


def clean(
    fasta_path: str | Path,
    output_path: str | Path | None = None,
    min_contig_len: int = DEFAULT_MIN_CONTIG_LEN,
    rename_headers: bool = True,
    prefix: str = "contig",
) -> CleanResult:
    """
    Clean an archaeal genome FASTA file for annotation.

    Parameters
    ----------
    fasta_path : str or Path
        Input FASTA file.
    output_path : str or Path or None
        Output path for cleaned FASTA. If None, writes to
        <input_stem>_cleaned.fasta in the same directory.
    min_contig_len : int
        Remove contigs shorter than this (default 200 bp).
    rename_headers : bool
        Rename sequence headers to simple IDs like contig_1, contig_2, ...
        (default True). Prokka works best with simple headers.
    prefix : str
        Prefix for renamed headers (default "contig").

    Returns
    -------
    CleanResult with statistics about the cleaning operation.
    """
    fasta_path = Path(fasta_path)

    if output_path is None:
        output_path = fasta_path.parent / f"{fasta_path.stem}_cleaned.fasta"
    output_path = Path(output_path)

    kept: list[tuple[str, str]] = []
    total_in = 0
    total_in_len = 0

    for header, seq in _parse_fasta_raw(fasta_path):
        total_in += 1
        total_in_len += len(seq)
        if len(seq) >= min_contig_len:
            kept.append((header, seq.upper()))

    total_out_len = sum(len(s) for _, s in kept)

    with output_path.open("w") as fh:
        for i, (header, seq) in enumerate(kept, 1):
            if rename_headers:
                new_header = f"{prefix}_{i}"
            else:
                new_header = _seq_id(header)
            fh.write(f">{new_header}\n")
            # Write sequence in 60-character lines
            for j in range(0, len(seq), 60):
                fh.write(seq[j:j + 60] + "\n")

    return CleanResult(
        input_path=str(fasta_path),
        output_path=str(output_path),
        contigs_input=total_in,
        contigs_output=len(kept),
        contigs_removed=total_in - len(kept),
        total_length_input=total_in_len,
        total_length_output=total_out_len,
        headers_renamed=rename_headers,
    )
