"""
archaea_genome_prep.annotator
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Wrapper around Prokka for archaeal genome annotation.

This module detects whether Prokka is installed and provides a
Python interface to run it with archaeal-appropriate settings.

Prokka reference
----------------
Seemann T (2014) Prokka: rapid prokaryotic genome annotation.
Bioinformatics 30(14):2068-2069. doi:10.1093/bioinformatics/btu153

The critical setting for archaea is --kingdom Archaea, which:
- Uses archaeal gene models for ORF calling
- Uses archaeal ribosomal RNA databases
- Applies archaeal start codon rules (GTG, TTG more common)
- Uses the correct genetic code for most archaea

Note on pyrrolysine-encoding methanogens
-----------------------------------------
Methanosarcina species encoding pyrrolysine (M. acetivorans, M. mazei,
M. barkeri) use UAG codons for both pyrrolysine and stop. Prokka does
not handle this ambiguity -- UAG-encoded pyrrolysine genes will be
annotated as truncated. Use amber-codon-scanner after annotation to
identify these sites.
"""

from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path


@dataclass
class AnnotationResult:
    """Result of running Prokka annotation."""
    input_fasta: str
    output_dir: str
    prefix: str
    gbk_path: str | None
    gff_path: str | None
    success: bool
    prokka_available: bool
    command_run: str
    stdout: str = ""
    stderr: str = ""
    error_message: str = ""

    def summary(self) -> str:
        lines = []
        if not self.prokka_available:
            lines += [
                "Prokka is not installed or not on PATH.",
                "",
                "To install Prokka:",
                "  conda install -c bioconda prokka",
                "",
                "Then run annotation with:",
                f"  prokka --kingdom Archaea --outdir {self.output_dir} "
                f"--prefix {self.prefix} {self.input_fasta}",
                "",
                "The --kingdom Archaea flag is required for correct gene models.",
            ]
            return "\n".join(lines)

        if self.success:
            lines += [
                "Annotation complete:",
                f"  Input FASTA: {self.input_fasta}",
                f"  Output dir:  {self.output_dir}",
                f"  GenBank:     {self.gbk_path}",
                f"  GFF:         {self.gff_path}",
                "",
                "Next steps:",
                "  - Review annotations in a genome browser (Artemis, IGV)",
                "  - The .gbk file can be used directly with archaeal-att-scanner",
            ]
        else:
            lines += [
                "Annotation failed:",
                f"  Error: {self.error_message}",
                f"  Command: {self.command_run}",
                "  Check Prokka is correctly installed and the input FASTA is valid.",
            ]
        return "\n".join(lines)


def prokka_available() -> bool:
    """Return True if Prokka is installed and on PATH."""
    return shutil.which("prokka") is not None


def annotate(
    fasta_path: str | Path,
    output_dir: str | Path | None = None,
    prefix: str | None = None,
    kingdom: str = "Archaea",
    genus: str = "",
    species: str = "",
    strain: str = "",
    cpus: int = 4,
    run: bool = True,
) -> AnnotationResult:
    """
    Annotate an archaeal genome FASTA file using Prokka.

    Parameters
    ----------
    fasta_path : str or Path
        Path to the cleaned genome FASTA file.
    output_dir : str or Path or None
        Directory for Prokka output. Defaults to <fasta_stem>_prokka/
    prefix : str or None
        Prefix for output files. Defaults to FASTA stem.
    kingdom : str
        Kingdom to use for gene models (default "Archaea").
        Do not change this unless you know what you are doing.
    genus : str
        Genus name (optional, improves annotation quality).
    species : str
        Species name (optional).
    strain : str
        Strain name (optional).
    cpus : int
        Number of CPUs for Prokka (default 4).
    run : bool
        If True (default), run Prokka. If False, return the command
        that would be run without executing it. Useful for checking
        the command before running on a cluster.

    Returns
    -------
    AnnotationResult with paths to output files and run status.
    """
    fasta_path = Path(fasta_path)

    if output_dir is None:
        output_dir = fasta_path.parent / f"{fasta_path.stem}_prokka"
    output_dir = Path(output_dir)

    if prefix is None:
        prefix = fasta_path.stem

    # Build Prokka command
    cmd = [
        "prokka",
        "--kingdom", kingdom,
        "--outdir", str(output_dir),
        "--prefix", prefix,
        "--cpus", str(cpus),
        "--force",
    ]

    if genus:
        cmd += ["--genus", genus]
    if species:
        cmd += ["--species", species]
    if strain:
        cmd += ["--strain", strain]

    cmd.append(str(fasta_path))
    cmd_str = " ".join(cmd)

    # Check Prokka availability
    if not prokka_available():
        return AnnotationResult(
            input_fasta=str(fasta_path),
            output_dir=str(output_dir),
            prefix=prefix,
            gbk_path=None,
            gff_path=None,
            success=False,
            prokka_available=False,
            command_run=cmd_str,
            error_message="Prokka not found on PATH",
        )

    if not run:
        return AnnotationResult(
            input_fasta=str(fasta_path),
            output_dir=str(output_dir),
            prefix=prefix,
            gbk_path=str(output_dir / f"{prefix}.gbk"),
            gff_path=str(output_dir / f"{prefix}.gff"),
            success=True,
            prokka_available=True,
            command_run=cmd_str,
        )

    # Run Prokka
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True,
        )

        gbk_path = output_dir / f"{prefix}.gbk"
        gff_path = output_dir / f"{prefix}.gff"

        return AnnotationResult(
            input_fasta=str(fasta_path),
            output_dir=str(output_dir),
            prefix=prefix,
            gbk_path=str(gbk_path) if gbk_path.exists() else None,
            gff_path=str(gff_path) if gff_path.exists() else None,
            success=True,
            prokka_available=True,
            command_run=cmd_str,
            stdout=result.stdout,
            stderr=result.stderr,
        )

    except subprocess.CalledProcessError as e:
        return AnnotationResult(
            input_fasta=str(fasta_path),
            output_dir=str(output_dir),
            prefix=prefix,
            gbk_path=None,
            gff_path=None,
            success=False,
            prokka_available=True,
            command_run=cmd_str,
            stdout=e.stdout or "",
            stderr=e.stderr or "",
            error_message=f"Prokka exited with code {e.returncode}",
        )
