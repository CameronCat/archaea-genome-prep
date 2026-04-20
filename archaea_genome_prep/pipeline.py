"""
archaea_genome_prep.pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
End-to-end pipeline: validate -> clean -> annotate.

This is the main entry point for most users.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .validator import validate, ValidationResult, DEFAULT_MIN_CONTIG_LEN
from .cleaner import clean, CleanResult
from .annotator import annotate, AnnotationResult, prokka_available


@dataclass
class PipelineResult:
    """Combined result of the full prepare pipeline."""
    validation: ValidationResult
    cleaning: CleanResult | None
    annotation: AnnotationResult | None

    @property
    def ready_for_annotation(self) -> bool:
        return self.cleaning is not None

    @property
    def gbk_path(self) -> str | None:
        if self.annotation:
            return self.annotation.gbk_path
        return None

    def summary(self) -> str:
        lines = ["=" * 60, "  Archaeal Genome Prep Pipeline", "=" * 60, ""]
        lines.append(self.validation.summary())

        if self.cleaning:
            lines += ["", self.cleaning.summary()]

        if self.annotation:
            lines += ["", self.annotation.summary()]
        elif self.cleaning and not prokka_available():
            lines += [
                "",
                "Prokka not found. To annotate your cleaned genome:",
                "  conda install -c bioconda prokka",
                f"  prokka --kingdom Archaea --outdir annotation/ "
                f"{self.cleaning.output_path}",
            ]

        lines += ["", "=" * 60]
        return "\n".join(lines)


def prepare(
    fasta_path: str | Path,
    output_dir: str | Path | None = None,
    min_contig_len: int = DEFAULT_MIN_CONTIG_LEN,
    rename_headers: bool = True,
    prefix: str | None = None,
    genus: str = "",
    species: str = "",
    strain: str = "",
    cpus: int = 4,
    run_annotation: bool = True,
    stop_on_error: bool = True,
) -> PipelineResult:
    """
    Prepare an archaeal genome FASTA file for annotation.

    This function runs the full pipeline:
    1. Validate the FASTA file
    2. Clean it (remove short contigs, standardise headers)
    3. Annotate with Prokka if available

    Parameters
    ----------
    fasta_path : str or Path
        Raw genome FASTA file.
    output_dir : str or Path or None
        Directory for all output files. Defaults to same directory
        as the input FASTA.
    min_contig_len : int
        Remove contigs shorter than this (default 200 bp).
    rename_headers : bool
        Rename sequence headers to simple IDs (default True).
    prefix : str or None
        Prefix for output files. Defaults to FASTA stem.
    genus : str
        Genus name for Prokka (optional).
    species : str
        Species name for Prokka (optional).
    strain : str
        Strain name for Prokka (optional).
    cpus : int
        CPUs for Prokka (default 4).
    run_annotation : bool
        Run Prokka if available (default True).
        Set to False to only validate and clean.
    stop_on_error : bool
        If True (default), stop pipeline if validation finds errors.
        If False, attempt to clean and annotate anyway.

    Returns
    -------
    PipelineResult with results from each stage.
    """
    fasta_path = Path(fasta_path)

    if output_dir is None:
        output_dir = fasta_path.parent
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if prefix is None:
        prefix = fasta_path.stem

    # Stage 1 - Validate
    print(f"[1/3] Validating {fasta_path.name}...")
    validation = validate(fasta_path, min_contig_len=min_contig_len)
    print(validation.summary())

    if not validation.valid and stop_on_error:
        print("\nValidation failed. Fix errors before proceeding.")
        return PipelineResult(
            validation=validation,
            cleaning=None,
            annotation=None,
        )

    # Stage 2 - Clean
    print(f"\n[2/3] Cleaning {fasta_path.name}...")
    cleaned_path = output_dir / f"{prefix}_cleaned.fasta"
    cleaning = clean(
        fasta_path,
        output_path=cleaned_path,
        min_contig_len=min_contig_len,
        rename_headers=rename_headers,
        prefix="contig",
    )
    print(cleaning.summary())

    # Stage 3 - Annotate
    annotation = None
    if run_annotation:
        print(f"\n[3/3] Annotating with Prokka...")
        if not prokka_available():
            print("  Prokka not found on PATH.")
            print("  Install with: conda install -c bioconda prokka")
            print(f"  Then run: prokka --kingdom Archaea --outdir "
                  f"{output_dir / 'annotation'} {cleaned_path}")
        else:
            annotation = annotate(
                cleaned_path,
                output_dir=output_dir / "annotation",
                prefix=prefix,
                genus=genus,
                species=species,
                strain=strain,
                cpus=cpus,
            )
            print(annotation.summary())
    else:
        print("\n[3/3] Skipping annotation (run_annotation=False)")

    return PipelineResult(
        validation=validation,
        cleaning=cleaning,
        annotation=annotation,
    )
