"""
archaea-genome-prep
~~~~~~~~~~~~~~~~~~~
Validate, clean, and annotate archaeal genome FASTA files for
downstream bioinformatic analysis.

Quick start
-----------
>>> from archaea_genome_prep import prepare
>>> result = prepare("my_genome.fasta", genus="Methanosarcina", species="acetivorans")
>>> print(result.gbk_path)  # path to annotated GenBank file

Individual steps
----------------
>>> from archaea_genome_prep import validate, clean, annotate
>>> validation = validate("my_genome.fasta")
>>> print(validation.summary())
"""

from .validator import validate, ValidationResult, ContigInfo
from .cleaner import clean, CleanResult
from .annotator import annotate, AnnotationResult, prokka_available
from .pipeline import prepare, PipelineResult

__version__ = "0.1.0"
__all__ = [
    "validate", "ValidationResult", "ContigInfo",
    "clean", "CleanResult",
    "annotate", "AnnotationResult", "prokka_available",
    "prepare", "PipelineResult",
]
