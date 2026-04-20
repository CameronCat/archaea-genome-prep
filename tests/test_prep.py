"""
Tests for archaea-genome-prep.
"""

import pytest
import tempfile
from pathlib import Path
from archaea_genome_prep import validate, clean, prokka_available
from archaea_genome_prep.validator import ValidationResult


# ── Fixtures ──────────────────────────────────────────────────────────────

def write_fasta(tmp_path: Path, records: list[tuple[str, str]]) -> Path:
    p = tmp_path / "test.fasta"
    with p.open("w") as f:
        for header, seq in records:
            f.write(f">{header}\n{seq}\n")
    return p


VALID_SEQ = "ATGCATGCATGCATGCATGC" * 50   # 1000 bp, ~50% GC
LARGE_SEQ = "ATGCATGCATGCATGCATGC" * 100000  # 2 Mb


# ── Validator tests ────────────────────────────────────────────────────────

class TestValidator:

    def test_valid_single_contig(self, tmp_path):
        fa = write_fasta(tmp_path, [("seq1", LARGE_SEQ)])
        result = validate(fa)
        assert isinstance(result, ValidationResult)
        assert len(result.contigs) == 1
        assert result.contig_count == 1

    def test_valid_multiple_contigs(self, tmp_path):
        fa = write_fasta(tmp_path, [
            ("seq1", LARGE_SEQ),
            ("seq2", LARGE_SEQ),
        ])
        result = validate(fa)
        assert result.contig_count == 2

    def test_duplicate_header_is_error(self, tmp_path):
        fa = write_fasta(tmp_path, [
            ("seq1", VALID_SEQ),
            ("seq1", VALID_SEQ),
        ])
        result = validate(fa)
        assert not result.valid
        assert any("Duplicate" in e for e in result.errors)

    def test_invalid_bases_detected(self, tmp_path):
        fa = write_fasta(tmp_path, [("seq1", "ATGCXYZATGC" * 200)])
        result = validate(fa)
        assert not result.valid
        assert any("Invalid characters" in e for e in result.errors)

    def test_empty_file_is_error(self, tmp_path):
        p = tmp_path / "empty.fasta"
        p.write_text("")
        result = validate(p)
        assert not result.valid

    def test_missing_file_is_error(self, tmp_path):
        result = validate(tmp_path / "nonexistent.fasta")
        assert not result.valid
        assert any("not found" in e for e in result.errors)

    def test_short_contig_warning(self, tmp_path):
        fa = write_fasta(tmp_path, [
            ("seq1", LARGE_SEQ),
            ("seq2", "ATGC" * 10),  # 40 bp - too short
        ])
        result = validate(fa, min_contig_len=200)
        assert any("shorter" in w for w in result.warnings)

    def test_gc_content_calculated(self, tmp_path):
        all_gc = "GCGCGCGCGC" * 200000  # 100% GC
        fa = write_fasta(tmp_path, [("seq1", all_gc)])
        result = validate(fa)
        assert result.gc_overall == pytest.approx(1.0)

    def test_n50_calculated(self, tmp_path):
        fa = write_fasta(tmp_path, [
            ("seq1", LARGE_SEQ),
            ("seq2", LARGE_SEQ[:500_000]),
        ])
        result = validate(fa)
        assert result.n50 > 0

    def test_summary_contains_status(self, tmp_path):
        fa = write_fasta(tmp_path, [("seq1", LARGE_SEQ)])
        result = validate(fa)
        summary = result.summary()
        assert "VALID" in summary or "INVALID" in summary


# ── Cleaner tests ──────────────────────────────────────────────────────────

class TestCleaner:

    def test_removes_short_contigs(self, tmp_path):
        fa = write_fasta(tmp_path, [
            ("seq1", LARGE_SEQ),
            ("seq2", "ATGC" * 10),  # 40 bp
        ])
        result = clean(fa, min_contig_len=200)
        assert result.contigs_output == 1
        assert result.contigs_removed == 1

    def test_renames_headers(self, tmp_path):
        fa = write_fasta(tmp_path, [
            ("my_complex_header with spaces", LARGE_SEQ),
            ("another|header", LARGE_SEQ),
        ])
        result = clean(fa, rename_headers=True, prefix="contig")
        out_path = Path(result.output_path)
        content = out_path.read_text()
        assert ">contig_1" in content
        assert ">contig_2" in content

    def test_output_file_created(self, tmp_path):
        fa = write_fasta(tmp_path, [("seq1", LARGE_SEQ)])
        result = clean(fa)
        assert Path(result.output_path).exists()

    def test_custom_output_path(self, tmp_path):
        fa = write_fasta(tmp_path, [("seq1", LARGE_SEQ)])
        out = tmp_path / "custom_output.fasta"
        result = clean(fa, output_path=out)
        assert Path(result.output_path) == out
        assert out.exists()

    def test_keeps_all_valid_contigs(self, tmp_path):
        fa = write_fasta(tmp_path, [
            ("seq1", LARGE_SEQ),
            ("seq2", LARGE_SEQ),
            ("seq3", LARGE_SEQ),
        ])
        result = clean(fa, min_contig_len=200)
        assert result.contigs_output == 3
        assert result.contigs_removed == 0

    def test_summary_string(self, tmp_path):
        fa = write_fasta(tmp_path, [("seq1", LARGE_SEQ)])
        result = clean(fa)
        assert "Cleaning complete" in result.summary()


# ── Prokka availability test ───────────────────────────────────────────────

class TestProkkaAvailability:

    def test_returns_bool(self):
        result = prokka_available()
        assert isinstance(result, bool)
