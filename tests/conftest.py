"""
Test configuration and fixtures.
"""

import pytest
from pathlib import Path
import tempfile
import shutil


@pytest.fixture
def temp_dir():
    """Create a temporary directory for test outputs."""
    tmpdir = tempfile.mkdtemp()
    yield Path(tmpdir)
    shutil.rmtree(tmpdir)


@pytest.fixture
def sample_fasta(temp_dir):
    """Create a sample FASTA file."""
    fasta_path = temp_dir / "test.fa"
    fasta_path.write_text(
        ">chr1\n"
        "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n"
        "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n"
        ">chr2\n"
        "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA\n"
        "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA\n"
    )
    return fasta_path


@pytest.fixture
def sample_snp_tsv(temp_dir):
    """Create a sample SNP TSV file."""
    tsv_path = temp_dir / "snps.tsv"
    tsv_path.write_text(
        "snp_id\tchrom\tpos\tref\talt\tleft_flank\tright_flank\n"
        "SNP001\tchr1\t100\tA\tG\tACGTACGTACGT\tGCTAGCTAGCTA\n"
        "SNP002\tchr1\t200\tC\tT\tTGCATGCATGCA\tACGTACGTACGT\n"
        "SNP003\tchr2\t150\tG\tA\tCGTACGTACGTA\tTAGCTAGCTAGC\n"
    )
    return tsv_path


@pytest.fixture
def sample_sequences():
    """Sample sequences for testing."""
    return {
        "seq1": "ATCGATCGATCGATCGATCG",
        "seq2": "GCTAGCTAGCTAGCTAGCTA",
        "seq3": "AAAAAAAAAAAAAAAAAAAA",
        "seq4": "GCGCGCGCGCGCGCGCGCGC",
    }
