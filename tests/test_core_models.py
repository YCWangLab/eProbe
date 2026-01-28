"""
Tests for eprobe.core.models module.
"""

import pytest
from eprobe.core.models import SNP, SNPDataFrame, Probe, ProbeSet


class TestSNP:
    """Tests for SNP dataclass."""
    
    def test_snp_creation_from_vcf_record(self):
        """Test creating a SNP from VCF-style record."""
        snp = SNP.from_vcf_record(
            chrom="chr1",
            pos=100,
            ref="A",
            alt="G",
        )
        
        assert snp.chrom == "chr1"
        assert snp.pos == 100
        assert snp.ref == "A"
        assert snp.alt == "G"
        assert snp.mutation_type == "ts"  # A->G is transition
    
    def test_snp_immutable(self):
        """Test that SNP is immutable (frozen)."""
        snp = SNP.from_vcf_record(
            chrom="chr1",
            pos=100,
            ref="A",
            alt="G",
        )
        
        with pytest.raises(AttributeError):
            snp.pos = 200
    
    def test_snp_id(self):
        """Test SNP ID generation."""
        snp = SNP.from_vcf_record(
            chrom="chr1",
            pos=100,
            ref="A",
            alt="G",
        )
        
        # ID format: chrom:pos_ref>alt
        assert snp.id == "chr1:100_A>G"
    
    def test_snp_invalid_ref(self):
        """Test SNP rejects invalid reference allele."""
        with pytest.raises(ValueError):
            SNP(chrom="chr1", pos=100, ref="X", alt="G", mutation_type="tv")
    
    def test_snp_invalid_position(self):
        """Test SNP rejects invalid position."""
        with pytest.raises(ValueError):
            SNP(chrom="chr1", pos=0, ref="A", alt="G", mutation_type="ts")
    
    def test_mutation_type_transition(self):
        """Test transition mutation type detection."""
        assert SNP.determine_mutation_type("A", "G") == "ts"
        assert SNP.determine_mutation_type("G", "A") == "ts"
        assert SNP.determine_mutation_type("C", "T") == "ts"
        assert SNP.determine_mutation_type("T", "C") == "ts"
    
    def test_mutation_type_transversion(self):
        """Test transversion mutation type detection."""
        assert SNP.determine_mutation_type("A", "C") == "tv"
        assert SNP.determine_mutation_type("A", "T") == "tv"
        assert SNP.determine_mutation_type("G", "C") == "tv"
        assert SNP.determine_mutation_type("G", "T") == "tv"


class TestSNPDataFrame:
    """Tests for SNPDataFrame class."""
    
    def test_from_snps(self):
        """Test creating DataFrame from SNP list."""
        snps = [
            SNP.from_vcf_record("chr1", 100, "A", "G"),
            SNP.from_vcf_record("chr1", 200, "C", "T"),
        ]
        
        df = SNPDataFrame.from_snps(snps)
        assert len(df) == 2
    
    def test_dataframe_iteration(self):
        """Test iterating over SNPs via DataFrame."""
        original_snps = [
            SNP.from_vcf_record("chr1", 100, "A", "G"),
            SNP.from_vcf_record("chr2", 200, "C", "T"),
        ]
        
        df = SNPDataFrame.from_snps(original_snps)
        
        # Check DataFrame length matches
        assert len(df) == len(original_snps)
        
        # Check data is stored correctly
        assert "chr1" in df._df["chr"].values
        assert "chr2" in df._df["chr"].values
    
    def test_to_tsv_and_from_tsv(self, temp_dir):
        """Test saving and loading TSV."""
        snps = [
            SNP.from_vcf_record("chr1", 100, "A", "G"),
        ]
        
        df = SNPDataFrame.from_snps(snps)
        tsv_path = temp_dir / "test.tsv"
        
        save_result = df.to_tsv(tsv_path)
        assert save_result.is_ok()
        
        load_result = SNPDataFrame.from_tsv(tsv_path)
        assert load_result.is_ok()
        
        loaded_df = load_result.unwrap()
        assert len(loaded_df) == 1
    
    def test_empty_dataframe(self):
        """Test creating empty SNPDataFrame."""
        df = SNPDataFrame()
        assert len(df) == 0
    
    def test_from_nonexistent_file(self, temp_dir):
        """Test loading from nonexistent file."""
        result = SNPDataFrame.from_tsv(temp_dir / "nonexistent.tsv")
        assert result.is_err()


class TestProbe:
    """Tests for Probe dataclass."""
    
    def test_probe_creation(self):
        """Test creating a Probe."""
        probe = Probe(
            id="PROBE001",
            sequence="ATCGATCGATCGATCG",
            source_chrom="chr1",
            source_start=100,
            source_end=115,
        )
        
        assert probe.id == "PROBE001"
        assert len(probe.sequence) == 16
        assert probe.length == 16
    
    def test_probe_gc_content(self):
        """Test probe GC content calculation."""
        probe = Probe(id="test", sequence="ATCGATCGATCGATCG")
        assert 45 <= probe.gc_content <= 55  # Approximately 50%
    
    def test_probe_fasta_record(self):
        """Test probe FASTA record formatting."""
        probe = Probe(id="test", sequence="ATCG")
        fasta = probe.to_fasta_record()
        assert fasta == ">test\nATCG"
    
    def test_probe_from_snp_region(self):
        """Test creating probe from SNP region."""
        snp = SNP.from_vcf_record("chr1", 100, "A", "G")
        probe = Probe.from_snp_region(
            snp=snp,
            sequence="ATCGATCGATCGATCGATCG",
            start=90,
            end=110,
        )
        
        assert probe.snp_pos == 100
        assert probe.snp_ref == "A"
        assert probe.snp_alt == "G"
        assert "chr1" in probe.id


class TestProbeSet:
    """Tests for ProbeSet class."""
    
    def test_probeset_creation(self):
        """Test creating ProbeSet."""
        probe_set = ProbeSet(name="test_set")
        
        probe1 = Probe(id="P1", sequence="ATCGATCG")
        probe2 = Probe(id="P2", sequence="GCTAGCTA")
        
        probe_set.add(probe1)
        probe_set.add(probe2)
        
        assert len(probe_set) == 2
    
    def test_probeset_get(self):
        """Test getting probe by ID."""
        probe_set = ProbeSet()
        probe = Probe(id="P1", sequence="ATCGATCG")
        probe_set.add(probe)
        
        retrieved = probe_set.get("P1")
        assert retrieved is not None
        assert retrieved.sequence == "ATCGATCG"
    
    def test_probeset_contains(self):
        """Test probe containment check."""
        probe_set = ProbeSet()
        probe = Probe(id="P1", sequence="ATCGATCG")
        probe_set.add(probe)
        
        assert "P1" in probe_set
        assert "P2" not in probe_set
    
    def test_probeset_iteration(self):
        """Test iterating over probes."""
        probe_set = ProbeSet()
        probe_set.add(Probe(id="P1", sequence="ATCGATCG"))
        probe_set.add(Probe(id="P2", sequence="GCTAGCTA"))
        
        ids = [p.id for p in probe_set]
        assert "P1" in ids
        assert "P2" in ids
