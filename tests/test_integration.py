"""
Integration tests using real test data files.

Tests the complete pipeline functionality without external dependencies
like Kraken2 or ngsLCA.
"""

import pytest
from pathlib import Path
from eprobe.core.fasta import read_fasta, write_fasta
from eprobe.core.models import SNP, SNPDataFrame
from eprobe.funcgen.from_fasta import tile_sequence, TilingConfig, parse_haplotype_id, group_haplotypes, HaplotypeConfig
from eprobe.funcgen.from_bed import parse_bed_file
from eprobe.biophysics import calculate_gc, calculate_tm, calculate_complexity, calculate_hairpin_score
from eprobe.biophysics.entropy import calculate_entropy


# Test data directory
TEST_DATA_DIR = Path(__file__).parent.parent / "test_data"


class TestFastaOperations:
    """Test FASTA file operations with real test data."""
    
    def test_read_genome_fasta(self):
        """Test reading genome FASTA."""
        fasta_path = TEST_DATA_DIR / "test.genome.fasta"
        if not fasta_path.exists():
            pytest.skip("Test data not available")
        
        result = read_fasta(fasta_path)
        assert result.is_ok()
        
        sequences = result.unwrap()
        assert "chr1" in sequences
        assert "chr2" in sequences
        assert len(sequences["chr1"]) > 0
    
    def test_read_gene_fasta(self):
        """Test reading gene FASTA."""
        fasta_path = TEST_DATA_DIR / "test.gene.fasta"
        if not fasta_path.exists():
            pytest.skip("Test data not available")
        
        result = read_fasta(fasta_path)
        assert result.is_ok()
        
        sequences = result.unwrap()
        assert len(sequences) == 3  # GENE1, GENE2, GENE3


class TestBedOperations:
    """Test BED file operations with real test data."""
    
    def test_parse_exon_bed(self):
        """Test parsing exon BED file."""
        bed_path = TEST_DATA_DIR / "test.exon.bed"
        if not bed_path.exists():
            pytest.skip("Test data not available")
        
        result = parse_bed_file(bed_path)
        assert result.is_ok()
        
        regions = result.unwrap()
        assert len(regions) == 5
    
    def test_parse_keep_bed(self):
        """Test parsing keep regions BED file."""
        bed_path = TEST_DATA_DIR / "test.keep.bed"
        if not bed_path.exists():
            pytest.skip("Test data not available")
        
        result = parse_bed_file(bed_path)
        assert result.is_ok()
        
        regions = result.unwrap()
        assert len(regions) == 3


class TestSNPDataFrameIntegration:
    """Test SNPDataFrame with real test data."""
    
    def test_load_snp_tsv(self):
        """Test loading SNP TSV file."""
        tsv_path = TEST_DATA_DIR / "test.snps.tsv"
        if not tsv_path.exists():
            pytest.skip("Test data not available")
        
        result = SNPDataFrame.from_tsv(tsv_path)
        assert result.is_ok()
        
        df = result.unwrap()
        assert len(df) == 19


class TestFuncgenFromFasta:
    """Test FUNCGEN probe generation from FASTA."""
    
    def test_tile_gene_sequences(self):
        """Test tiling gene sequences."""
        fasta_path = TEST_DATA_DIR / "test.gene.fasta"
        if not fasta_path.exists():
            pytest.skip("Test data not available")
        
        result = read_fasta(fasta_path)
        assert result.is_ok()
        sequences = result.unwrap()
        
        config = TilingConfig(probe_length=81, step_size=30)
        
        total_probes = []
        for seq_id, sequence in sequences.items():
            probes = tile_sequence(seq_id, sequence, config)
            total_probes.extend(probes)
        
        assert len(total_probes) > 0
        # Check all probes are 81bp
        for probe in total_probes:
            assert probe.length == 81
    
    def test_haplotype_grouping(self):
        """Test haplotype-aware probe generation."""
        fasta_path = TEST_DATA_DIR / "test.allele.fasta"
        if not fasta_path.exists():
            pytest.skip("Test data not available")
        
        result = read_fasta(fasta_path)
        assert result.is_ok()
        sequences = result.unwrap()
        
        # Check haplotype ID parsing
        gene_id, allele_id = parse_haplotype_id("BRCA1_1", "_")
        assert gene_id == "BRCA1"
        assert allele_id == "1"
        
        # Group by gene
        config = HaplotypeConfig(enabled=True, separator="_")
        groups = group_haplotypes(sequences, config)
        
        assert "BRCA1" in groups
        assert "BRCA2" in groups
        assert len(groups["BRCA1"]) == 2  # Two alleles


class TestBiophysicsWithRealSequences:
    """Test biophysics calculations on real sequences."""
    
    def test_calculate_biophysics_on_probes(self):
        """Calculate biophysics properties on generated probes."""
        fasta_path = TEST_DATA_DIR / "test.gene.fasta"
        if not fasta_path.exists():
            pytest.skip("Test data not available")
        
        result = read_fasta(fasta_path)
        assert result.is_ok()
        sequences = result.unwrap()
        
        config = TilingConfig(probe_length=81, step_size=30)
        
        for seq_id, sequence in sequences.items():
            probes = tile_sequence(seq_id, sequence, config)
            
            for probe in probes:
                # Calculate all biophysics properties
                gc_result = calculate_gc(probe.sequence)
                assert gc_result.is_ok()
                gc = gc_result.unwrap()
                assert 0 <= gc <= 100
                
                tm_result = calculate_tm(probe.sequence)
                assert tm_result.is_ok()
                tm = tm_result.unwrap()
                assert tm > 0
                
                complexity_result = calculate_complexity(probe.sequence)
                assert complexity_result.is_ok()
                
                hairpin_result = calculate_hairpin_score(probe.sequence)
                assert hairpin_result.is_ok()
                
                entropy_result = calculate_entropy(probe.sequence)
                assert entropy_result.is_ok()


class TestEndToEndWorkflow:
    """Test complete end-to-end workflow without external dependencies."""
    
    def test_funcgen_workflow(self, temp_dir):
        """Test complete FUNCGEN workflow."""
        fasta_path = TEST_DATA_DIR / "test.gene.fasta"
        if not fasta_path.exists():
            pytest.skip("Test data not available")
        
        # Step 1: Read input
        result = read_fasta(fasta_path)
        assert result.is_ok()
        sequences = result.unwrap()
        
        # Step 2: Generate probes
        config = TilingConfig(probe_length=81, step_size=30)
        all_probes = []
        for seq_id, sequence in sequences.items():
            probes = tile_sequence(seq_id, sequence, config)
            all_probes.extend(probes)
        
        assert len(all_probes) > 0
        
        # Step 3: Calculate biophysics
        probe_data = []
        for probe in all_probes:
            gc_result = calculate_gc(probe.sequence)
            tm_result = calculate_tm(probe.sequence)
            
            probe_data.append({
                "id": probe.id,
                "sequence": probe.sequence,
                "gc": gc_result.unwrap() if gc_result.is_ok() else None,
                "tm": tm_result.unwrap() if tm_result.is_ok() else None,
            })
        
        # Step 4: Write output
        output_sequences = {p["id"]: p["sequence"] for p in probe_data}
        output_path = temp_dir / "probes.fasta"
        
        write_result = write_fasta(output_sequences, output_path)
        assert write_result.is_ok()
        assert output_path.exists()
        
        # Verify output
        read_back_result = read_fasta(output_path)
        assert read_back_result.is_ok()
        assert len(read_back_result.unwrap()) == len(all_probes)
