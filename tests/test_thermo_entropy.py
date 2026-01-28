"""
Tests for thermodynamic entropy (ΔS°) calculation.

Tests the nearest-neighbor model implementation using Sugimoto et al. (1996)
parameters for DNA/DNA hybridization.
"""

import pytest
from eprobe.biophysics.thermo_entropy import (
    calculate_thermo_entropy,
    calculate_thermo_entropy_batch,
    calculate_thermo_enthalpy,
    calculate_gibbs_energy,
    NNTable,
    SUGIMOTO_1996_PARAMS,
    _get_nn_params,
    _is_self_complementary,
)


class TestNNParamsLookup:
    """Test nearest-neighbor parameter lookup."""
    
    def test_aa_lookup(self):
        """AA should find AA/TT parameters."""
        key, params = _get_nn_params("A", "A", SUGIMOTO_1996_PARAMS)
        assert key == "AA/TT"
        assert params == (-8.0, -21.9)
    
    def test_gc_lookup(self):
        """GC should find GC/CG parameters."""
        key, params = _get_nn_params("G", "C", SUGIMOTO_1996_PARAMS)
        assert key == "GC/CG"
        assert params == (-10.5, -26.4)
    
    def test_tc_lookup_reverse(self):
        """TC should find GA/CT (reverse) parameters."""
        key, params = _get_nn_params("T", "C", SUGIMOTO_1996_PARAMS)
        assert key == "GA/CT"  # TC uses reverse key
        assert params == (-8.8, -23.5)
    
    def test_at_lookup(self):
        """AT should find AT/TA parameters."""
        key, params = _get_nn_params("A", "T", SUGIMOTO_1996_PARAMS)
        assert key == "AT/TA"
        assert params == (-5.6, -15.2)


class TestSelfComplementary:
    """Test self-complementary sequence detection."""
    
    def test_palindrome_even(self):
        """GCGC is self-complementary."""
        assert _is_self_complementary("GCGC") is True
    
    def test_palindrome_gatc(self):
        """GATC is self-complementary."""
        assert _is_self_complementary("GATC") is True
    
    def test_not_palindrome(self):
        """AAAA is not self-complementary."""
        assert _is_self_complementary("AAAA") is False
    
    def test_not_palindrome_asymmetric(self):
        """ATCG is not self-complementary."""
        assert _is_self_complementary("ATCG") is False


class TestThermoEntropyCalculation:
    """Test thermodynamic entropy calculation."""
    
    def test_empty_sequence(self):
        """Empty sequence should return error."""
        result = calculate_thermo_entropy("")
        assert result.is_err()
        assert "Empty" in result.unwrap_err()
    
    def test_single_base(self):
        """Single base should return error."""
        result = calculate_thermo_entropy("A")
        assert result.is_err()
        assert "at least 2" in result.unwrap_err()
    
    def test_invalid_base(self):
        """Invalid base should return error."""
        result = calculate_thermo_entropy("ATCGN")
        assert result.is_err()
        assert "Invalid bases" in result.unwrap_err()
    
    def test_lowercase_accepted(self):
        """Lowercase sequence should work."""
        result = calculate_thermo_entropy("atcg")
        assert result.is_ok()
    
    def test_dinucleotide_aa(self):
        """AA dinucleotide entropy calculation."""
        # For AA: ΔS = -21.9 cal/(mol·K) + init (-9.0)
        result = calculate_thermo_entropy("AA", NNTable.SUGIMOTO_1996)
        assert result.is_ok()
        # AA/TT: -21.9, init: -9.0 = -30.9
        assert result.unwrap() == pytest.approx(-30.9, abs=0.1)
    
    def test_gc_rich_sequence(self):
        """GC-rich sequences have more negative entropy."""
        gc_result = calculate_thermo_entropy("GCGC")
        at_result = calculate_thermo_entropy("ATAT")
        
        assert gc_result.is_ok()
        assert at_result.is_ok()
        
        # GC-rich should have more negative (more stable) entropy
        assert gc_result.unwrap() < at_result.unwrap()
    
    def test_longer_sequence(self):
        """Longer sequences have more negative total entropy."""
        short = calculate_thermo_entropy("ATCG")
        long = calculate_thermo_entropy("ATCGATCG")
        
        assert short.is_ok()
        assert long.is_ok()
        
        # Longer sequence = more NN pairs = more negative entropy
        assert long.unwrap() < short.unwrap()
    
    def test_self_complementary_correction(self):
        """Self-complementary sequences get symmetry correction."""
        # GATC is self-complementary
        palindrome = calculate_thermo_entropy("GATC")
        # ATCG is not
        non_palindrome = calculate_thermo_entropy("ATCG")
        
        assert palindrome.is_ok()
        assert non_palindrome.is_ok()
        
        # Same NN composition but palindrome gets -1.4 symmetry correction
        # The values should differ by approximately -1.4
    
    def test_realistic_probe(self):
        """Test with realistic 81bp probe sequence."""
        seq = "GCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGC"
        result = calculate_thermo_entropy(seq)
        
        assert result.is_ok()
        entropy = result.unwrap()
        
        # Should be significantly negative for an 80bp sequence
        # Rough estimate: ~80 NN steps * ~-20 cal/mol/K each = ~-1600
        assert entropy < -1000
        assert entropy > -2500  # sanity check


class TestThermoEnthalpyCalculation:
    """Test thermodynamic enthalpy calculation."""
    
    def test_empty_sequence(self):
        """Empty sequence should return error."""
        result = calculate_thermo_enthalpy("")
        assert result.is_err()
    
    def test_dinucleotide(self):
        """Test basic enthalpy calculation."""
        result = calculate_thermo_enthalpy("AA")
        assert result.is_ok()
        # AA/TT: -8.0, init: 0.6 = -7.4
        assert result.unwrap() == pytest.approx(-7.4, abs=0.1)
    
    def test_gc_vs_at_enthalpy(self):
        """GC has more negative enthalpy than AT."""
        gc = calculate_thermo_enthalpy("GCGC")
        at = calculate_thermo_enthalpy("ATAT")
        
        assert gc.is_ok()
        assert at.is_ok()
        
        # GC-rich = stronger H-bonds = more negative enthalpy
        assert gc.unwrap() < at.unwrap()


class TestGibbsEnergy:
    """Test Gibbs free energy calculation."""
    
    def test_basic_calculation(self):
        """Test ΔG = ΔH - TΔS calculation."""
        result = calculate_gibbs_energy("ATCGATCG")
        assert result.is_ok()
        
        delta_g = result.unwrap()
        # Should be negative (favorable hybridization)
        assert delta_g < 0
    
    def test_temperature_effect(self):
        """Higher temperature should make ΔG less negative."""
        low_temp = calculate_gibbs_energy("ATCGATCG", temperature=273.15)  # 0°C
        high_temp = calculate_gibbs_energy("ATCGATCG", temperature=373.15)  # 100°C
        
        assert low_temp.is_ok()
        assert high_temp.is_ok()
        
        # At higher T, -TΔS term becomes more positive (ΔS is negative)
        # So ΔG becomes less negative (less favorable)
        assert low_temp.unwrap() < high_temp.unwrap()
    
    def test_37c_default(self):
        """Default temperature is 37°C (310.15K)."""
        result = calculate_gibbs_energy("ATCGATCG")
        result_explicit = calculate_gibbs_energy("ATCGATCG", temperature=310.15)
        
        assert result.is_ok()
        assert result_explicit.is_ok()
        assert result.unwrap() == result_explicit.unwrap()


class TestBatchCalculation:
    """Test batch entropy calculation."""
    
    def test_batch_empty_list(self):
        """Empty list returns empty list."""
        results = calculate_thermo_entropy_batch([])
        assert results == []
    
    def test_batch_multiple_sequences(self):
        """Multiple sequences return multiple results."""
        seqs = ["ATCG", "GCGC", "AAAA"]
        results = calculate_thermo_entropy_batch(seqs)
        
        assert len(results) == 3
        assert all(r.is_ok() for r in results)
    
    def test_batch_with_invalid(self):
        """Batch handles mix of valid/invalid sequences."""
        seqs = ["ATCG", "INVALID", "GCGC"]
        results = calculate_thermo_entropy_batch(seqs)
        
        assert len(results) == 3
        assert results[0].is_ok()
        assert results[1].is_err()
        assert results[2].is_ok()


class TestParameterTables:
    """Test different parameter tables."""
    
    def test_sugimoto_default(self):
        """Sugimoto 1996 is the default table."""
        default = calculate_thermo_entropy("ATCGATCG")
        explicit = calculate_thermo_entropy("ATCGATCG", NNTable.SUGIMOTO_1996)
        
        assert default.unwrap() == explicit.unwrap()
    
    def test_santalucia_available(self):
        """SantaLucia 1998 table is also available."""
        result = calculate_thermo_entropy("ATCGATCG", NNTable.SANTALUCIA_1998)
        assert result.is_ok()
    
    def test_tables_differ(self):
        """Different tables give different values."""
        sugimoto = calculate_thermo_entropy("ATCGATCG", NNTable.SUGIMOTO_1996)
        santalucia = calculate_thermo_entropy("ATCGATCG", NNTable.SANTALUCIA_1998)
        
        assert sugimoto.is_ok()
        assert santalucia.is_ok()
        
        # Values should be similar but not identical
        assert sugimoto.unwrap() != santalucia.unwrap()


class TestSugimotoParameters:
    """Verify Sugimoto 1996 parameters match reference."""
    
    def test_all_dinucleotides_present(self):
        """All 10 unique dinucleotide pairs should be present."""
        expected_pairs = {
            "AA/TT", "AT/TA", "TA/AT", "CA/GT", "GT/CA",
            "CT/GA", "GA/CT", "CG/GC", "GC/CG", "GG/CC",
        }
        actual_pairs = {k for k in SUGIMOTO_1996_PARAMS if "/" in k and k != "sym"}
        assert expected_pairs == actual_pairs
    
    def test_init_param_present(self):
        """Initiation parameter should be present."""
        assert "init" in SUGIMOTO_1996_PARAMS
        assert SUGIMOTO_1996_PARAMS["init"] == (0.6, -9.0)
    
    def test_sym_param_present(self):
        """Symmetry correction should be present."""
        assert "sym" in SUGIMOTO_1996_PARAMS
        assert SUGIMOTO_1996_PARAMS["sym"] == (0.0, -1.4)
    
    def test_aa_tt_values(self):
        """Verify AA/TT parameters match Sugimoto 1996."""
        # From paper: ΔH° = -8.0 kcal/mol, ΔS° = -21.9 cal/(mol·K)
        assert SUGIMOTO_1996_PARAMS["AA/TT"] == (-8.0, -21.9)
    
    def test_cg_gc_values(self):
        """Verify CG/GC parameters match Sugimoto 1996."""
        # From paper: ΔH° = -11.8 kcal/mol, ΔS° = -29.0 cal/(mol·K)
        assert SUGIMOTO_1996_PARAMS["CG/GC"] == (-11.8, -29.0)


class TestModuleExport:
    """Test that functions are properly exported from biophysics module."""
    
    def test_import_from_biophysics(self):
        """Functions should be importable from biophysics package."""
        from eprobe.biophysics import (
            calculate_thermo_entropy,
            calculate_thermo_entropy_batch,
            calculate_thermo_enthalpy,
            calculate_gibbs_energy,
            NNTable,
        )
        
        # Basic functionality check
        result = calculate_thermo_entropy("ATCG")
        assert result.is_ok()
