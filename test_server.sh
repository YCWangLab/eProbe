#!/bin/bash
# Server test script for eprobe
# Captures all output for debugging

set -e

# Configuration - modify these paths for your server
KRAKEN2_DB="/path/to/kraken2_db"
NGSLCA_DB="/path/to/ngslca_db"
BOWTIE2_DB="/path/to/bowtie2_db"
TEST_VCF="/path/to/test.vcf.gz"
TEST_REF="/path/to/reference.fa"
OUTPUT_DIR="./test_output_$(date +%Y%m%d_%H%M%S)"

# Create output directory
mkdir -p "$OUTPUT_DIR"
LOG_FILE="$OUTPUT_DIR/test_run.log"

echo "=== eprobe Server Test ===" | tee "$LOG_FILE"
echo "Date: $(date)" | tee -a "$LOG_FILE"
echo "Python: $(which python)" | tee -a "$LOG_FILE"
echo "Output dir: $OUTPUT_DIR" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

# Function to run and log commands
run_test() {
    local test_name="$1"
    local cmd="$2"
    
    echo "----------------------------------------" | tee -a "$LOG_FILE"
    echo "TEST: $test_name" | tee -a "$LOG_FILE"
    echo "CMD: $cmd" | tee -a "$LOG_FILE"
    echo "----------------------------------------" | tee -a "$LOG_FILE"
    
    # Run command and capture both stdout and stderr
    if eval "$cmd" >> "$LOG_FILE" 2>&1; then
        echo "✓ PASSED: $test_name" | tee -a "$LOG_FILE"
    else
        echo "✗ FAILED: $test_name (exit code: $?)" | tee -a "$LOG_FILE"
    fi
    echo "" | tee -a "$LOG_FILE"
}

# Test 1: Check installation
run_test "Version check" "python -m eprobe --version"

# Test 2: Check info
run_test "Info command" "python -m eprobe info"

# Test 3: FUNCGEN from-fasta (no external tools needed)
run_test "funcgen from_fasta" \
    "python -m eprobe funcgen from_fasta -f test_data/test.gene.fasta -o $OUTPUT_DIR/funcgen_test -l 81 -s 30"

# Test 4: POPGEN extract (if VCF available)
if [[ -f "$TEST_VCF" && -f "$TEST_REF" ]]; then
    run_test "popgen extract" \
        "python -m eprobe popgen extract -v $TEST_VCF -r $TEST_REF -o $OUTPUT_DIR/popgen_extract -t 4"
fi

# Test 5: POPGEN filter with Kraken2 (if available)
if [[ -d "$KRAKEN2_DB" && -f "$OUTPUT_DIR/popgen_extract.snps.tsv" ]]; then
    run_test "popgen filter (kraken2)" \
        "python -m eprobe popgen filter -i $OUTPUT_DIR/popgen_extract.snps.tsv -r $TEST_REF -o $OUTPUT_DIR/popgen_filter --bg_db $KRAKEN2_DB -t 4"
fi

# Test 6: Biophysics calculations
run_test "biophysics test" \
    "python -c \"
from eprobe.biophysics import calculate_thermo_entropy, calculate_tm, calculate_gc
seq = 'GCGATCGATCGATCGATCGATCGATCGATCGATCGATCG'
print(f'Sequence: {seq}')
print(f'GC: {calculate_gc(seq).unwrap():.2f}%')
print(f'Tm: {calculate_tm(seq).unwrap():.2f}°C')
print(f'ΔS°: {calculate_thermo_entropy(seq).unwrap():.2f} cal/mol/K')
\""

echo "========================================" | tee -a "$LOG_FILE"
echo "Tests completed. Log saved to: $LOG_FILE" | tee -a "$LOG_FILE"
echo "========================================" | tee -a "$LOG_FILE"

# Create a summary for easy sharing
echo ""
echo "To share results, copy this file:"
echo "  $LOG_FILE"
echo ""
echo "Or run: cat $LOG_FILE | pbcopy  (on macOS)"
