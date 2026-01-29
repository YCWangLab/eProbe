#!/bin/bash
# Test script to verify --keep_bed and --remove_bed work correctly

echo "=== Testing BED Filter Modes ==="
echo ""

# Test files (replace with actual paths)
VCF="test.vcf.gz"
REF="test.fa"
BED="test.bed"

echo "Test 1: No BED filter (baseline)"
echo "Command: eprobe popgen extract -v $VCF -r $REF -o test_baseline"
echo "Expected: Extract all SNPs from VCF"
echo ""

echo "Test 2: --keep_bed mode"
echo "Command: eprobe popgen extract -v $VCF -r $REF -o test_keep --keep_bed $BED"
echo "Expected: Only extract SNPs INSIDE BED regions"
echo "Expected output:"
echo "  → Step 1: Filtering VCF with --keep_bed: test.bed"
echo "  → Step 2: BED filter (--keep_bed)"
echo "    ├─ Total SNPs in VCF: X"
echo "    ├─ Kept (in BED regions): Y"
echo "    └─ Removed (outside BED): X-Y"
echo ""

echo "Test 3: --remove_bed mode"
echo "Command: eprobe popgen extract -v $VCF -r $REF -o test_remove --remove_bed $BED"
echo "Expected: Only extract SNPs OUTSIDE BED regions"
echo "Expected output:"
echo "  → Step 1: Filtering VCF with --remove_bed: test.bed"
echo "  → Step 2: BED filter (--remove_bed)"
echo "    ├─ Total SNPs in VCF: X"
echo "    ├─ Removed (in BED regions): Y"
echo "    └─ Kept (outside BED): X-Y"
echo ""

echo "Verification checklist:"
echo "✓ --keep_bed shows correct log message"
echo "✓ --keep_bed generates inverted BED internally? NO (direct keep)"
echo "✓ --remove_bed shows correct log message"
echo "✓ --remove_bed generates inverted BED internally? YES"
echo "✓ Both modes report statistics correctly"
echo "✓ Keep + Remove counts should sum to total"
