# Call Chain Verification for --keep_bed and --remove_bed

## Complete Call Chain

### 1. CLI Layer (popgen.py)
```
User runs: eprobe popgen extract --keep_bed regions.bed
         OR: eprobe popgen extract --remove_bed regions.bed

↓

Lines 166-174: Set bed_filter_type and bed_path_used
- keep_bed  → bed_filter_type = "keep", bed_path_used = keep_bed
- remove_bed → bed_filter_type = "remove", bed_path_used = remove_bed

↓

Lines 180-190: Call run_extract()
run_extract(
    bed_path=bed_path_used,           # ✓ Path passed
    bed_mode=bed_filter_type,         # ✓ Mode passed ("keep" or "remove")
)
```

### 2. Extract Layer (extract.py)
```
Lines 260-271: run_extract() receives parameters
def run_extract(
    bed_path: Optional[Path] = None,   # ✓ Received
    bed_mode: str = "keep",            # ✓ Received
)

↓

Lines 323-328: Call extract_snps_from_vcf()
vcf_result = extract_snps_from_vcf(
    vcf_path, 
    bed_path,      # ✓ Path passed
    bed_mode,      # ✓ Mode passed
    fai_path=fai_path, 
    threads=threads
)

↓

Lines 334-337: Extract results
raw_snps = vcf_data['snps']
bed_applied = vcf_data['bed_applied']    # ✓ Returned
bed_mode_applied = vcf_data['bed_mode']  # ✓ Returned
```

### 3. VCF Layer (vcf.py)
```
Lines 263-268: extract_snps_from_vcf() receives parameters
def extract_snps_from_vcf(
    bed_path: Optional[Path] = None,   # ✓ Received
    bed_mode: str = "keep",            # ✓ Received
    fai_path: Optional[Path] = None,   # ✓ Received
)

↓

Lines 303-307: Track BED application
bed_applied = bed_path is not None     # ✓ Set
bed_mode_used = bed_mode if bed_applied else None  # ✓ Set

↓

Lines 308-325: Process BED based on mode
if bed_path:
    if bed_mode == "remove":
        # ✓ Generate inverted BED
        bed_regions = invert_bed_regions(bed_path, fai_path)
        logger.info("Generated {len} inverted regions")
    else:
        # ✓ Read BED directly (keep mode)
        bed_regions = [(chrom, start, end), ...]
        logger.info("Loaded {len} BED regions (keep mode)")

↓

Lines 342-344 OR 348-350: Pass to chromosome processing
_extract_snps_for_chromosome(
    vcf_path, chrom, bed_regions  # ✓ Regions passed (inverted or direct)
)

↓

Lines 365-368: Return results with metadata
result = {
    'snps': all_snps,               # ✓ Filtered SNPs
    'bed_applied': bed_applied,     # ✓ True/False
    'bed_mode': bed_mode_used,      # ✓ "keep"/"remove"/None
}
return Ok(result)
```

### 4. Chromosome Processing Layer (vcf.py)
```
Lines 88-95: _extract_snps_for_chromosome() receives regions
def _extract_snps_for_chromosome(
    bed_regions: Optional[List[Tuple[str, int, int]]],  # ✓ Received
)

↓

Lines 117-132: Use regions for indexed VCF access
if bed_regions:
    # ✓ Iterate through regions (works for both keep and inverted-remove)
    for region_chrom, region_start, region_end in bed_regions:
        for variant in vcf(f"{region_chrom}:{region_start}-{region_end}"):
            # Extract SNPs from this region
```

## Verification Summary

✅ **--keep_bed path:**
1. CLI sets bed_mode="keep"
2. vcf.py reads BED directly
3. Uses regions as-is for VCF indexed query
4. Returns SNPs inside BED regions

✅ **--remove_bed path:**
1. CLI sets bed_mode="remove"
2. vcf.py calls invert_bed_regions(bed_path, fai_path)
3. Uses inverted regions for VCF indexed query
4. Returns SNPs outside original BED regions

✅ **Statistics path:**
1. extract.py calls count_snps_in_vcf() for total
2. Calculates kept = len(raw_snps)
3. Calculates removed = total - kept
4. CLI displays based on bed_mode

## Key Insight

Both modes use the SAME code path in _extract_snps_for_chromosome().
The difference is:
- keep: bed_regions = original BED
- remove: bed_regions = inverted BED

This ensures identical performance for both modes!
