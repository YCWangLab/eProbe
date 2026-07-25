# eProbe test data

A small benchmark dataset sampled from a rice population VCF (320 accessions, chr1–3).

## Files

| File | Description |
|------|-------------|
| `test.vcf.gz` / `test.vcf.gz.tbi` | Population VCF (bgzipped + tabix indexed) |
| `test.genome.fasta` / `*.fai` | Reference genome |
| `test.keep.bed` | Example keep-regions BED |
| `test.remove.bed` | Example remove-regions BED |
| `test.feature.bed` | CDS / gene feature BED |
| `test.gene.fasta` | Example gene sequences (FUNCGEN) |
| `test.allele.fasta` | Example allele sequences (FUNCGEN) |
| `test.exon.bed` | Exon BED (FUNCGEN) |

## POPGEN panel walkthrough

```bash
vcf=test_data/test.vcf.gz
ref=test_data/test.genome.fasta
keep_bed=test_data/test.keep.bed
remove_bed=test_data/test.remove.bed
feature_bed=test_data/test.feature.bed

# 1. Extract SNPs
eprobe popgen extract \
    -v $vcf -r $ref \
    --keep_bed $keep_bed \
    --remove_bed $remove_bed \
    -o test_out/step1

# 2. Filter (biophysical only; omit --BG_db / --AC_db if databases not available)
eprobe popgen filter \
    -i test_out/step1.snps.tsv \
    -r $ref \
    -o test_out/step2

# 3. Select – random (one SNP per 10 kb window)
eprobe popgen select \
    -i test_out/step2.filtered.tsv \
    --window_size 10000 \
    -o test_out/step3

# 3b. Select – biophysical-weighted
eprobe popgen select \
    -i test_out/step2.filtered.tsv \
    --window_size 10000 \
    --strategy weighted \
    --weights 0.2,0.2,0.2,0.2,0.2 \
    -o test_out/step3w

# 3c. Select – prioritise SNPs in CDS
eprobe popgen select \
    -i test_out/step2.filtered.tsv \
    --window_size 10000 \
    --priority_bed $feature_bed \
    -o test_out/step3p

# 4. Assess
eprobe popgen assess \
    -i test_out/step3.selected.tsv \
    --vcf $vcf \
    -o test_out/step4

# 5. Build probe FASTA
eprobe popgen build \
    -i test_out/step3.selected.tsv \
    -r $ref \
    --length 81 --snp_offset 0 \
    -o test_out/probes
```

## FUNCGEN panel walkthrough

```bash
gene=test_data/test.gene.fasta
allele=test_data/test.allele.fasta
bed=test_data/test.exon.bed
ref=test_data/test.genome.fasta
vcf=test_data/test.vcf.gz

# From FASTA (simple tiling)
eprobe funcgen from-fasta -f $gene -l 81 -s 30 -o test_out/funcgen_gene

# From FASTA (haplotype-aware; allele headers share prefix split by '_')
eprobe funcgen from-fasta -f $allele -l 81 -s 30 --haplotyping --sep _ -o test_out/funcgen_allele

# From BED (simple)
eprobe funcgen from-bed -b $bed -r $ref -l 81 -s 30 -o test_out/funcgen_bed

# From BED (haplotype-aware with VCF)
eprobe funcgen from-bed -b $bed -r $ref -v $vcf -l 81 -s 30 -o test_out/funcgen_hap
```
