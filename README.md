# eProbes: a one-stop target genome capture baits design toolkit for ancient environmental DNA

This is the official repository for eProbes, which is a comprehensive python toolkit for designing target genome capture probe set for (ancient) environmental DNA. The core concept of this toolkit involves selecting an optimal probe set based on sequences containing informative SNPs and considering factors such as sequence accessibility across the species and the environment, analytical reliability, and expected performance in capture experiments.

# Installation
The eProbes packages and environments can be managed using conda (https://docs.conda.io/en/latest/) .

## Installation using conda
under developing

```
upcoming
```

## Build eProbes from development version by cloning the github repo
eProbes dependencies describe inserted here 

```
git clone --recursive https://github.com/YCWangLab/eProbes
cd eProbes
python setup.py install
```
## Dependency
```
Python >= v3.6
PyVCF>=0.6.8
Pysam>=0.11.2.2
Biopython >= v1.79
pandas>=1.1.5
numpy>=1.19.5
```

# Data preparation

The first step in probe generation is to prepare a Variant Call Format (VCF) file representing genetic variations. This requires users to sample the genetic diversity of the species as comprehensively as possible to meet the needs of population genetic analysis. Below is a simple and fast pipeline for generating the input VCF file. The main software tools used are BWA, samtools, picard, bcftools, and vcftools.

## Step1
Map reads from all individual against reference genome

```
# Indexing genome
bwa index ref_genome.fasta -p genome 
# Mapping 
bwa mem -t 24 -M genome -R "@RG\tID:your_indiv_id\tPL:ILLUMINA\tSM:your_indiv_id" R1.fastq R2.fastq > indiv.sam 
# Sorting
samtools view -bhF 4 indiv.sam > indiv.F4.bam && samtools sort indiv.F4.bam > indiv.sorted.bam 
# Deduplicating
java -jar picard.jar MarkDuplicates I=indiv.sorted.bam O=indiv.dd.bam M=indiv.dd.metrics REMOVE_DUPLICATES=true
```


## Step2
Perform SNP-calling (Tips: This can be done seperately on different chromosomes and different regions, and finally use ```bcftools concat``` to merge them into a VCF)

```
bcftools mpileup --threads 24 --fasta-ref ref_genome.fasta ./*.dd.bam | bcftools call -mv --threads 24 -o bcftools_raw.vcf &
```


## Step3
Filter, compress, and index VCF.

```
# keep bi-allelic variants only
vcftools --vcf  bcftools_raw.vcf  --max-missing 0.80 --minDP 3 --maf 0.05 --mac 3 --minQ 30 --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out MM80_filtered && bgzip -c MM80_filtered.recode.vcf > MM80.vcf.gz && tabix -p vcf MM80.vcf.gz &
```

## Step4 (optional)

We recommand users to conduct population genetic analyses with the VCF, such as distance-based tree, PCA, Admixture, computing site frequency spectrum, see whether it captures enough genetic diversity, and compare these with results conducted with a subset of SNPs covered by probes to evalute the performance of probes.


# Tutorial for running eProbes

Next, probe generation will be based on the VCF file, mainly involving preprocess, filtering, subsampling, and evaluation.

## Step1: eProbe_SNP_preprocessor
The eProbe_SNP_preprocessor module retrieves qualified SNPs from the VCF file (gz compressed and indexed). By default, eProbe exclude all SNPs near structural variants (InDels or SVs). Users can input a BED file to either retain or exclude SNPs in these regions, such as repetitive regions. Then, based on user-defined criteria, SNPs clustered together are filtered out, as these segments can increase  template bias of probes designed on them, reducing specificity when capturing non-template target sequences.

```
python eProbe_SNP_preprocessor.py -v vcf.gz -r ref_genome.fa -t threads 
```

## Step2: eProbe_SNP_filter
The eProbe_SNP_filter module performs default filtering processes, including Background filtering, Accessibility filtering, Taxonomic filtering, and Biophysical filtering.

Background filtering: Utilizes a kmer-based algorithm to remove probes that overlap with sequences in user-input databases. Users can specify dominant species in environmental samples, such as microbes or dominant plants, to reduce DNA from these sources overly occupying probes.

Accessibility filtering: Filters probes based on their accessibility within the intra-specific level.

Taxonomic filtering: Retains probes that can still be successfully assigned to the target species within allowable mismatch ranges, ensuring that sequences enriched in subsequent analyses can be reliably assigned to the target species.

Biophysical filtering: Filters probes based on GC content, melting temperature, complexity, and the likelihood of forming secondary structures.

```
python eProbe_SNP_filter.py -f SNPs_df.tsv(from step1) -r ref_genome.fa -t threads 
--BG_db Kraken2_databases --AC_db Bowtie2_databases --TX_db Bowtie2_databases --TX_taxa taxon_id  --TX_names_dmp names.dmp --TX_names_dmp nodes.dmp
```

## Step3: eProbe_SNP_subsampler
eProbe_SNP_subsampler allows for uniform sampling of qualified SNPs from the genome based on user-defined window sizes. By default, the program randomly samples within the window. Users can set weights for different biophysical tags to select the optimal SNPs within the window, or input a BED file to prioritize SNPs within the BED regions. Users can also specify a certain number of SNPs based on practical application needs.
```
python eProbe_SNP_subsampler.py -f SNPs_df.tsv(from step2) -r ref_genome.fa -t threads -o output
```

## Step4: eProbe_SNP_evaluator
eProbe_SNP_evaluator evaluates the final probe set in three main aspects:
Evaluate_by_tags: Generates distribution plots of different biophysical metrics after filtering.
Evaluate_by_distance: Computes and compares pairwise distances calculated based on original data and probe-covered SNPs. It visualizes distance matrices as heatmaps and hierarchical clustering, and computes matrix correlation and Manhattan distance.
Evaluate_by_missing: Simulates various degrees of missing data by randomly adding missing values to individuals. It assesses the probe set's tolerance to missing data on genetic typing.
```
python eProbe_SNP_evaluator.py -f SNPs_df.tsv(from step3) -r ref_genome.fa -v vcf.gz -t threads -o output
```


## Step5: eProbe_SNP_generator
eProbe_SNP_generator takes as input the SNPs_df.tsv provided by the user and generates probes accordingly. Users can specify the length of the probes, the positions of SNPs on the probes, whether to replace SNPs with bases other than REF and ALT, and whether to retain only transversions (tv) SNPs.
```
python eProbe_Seq_generator.py -f SNPs_df.tsv(from step3) -r ref_genome.fa --probe_length 81 --replace on --type tv -t threads -o output
```

## Step6: eProbe_Seq_generator
eProbe also allows users to input sequences of interest for probe generation. Users can input either a gene fasta file or a genome fasta file along with the positions where probes are needed (BED). Additionally, users can input a vcf file, and eProbe_Seq_generator can phase the regions within the BED file and obtain sequences of different alleles for probe design.
```
python eProbe_Seq_generator -f seq.fa --probe_length 81 --step 30
```
