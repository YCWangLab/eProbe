# eProbes
This is the official developing repository for eProbes, a comprehensive python toolkit for designing target genome capture probe set for (ancient) environmental DNA.

# Installation
The eProbes packages and environments can be managed using conda (https://docs.conda.io/en/latest/) .

## Installation using conda
under developing

```
conda installation code inserted here
```

## Build eProbes from development version by cloning the github repo
eProbes dependencies describe inserted here 

```
git clone --recursive https://github.com/YCWangLab/eProbes
cd eProbes
make
```
## Dependency
```
Python v3.6
Biopython >= v1.79
PyVCF v0.6.0 (https://github.com/jamescasbon/PyVCF/)
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

Next, probe generation will be based on the VCF file, mainly involving candidate probe acquisition, calculation of various metrics (so-called tag in the following steps) for candidate probes, filtering, selection, and assessment.

## Test dataset
upcoming

## Step1
Retrieve location of all SNPs and filter out SNPs around InDels and SVs (if they are called in the same VCF), since probes targeted these SNPs might not be accessed from some population. Users can also input a BED format file to indicate SNPs in which genome regions need to be kept. The initial SNPs dataframe has  chromosome (chr), positon (pos), and type (transition or transversion) as the first three columns.

```
python Get_Accessible_SNPs.py -f ref_genome.fasta -v MM80.vcf.gz -d 100 -b genome_blocks.bed -o Accessible_SNPs_df.txt
```

## Step2
Draw distributions of GC content and melting temperature (Tm) with a batch of SNPs and probe length. The means and standard deviation of GC and Tm are useful in determing the length of probes and in selection step.

```
python Assess_Probes_By_GC_Tm.py -s Accessible_SNPs_df.txt -g ref_genome.fasta -k 71,121,10 -m Tm_NN > GC_Tm_report.txt
```
## Step3
Extract sequences of candidate probes.
```
python Get_Probes_Seq.py -s Accessible_SNPs_df.txt -g ref_genome.fasta -l 81 -o P81_pre_probes.fasta
```
## Step4
Attach tags for each probe.
```
# GC and Tm
python Get_Tag_GC_Tm.py -p P81_pre_probes.fasta -m "Tm_NN" -t 24 -o P81_Tag_GC_Tm.txt
# Hairpin structure
python Get_Tag_Hairpin.py -p P81_pre_probes.fasta -s 24 -t 24 -o P81_Tag_Hairpin.txt
# Dimer structure
python Get_Tag_Dimer.py -p P81_pre_probes.fasta -s 24 -t 24 -o P81_Tag_Dimer.txt
# Complexity
python Get_Tag_DUST.py -p P81_pre_probes.fasta -s 24 -t 24 -o P81_Tag_DUST.txt
```
Before obtaining the taxid, users need to run Kraken2 (https://github.com/DerrickWood/kraken2), a fast sequence classification tool based on k-mers, and use the output of Kraken2 as input. Considering that the majority of environmental DNA is often contaminated by microbial interference, it is suggested to maximize the diversity of microbial genomes when constructing the database.

Of course, this process might be time-consuming. Therefore, if users are unwilling to perform this step, they can generate a pseudo taxid file where the first and third columns are both 0, and the second column contains the IDs of all pre-probes.
``` 
# Running Kraken2
kraken2 --db your_database --output P81_output --classified-out P81_classified.fasta --threads 6 --report P81_report P81_pre_probes.fasta
# Taxid
python Get_Tag_Taxid.py -k P81_output -o P81_Tag_Taxid.txt
# Merged all tag files
python Merge_Tags.py -p P81_pre_probes.fasta -g P81_Tag_GC_Tm.txt -a P81_Tag_Hairpin.txt -d P81_Tag_Dimer.txt -u P81_Tag_DUST.txt -t P81_Tag_Taxid.txt -o P81_Merged_Tags.txt
```
## Step5 
Filter pre-probes based on given thresholds.
``` 
Filter_Probes_By_Tags.py -m Merged_Tags.txt -t 68,78 -g 40,60 -a 0,30 -u 0,2 -d 0.85 -i 0,4050,4051 -o Filtered_Tags.txt
```

## Step6
Select optimal probe from each window according to weights for each tags and keep SNPs in desired regions (identifying based on feature, such as exon/gene/mRNA). When enabling the ```both```mode and there are SNPs within a window that match feature, the program will prioritize selecting those SNPs.  
``` 
python Select_Probes_By_Tags.py -l Filtered_Tags.txt -t 75(desire_Tm) -g 50(desire_GC) -w 0.8,0.2,0,0,0(weights) -c Chr_length.txt -s 10000(win_size) -b bed.txt -n 8 -f exon -m both -r 24 -o P81
```

## Step7
Generate VCF for evaluation performance of probes
```
# SNPs covered by probes
python Probes_to_VCF.py -v MM80.vcf.gz -o Evaluation.vcf -f P81_merged_SNPs.txt
# Set random missing to VCF
python Assess_By_Missing.py -v Evaluation.vcf -r 0.2,0.4,0.6 -g ./.:0,0,0 -o M246.simulate.vcf
```

## Step8
Generate probes
```
python Get_Probes_Seq.py -s P81_merged_SNPs.txt(your_final_SNPs_df.txt) -g ref_genome.fasta -l 81 -o P81_probes.fasta
```
