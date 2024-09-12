# Install 
```
git clone https://github.com/YCWangLab/eProbe.git
cd eProbe
python setup.py install
chmod 755 *py 
export PATH=/path_to_eProbe/:$PATH
```

# Dependencies
```
conda create -n eProbe python>=3.8 bowtie2 samtools kraken2 bedtools plink bcftools shapeit clustal htslib -c bioconda
```
If there are conflicts (which is very common), install separately with `conda` or manually and configure them to `$PATH`

# POPGEN panel
## test data
```
wget -O test_data.zip https://sid.erda.dk/share_redirect/GuARqfEdiI
```
```
vcf=/test_data/test.vcf.gz
ref=/test_data/test.genome.fasta
keep_bed=/test_data/test.keep.bed
remove_bed=/test_data/test.remove.bed
feature_bed=/test_data/test.feature.bed
```
## step 1: preprocessing

```
python SNP_preprocessor.py -v $vcf -r $ref --max_clump_snp 3 --cluster_flank 60 --keep_bed $keep_bed --keep_distance -10 --remove_bed $remove_bed --rm_distance 10 -t 6 -o test_step1
```

## step 2: filtering
### 2.1 background noise filtering
```
# kraken2 database for background noise filtering
BG_db1=/database/kraken2/bacteria
BG_db2=/database/kraken2/prok
BG_db3=/database/kraken2/dominent_plant
```
configure kraken2 database: https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown

```
python SNP_filter.py -f test_step1.preprocessed_SNPs.tsv -r $ref -t 6 --BG_db $BG_db1,$BG_db2,$BG_db3 --BG_min_hits 3 -o test_step2.1 
```
### 2.2 accessibility filtering
```
# bowtie2 index
AC_db1=Nipponbare
AC_db2=Indica
AC_db3=Rufipogen
python SNP_filter.py -f test_step1.preprocessed_SNPs.tsv/test_step2.1.filtered_probes.tsv -r $ref -t 6 --AC_db $AC_db1,$AC_db2,$AC_db3 -o test_step2.2
```

### 2.3 test taxonomic filtering
download **Poaceae** assemblies from NCBI Refseq using **ncbi-genome-download**  (doi: 10.5281/zenodo.8192432) and index using **Bowtie2**.
taxonomy files downloaded from: https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump
```
# bowtie2 index
db=Poaceae
names_dmp=names.dmp
nodes_dmp=nodes.dmp
acc2tax=accession2taxid.txt
python SNP_filter.py -f test_step1.preprocessed_SNPs.tsv/test_step2.2.filtered_probes.tsv -r $ref -t 6 --TX_db $db --TX_taxa 4530(rice taxid) --TX_names_dmp $names_dmp --TX_nodes_dmp $nodes_dmp --TX_acc2tax $acc2tax --TX_minedit 0 --TX_maxedit 2 -o test_step2.3
```

### 2.4 biophysical/biochemical filtering
```
python SNP_filter.py -f test_step1.preprocessed_SNPs.tsv/test_step2.3.filtered_probes.tsv -r $ref -t 6 --BFfilter on -o test_step2.4
```

### Run all filtering in one go
```
python SNP_filter.py -f test_step1.preprocessed_SNPs.tsv -r $ref -t 6 --BG_db $BG_db1,$BG_db2,$BG_db3 --BG_min_hits 2 --AC_db $AC_db1,$AC_db2,$AC_db3 --TX_db $db --TX_taxa 4530 --TX_names_dmp $names_dmp --TX_nodes_dmp $nodes_dmp --TX_acc2tax $acc2tax --TX_minedit 0 --TX_maxedit 2 --BFfilter on -o test_step2
```
## step 3: subsampling
### 3.1 random selection
```
python SNP_subsampler.py -f test_step2.filtered_probes.tsv -r $ref -t 6 --window_size 10000 -o test_step3
```
### 3.2 optimization
```
python SNP_subsampler.py -f test_step2.filtered_probes.tsv -r $ref -t 6 --window_size 10000 --select_weights 0.2,0.2,0.2,0.2,0.2 -o test_step3 
```
### 3.3 prioritizing SNPs in CDS
```
python SNP_subsampler.py -f test_step2.filtered_probes.tsv -r $ref -t 6 --window_size 10000 --keep_bed $feature_bed --keep_bed_column 4 --keep_bed_feature CDS --output_mode both -o test_step3
```
### 3.4 downsize
```
python SNP_subsampler.py -f test_step2.filtered_probes.tsv -r $ref -t 6 --window_size 10000 --probe_number 1000 -o test_step3 
```
### Run all subsampling function can be run in one go

## step 4: assessing
### 4.1 biophysical/biochemical features
```
python SNP_accessor.py -d snp_df.tsv --assessor tag --reference $ref --tag all -o test_step4.1
```
### 4.2 pair-wise distance
```
python SNP_accessor.py -d snp_df.tsv --assessor distance --vcf $vcf -o test_step4.2
```
### 4.3 missing rate
```
python SNP_accessor.py -d snp_df.tsv --assessor missing --vcf $vcf --pop_data pop_data.txt -o test_step4.3
```

## step 5: generating
### 6.1 pos (default)
* use when you want to adjust the lenght or the position of SNPs on the probe
```
python SNP_generator.py -d snp_df.tsv -r $ref --method pos -l 81 -s 0 --replace on -o test_step5.1
```
### 6.2 edge
```
python SNP_generator.py -d snp_df.tsv -r $ref --method edge --replace on -o test_step5.2
```

# FUNCGEN panel
## test data
```
gene=/test_data/test.gene.fasta
allele=/test_data/test.allele.fasta
bed=/test_data/test.exon.bed
ref=/test_data/test.genome.fasta
```
### 7.1 fasta file
#### gene sequences without alleles
```
python Seq_generator.py -m seq -l 81 -s 30 -f $gene -o test_step7.1
```
#### allele sequences with alleles
```
python Seq_generator.py -m seq --haplotyping on -l 81 -s 30 -f $allele --sep _ -o test_step7.1
```

### 2 bed file 
#### without inferring haplotype
```
python Seq_generator.py -m bed -b $bed -r $ref -l 81 -s 30 -o test_step7.2
```
#### inferring haplotype

```
python Seq_generator.py -m bed --haplotyping on -b $bed -r $ref -v $vcf -l 81 -s 30 -o test_step7.2
```