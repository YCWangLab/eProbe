import pybedtools
import pysam
import subprocess
import multiprocessing
import os
import shutil
import datetime
from pipe_controler import make_temp_dir
from fasta_operator import write_fasta


# phase an interval (i.e., a BED object) of a whole VCF
def phase_common(bed_obj, vcf_path, tmp_dir, name):
    chrom = bed_obj.chrom
    start = bed_obj.start
    end = bed_obj.end

    output_bcf = os.path.join(tmp_dir, f'{name}.bcf')
    output_vcf = os.path.join(tmp_dir, f'{name}.vcf')
    output_vcf_gz = f'{output_vcf}.gz'

    success = True

    if not os.path.exists(output_bcf):
        try:
            shapeit_command = f'phase_common --input {vcf_path} --output-format bcf --output {output_bcf} --region {chrom}:{start}-{end} --thread 1'
            subprocess.run(shapeit_command, shell=True, check=True)

            bcftools_command = f'bcftools view {output_bcf} -Ov -o {output_vcf}'
            subprocess.run(bcftools_command, shell=True, check=True)

            bgzip_command = f'bgzip -c {output_vcf} > {output_vcf_gz}'
            subprocess.run(bgzip_command, shell=True, check=True)

            tabix_command = f'tabix -f {output_vcf_gz}'
            subprocess.run(tabix_command, shell=True, check=True)
        except subprocess.CalledProcessError:
            print(f"Error occurred during phasing in region {chrom}:{start}-{end}. "
                  f"Possibly because no variants in this regions. Skipping...")
            success = False
    else:
        print(f'Duplicated genomic region detected: {name}. Skipping phasing...')

    return 1 if success else 0


def replace_dna_sequence(template, pos_r, pos_snp, ref, allele):
    # pos_r: relative position within bed; pos_snp: genome position
    ref_len = len(ref)
    allele_len = len(allele)

    # Handle substitution
    if ref_len == allele_len == 1:
        template[pos_r] = allele

    # Handle deletion (ref is longer than allele)
    elif ref_len > allele_len:
        start = pos_r
        end = pos_r + ref_len
        template[start:end] = [''] * (end - start) 
        template[pos_r] = allele 

    # Handle insertion (allele is longer than ref)
    elif allele_len > ref_len:
        start = pos_r
        end = pos_r + ref_len
        template[start:end] = [''] * (end - start)
        template.insert(pos_r, allele)

    else:
        print(f"Unhandled case at position {pos_snp}: ref='{ref}' allele='{allele}'")

    return template


def calculate_frequencies(input_list):
    frequency_dict = {}
    total_count = len(input_list)

    for item in input_list:
        frequency_dict[item] = frequency_dict.get(item, 0) + 1

    for item in frequency_dict:
        frequency_dict[item] /= total_count

    return frequency_dict


def get_bed_sequence(chrom, start, end, fasta_obj):
    sequence = fasta_obj.fetch(chrom, start, end)
    return sequence


def get_allele_sequences(bed_obj, fasta_obj, bed_vcf):
    vcf_reader = pysam.VariantFile(bed_vcf)

    bed_chrom = bed_obj.chrom
    bed_start = bed_obj.start
    bed_end = bed_obj.end
    ref_sequence = get_bed_sequence(bed_chrom, bed_start - 1, bed_end, fasta_obj) # use 0-based position

    sample_list = list(vcf_reader.header.samples)
    sample_variants = {sample: {} for sample in sample_list} 
    sample_sequences = {sample: {0: list(ref_sequence), 1: list(ref_sequence)} for sample in sample_list}

    for record in vcf_reader.fetch(bed_chrom, bed_start, bed_end + 1):
        pos_snp = record.pos
        # relative position in the interval
        pos_r = pos_snp - bed_start # 0-based position

        ref = record.ref
        alts = record.alts
        alleles = [ref] + list(alts)

        for sample in record.samples:
            GT1 = record.samples[sample]['GT'][0]
            GT2 = record.samples[sample]['GT'][1]
            if GT1 is not None and GT2 is not None:
                sample_variants[sample][f'{pos_r}_{pos_snp}'] = {'ref': ref, 0: alleles[GT1], 1: alleles[GT2]}

    for sample, variants in sample_variants.items():
        for pos, variant in variants.items():
            # for haplotype 0
            if variant['ref'] != variant[0]:
                sample_sequences[sample][0] = replace_dna_sequence(sample_sequences[sample][0], int(pos.split("_")[0]), int(pos.split("_")[1]),variant['ref'], variant[0])
            # for haplotype 1
            if variant['ref'] != variant[1]:
                sample_sequences[sample][1] = replace_dna_sequence(sample_sequences[sample][1], int(pos.split("_")[0]), int(pos.split("_")[1]),variant['ref'], variant[1])

    sample_sequences = {
        sample: {0: ''.join(haplotypes[0]), 1: ''.join(haplotypes[1])}
        for sample, haplotypes in sample_sequences.items()
    }

    bed_alleles = calculate_frequencies(
        [sample_allele for sample in sample_sequences.values() for sample_allele in sample.values()]
    )

    return bed_alleles


def lanch_process(bed_index, bed_file, whole_vcf, genome_file, tmp_dir):
    fasta_obj = pysam.FastaFile(genome_file)
    BED = pybedtools.BedTool(bed_file)
    bed_obj = BED[bed_index]
    # name of the interval
    if len(bed_obj.name) <= 1:
        bed_obj_name = f'{bed_obj.chrom}:{bed_obj.start}-{bed_obj.end}'
    else:
        bed_obj_name = bed_obj.name

    # has variants in this region?
    v = phase_common(bed_obj, whole_vcf, tmp_dir, bed_obj_name)
    if v == 1:
        bed_vcf = os.path.join(tmp_dir, f'{bed_obj_name}.vcf.gz')
        bed_alleles = get_allele_sequences(bed_obj, fasta_obj, bed_vcf)
    else:
        bed_chrom = bed_obj.chrom
        bed_start = bed_obj.start - 1
        bed_end = bed_obj.end - 1
        ref_sequence = get_bed_sequence(bed_chrom, bed_start, bed_end, fasta_obj)
        bed_alleles = calculate_frequencies([ref_sequence])
    return bed_obj_name, bed_alleles


if __name__ == '__main__':
    import argparse

    usage = '''python Get_allele_from_vcf.py -v vcf.gz -r genome.fasta -b genome_block.bed -t threads -o output.fasta'''

    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument('-v', '--vcf', dest='vcf', action='store', nargs='?',
                        help='Input vcf.gz (with .tbi) for phasing with shapeit (should be in your $PATH).',
                        metavar='FILE', required=True)
    parser.add_argument("-r", "--ref", dest="ref", action="store", nargs='?',
                        help="Input reference genome for subtracting sequences with BED.",
                        metavar="FILE", required=True)
    parser.add_argument('-b', '--bed', dest='bed', action='store', nargs='?',
                        help='Genomic regions for subtracting sequences and phasing (BED).',
                        metavar='FILE', required=True)
    parser.add_argument('-t', '--thread', dest='thread', action='store', nargs='?',
                        help='Number of threads (default: 1).',
                        default=1, metavar='STRING')
    parser.add_argument('-f', '--frequency', dest='freq', action='store', nargs='?',
                        help='Minimum frequency of retained alleles (default: 0.05).',
                        default=0.05, metavar='STRING')
    parser.add_argument('-o', '--output', dest='output', action='store', nargs='?',
                        help='Output prefix.',
                        metavar='STRING', required=True)
    args = parser.parse_args()

    BED = pybedtools.BedTool(args.bed)
    results = []
    BED_alleles_dict = {}

    tmp_dir = f'{args.output}_eProbe_temp'
    make_temp_dir(tmp_dir)
    print('Inferring haplotype sequences ...', datetime.datetime.now())
    with multiprocessing.Pool(processes=int(args.thread)) as pool:
        for bed_index in range(len(BED)):
            results.append(pool.apply_async(lanch_process, args=(bed_index, args.bed, args.vcf, args.ref, tmp_dir)))
        for result in results:
            BED_alleles_dict[result.get()[0]] = result.get()[1]

    output_probe_dict = {}
    for bed, alleles in BED_alleles_dict.items():
        alt_count = 1
        for allele, freq in alleles.items():
            if freq >= float(args.freq):
                output_probe_dict[f'{bed}_Allele{str(alt_count)}_Freq{round(freq, 4)}'] = allele
                alt_count = alt_count + 1

    print('Inferring haplotype sequences completed.', datetime.datetime.now())
    write_fasta(output_probe_dict, f'{args.output}.allele.fasta')
    shutil.rmtree(tmp_dir)
