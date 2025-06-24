#!/usr/bin/env python

# Standard library imports
import os
import sys
import re
import logging
import shutil

# Third party imports
import argparse
import subprocess
import datetime
import gzip
import pandas as pd
from tkinter import END
import concurrent.futures


# Local application imports

from misc_virion import (check_create_dir, check_file_exists, extract_read_list, extract_sample_list, execute_subprocess, check_reanalysis, file_to_list, obtain_group_cov_stats,
                         obtain_overal_stats, annotate_snpeff, user_annotation, user_annotation_aa, annotate_pangolin, annotation_to_html, report_samples_html, create_consensus, kraken, mash_screen)

from compare_virion import (
    ddbb_create_intermediate, remove_position_range, revised_df, ddtb_compare, extract_lowcov)


"""
=============================================================
HEADER
=============================================================
Institution: IiSGM
Author: Sergio Buenestado-Serrano (sergio.buenestado@gmail.com), Pedro J. Sola (pedroscampoy@gmail.com)
Version = 0
Created: 21 June 2021

TODO:
    Adapt check_reanalysis
    Check file with multiple arguments
    Check program is installed (dependencies)
================================================================
END_OF_HEADER
================================================================
"""

# COLORS AND AND FORMATTING

"""
http://ozzmaker.com/add-colour-to-text-in-python/
The above ANSI escape code will set the text colour to bright green. The format is;
\033[  Escape code, this is always the same
1 = Style, 1 for normal.
32 = Text colour, 32 for bright green.
40m = Background colour, 40 is for black.
"""

END_FORMATTING = '\033[0m'
WHITE_BG = '\033[0;30;47m'
BOLD = '\033[1m'
UNDERLINE = '\033[4m'
RED = '\033[31m'
GREEN = '\033[32m'
MAGENTA = '\033[35m'
BLUE = '\033[34m'
CYAN = '\033[36m'
YELLOW = '\033[93m'
DIM = '\033[2m'

logger = logging.getLogger()


# ARGUMENTS

def get_arguments():

    parser = argparse.ArgumentParser(
        prog='virion.py', description='Pipeline to call variants (SNVs) with any non model organism. Specialised in SARS-CoV-2')

    input_group = parser.add_argument_group('Input', 'Input parameters')

    input_group.add_argument('-i', '--input', dest='input_dir', metavar='input_directory',
                             type=str, required=True, help='REQUIRED. Input directory containing all fast5 files')

    input_group.add_argument('-o', '--output', type=str, required=True,
                             help='REQUIRED. Output directory to extract all results')

    input_group.add_argument('-t', '--threads', type=int, dest='threads',
                             required=False, default=30, help='Threads to use. Default: 30)')

    basecall_group = parser.add_argument_group('Basecall', 'Basecall parameters')

    basecall_group.add_argument('-s', '--samples', metavar='Samples', type=str,
                             required=False, help='Sample list for conversion from barcode to samples ID')

    basecall_group.add_argument('--basecall', type=str, default='dorado', required=True, help="REQUIRED. Program to use in data preprocessing (basecaller and barcoder)")

    basecall_group.add_argument('--model', type=str, default='~/Scripts/git_repos/Dorado/Models/dna_r10.4.1_e8.2_400bps_hac@v4.3.0', required=False, help='The basecaller model to run. Default dna_r10.4.1_e8.2_400bps_hac@v4.3.0 | dna_r10.4_e8.1_hac.cfg')

    basecall_group.add_argument('-b', '--require_barcodes_both_ends', required=False, action='store_true',
                             help='Require barcodes at both ends. By default it only requires the barcode at one end for the sequences identification')

    basecall_group.add_argument('--barcode_kit', type=str, required=False, default='SQK-RBK110-96',
                             help='Kit of barcodes used [SQK-RBK110-96|EXP-NBD196|SQK-NBD112-24]. Default: SQK-RBK110-96')

    basecall_group.add_argument('-g', '--gpu', dest='gpu', required=False, default=False,
                             action="store_true", help='Specify GPU device: "auto", or "cuda:<device_id>"')

    basecall_group.add_argument('--num_callers', type=int, dest='num_callers',
                             required=False, default=8, help="Number of parallel basecallers. Default: 8")

    basecall_group.add_argument('--chunks', type=int, dest='chunks', required=False,
                             default=1536, help='Maximum chunks per runner. Default: 1536')

    basecall_group.add_argument('--records_per_fastq', type=int, dest='records_per_fastq',
                             required=False, default=0, help='Maximum number of records per fastq')

    basecall_group.add_argument("-rq", "--min_read_quality", type=int, dest="min_read_quality",
                             required=False, default=8, help="Filter on a minimum average read quality score. Default: 8")

    basecall_group.add_argument("--trim", type=str, dest="trim",
                             required=False, default="adapters", help="Specify what to trim. Options are 'none', 'all', and 'adapters'. The default behaviour is to trim all detected adapters, primers, and barcodes. Choose 'adapters' to just trim adapters. The 'none' choice is equivelent to using --no-trim. Note that this only applies to DNA. RNA adapters are always trimmed. Default: adapters")

    basecall_group.add_argument("--headcrop", type=int, dest="headcrop", required=False,
                             default=20, help="Trim n nucleotides from start of read. Default: 20")

    basecall_group.add_argument("--tailcrop", type=int, dest="tailcrop", required=False,
                             default=20, help="Trim n nucleotides from end of read. Default: 20")

    varcal_group = parser.add_argument_group('Varcal', 'Varcal parameters')

    varcal_group.add_argument('-r', '--reference', metavar="reference",
                              type=str, required=True, help='REQUIRED. File to map against')

    varcal_group.add_argument('-a', '--annotation', metavar="annotation",
                              type=str, required=True, help='REQUIRED. GFF3 file to annotate variants')

    varcal_group.add_argument("-sample", "--sample", metavar="sample",
                              type=str, required=False, help="Sample to identify further files")

    varcal_group.add_argument("-L", "--sample_list", type=str, required=False,
                              help="Sample names to analyse only in the file supplied")

    species_group = parser.add_argument_group(
        "Species determination", "Species databases")

    species_group.add_argument("--kraken2", dest="kraken2_db", type=str,
                               default=False, required=False, help="Kraken2 database")

    species_group.add_argument("--mash_db", dest="mash_db", type=str, required=False,
                               default=False, help="MASH NCBI annotation containing bacterial database")

    variant_group = parser.add_argument_group(
        "Variant Calling", "Variant Calling parameters")

    variant_group.add_argument("-f", "--min_allele_frequency", type=int, dest="min_allele", required=False,
                               default=0.2, help="Minimum fraction of observations supporting an alternate allele. Default: 0.2")

    variant_group.add_argument("-q", "--min_base_quality", type=int, dest="min_quality", required=False,
                               default=15, help="Exclude alleles from analysis below threshold. Default: 15")

    variant_group.add_argument("-freq", "--min_frequency", type=int, dest="min_frequency", required=False,
                               default=0.6, help="Minimum fraction of observations to call a base. Default: 0.6")

    variant_group.add_argument("-d", "--min_depth", type=int, dest="min_depth",
                               required=False, default=12, help="Minimum depth to call a base. Default: 12")

    quality_group = parser.add_argument_group(
        'Quality parameters', "Parameters for diferent Quality conditions")

    quality_group.add_argument('-c', '--coverage30', type=int, default=90, required=False,
                               help='Minimum percentage of coverage at 30x to classify as uncovered. Default: 90')

    quality_group.add_argument('-n', '--min_snp', type=int, required=False,
                               default=1, help='SNP number to pass quality threshold')

    annot_group = parser.add_argument_group(
        'Annotation', 'Parameters for variant annotation')

    annot_group.add_argument('-B', '--annot_bed', type=str, default=[],
                             required=False, action='append', help='BED file to annotate')

    annot_group.add_argument('-V', '--annot_vcf', type=str, default=[],
                             required=False, action='append', help='VCF file to annotate')

    annot_group.add_argument('-A', '--annot_aa', type=str, default=[],
                             required=False, action='append', help='Aminoacid file to annotate')

    annot_group.add_argument('-R', '--remove_bed', type=str, default=False,
                             required=False, help='BED file with positions to remove')

    annot_group.add_argument('--snpeff_database', type=str, required=False,
                             default=False, help='snpEFF annotation database')

    annot_group.add_argument('--pangolin', required=False, default=False,  action="store_true",
                             help='Use pangolin for classification of epidemiological lineages of SARS-CoV2')

    compare_group = parser.add_argument_group(
        'Compare', 'parameters for compare_snp')

    compare_group.add_argument('-S', '--only_snp', required=False,
                               action='store_true', help='Use INDELS while comparing')

    compare_group.add_argument("--min_threshold_discard_sample", required=False, type=float,
                               default=0.6, help="Minimum inaccuracies to discard a sample. Default: 0.6")

    compare_group.add_argument("--min_threshold_discard_position", required=False, type=float,
                               default=0.5, help="Minimum inaccuracies to discard a position. Default: 0.5")

    compare_group.add_argument('--distance', default=5, required=False,
                               help='Minimun distance to cluster groups after comparison')

    arguments = parser.parse_args()

    return arguments


### Functions from Guppy_virion.py ###

def pod5_conversion(input_dir, out_pod5):

    """
    https://github.com/nanoporetech/pod5-file-format
    """

    cmd_conversion = "pod5 convert fast5 {}/*.fast5 --output {} --one-to-one {} --threads {}".format(input_dir, out_pod5, input_dir, str(args.threads))

    print(cmd_conversion)
    execute_subprocess(cmd_conversion, isShell=True)


def basecalling_dorado(input_dir, out_basecalling_dir):

    """
    https://github.com/nanoporetech/dorado
    """

    fastq_basecalled = os.path.join(out_basecalling_dir, 'Dorado_basecalled.fastq')

    cmd_basecalling = "dorado basecaller {} {} --trim {} --emit-fastq > {}".format(str(args.model), input_dir, args.trim, fastq_basecalled)

    print(cmd_basecalling)
    execute_subprocess(cmd_basecalling, isShell=True)


def demux_dorado(fastq_basecalled, out_barcoding_dir, require_barcodes_both_ends=False, barcode_kit="EXP-NBD104"):

    if require_barcodes_both_ends:
        logger.info(
            GREEN + BOLD + "Barcodes are being used at both ends" + END_FORMATTING + "\n")
        require_barcodes_both_ends = "--barcode-both-ends"
    else:
        logger.info(
            YELLOW + BOLD + "Barcodes are being used on at least 1 of the ends" + END_FORMATTING + "\n")
        require_barcodes_both_ends = ""

    cmd_demux = "dorado demux {} --kit-name {} --output-dir {} --emit-fastq {} --emit-summary".format(fastq_basecalled, barcode_kit, out_barcoding_dir, require_barcodes_both_ends)

    print(cmd_demux)
    execute_subprocess(cmd_demux, isShell=True)


def reorganize_demux(out_barcoding_dir):

    for filename in os.listdir(out_barcoding_dir):
        if filename.endswith('.fastq'):
            fastq_file_path = os.path.join(out_barcoding_dir, filename)
            # with open(fastq_file_path, 'rb') as f_in:
            #     with gzip.open(os.path.join(out_barcoding_dir, filename + '.gz'), 'wb') as f_out:
            #         shutil.copyfileobj(f_in, f_out)
            # os.remove(fastq_file_path)

            if 'barcode' in filename:
                barcode_folder = os.path.join(out_barcoding_dir, filename.split('_')[2].split('.')[0])
                os.makedirs(barcode_folder, exist_ok=True)
                shutil.move(fastq_file_path, barcode_folder)


def basecalling_ion(input_dir, out_basecalling_dir, config='dna_r9.4.1_450bps_fast.cfg', records=0):

    # -i: Path to input fast5 files
    # -s: Path to save fastq files
    # -c: Config file to use > https://community.nanoporetech.com/posts/guppy-v5-0-7-release-note (fast // hac // sup)
    # -x: Specify GPU device: 'auto', or 'cuda:<device_id>'
    # --num_callers: Number of parallel basecallers to Basecaller, if supplied will form part
    # --gpu_runners_per_device: Number of runners per GPU device.
    # --cpu_threads_per_caller: Number of CPU worker threads per basecaller
    # --chunks_per_runner: Maximum chunks per runner
    # --compress_fastq: Compress fastq output files with gzip
    # --records_per_fastq: Maximum number of records per fastq file, 0 means use a single file (per worker, per run id)

    if args.gpu != False:
        logger.info(
            GREEN + 'Basecalling executing on GPU device' + END_FORMATTING)
        gpu_device = "auto"
        cmd = ['guppy_basecaller', '-i', input_dir, '-s', out_basecalling_dir, '-c',
               config, '-x', gpu_device, '--records_per_fastq', str(records), '--compress_fastq']
    else:
        logger.info(
            YELLOW + 'Basecalling executing on CPU device' + END_FORMATTING)
        cmd = ['guppy_basecaller', '-i', input_dir, '-s', out_basecalling_dir, '-c', config, '--num_callers', str(args.num_callers), '--chunks_per_runner', str(
            args.chunks), '--cpu_threads_per_caller', str(args.threads), '--records_per_fastq', str(records), '--compress_fastq']

    print(cmd)
    execute_subprocess(cmd, isShell=False)


def barcoding_ion(out_basecalling_dir, out_barcoding_dir, require_barcodes_both_ends=False, barcode_kit="EXP-NBD104", threads=30):

    # -i: Path to input files
    # -r: Search for input file recursively
    # -s: Path to save files
    # -t: Number of worker threads
    # --num_barcoding_threads: Number of worker threads to use for barcoding.
    # -x: Specify GPU device to accelerate barcode detection: 'auto', or 'cuda:<device_id>'.

    # --fastq_out: Output Fastq files
    # --compress_fastq: Compress fastq output files with gzip
    # --barcode_kits: Space separated list of barcoding kit(s) or expansion kit(s) to detect against. Must be in double quotes
    # --require_barcodes_both_ends: Reads will only be classified if there is a barcode above the min_score at both ends of the read
    # --records_per_fastq: Maximum number of records per fastq file, 0 means use a single file (per worker, per run id)
    # --allow_inferior_barcodes: Reads will still be classified even if both the barcodes at the front and rear (if applicable) were not the best scoring barcodes above the min_score.

    # --detect_barcodes: Detect barcode sequences at the front and rear of the read.
    # --detect_adapter: Detect adapter sequences at the front and rear of the read.
    # --detect_primer: Detect primer sequences at the front and rear of the read.
    # --enable_trim_barcodes: Enable trimming of barcodes from the sequences in the output files. By default is false, barcodes will not be trimmed.
    # --trim_adapters: Trim the adapters from the sequences in the output files.
    # --trim_primers: Trim the primers from the sequences in the output files.

    # --min_score_barcode_front: Minimum score to consider a front barcode to be a valid barcode alignment (Default: 60).
    # --min_score_barcode_rear: Minimum score to consider a rear barcode to be a valid alignment (and min_score_front will then be used for the front only when this is set).

    if args.gpu != False:
        gpu_device = "-x auto "
    else:
        gpu_device = ""

    if require_barcodes_both_ends:
        logger.info(
            GREEN + BOLD + "Barcodes are being used at both ends" + END_FORMATTING + "\n")
        require_barcodes_both_ends = "--require_barcodes_both_ends"
    else:
        logger.info(
            YELLOW + BOLD + "Barcodes are being used on at least 1 of the ends" + END_FORMATTING + "\n")
        require_barcodes_both_ends = ""

    cmd = ["guppy_barcoder", "-i", out_basecalling_dir, "-s", out_barcoding_dir, "-r", gpu_device, require_barcodes_both_ends,
           "--barcode_kits", barcode_kit, "-t", str(threads), '--detect_barcodes', '--enable_trim_barcodes', '--detect_primer', '--trim_primers', '--detect_adapter', '--trim_adapters', "--fastq_out"]

    print(cmd)
    execute_subprocess(cmd, isShell=False)


def rename_files(output_samples):

    with open(output_samples, "w") as bc_output:
        for fastq_file in sum_files:
            with open(fastq_file, "r") as f_in:
                shutil.copyfileobj(f_in, bc_output)


def ONT_QC_filtering(output_samples, filtered_samples):

    # -c: Write on standard output, keep the original files unchanged
    # -q: Filter on a minimum average read quality score
    # --headcrop: Trim n nucleotides from start of read
    # --tailcrop: Trim n nucleotides from end of read

    cmd_filtering = "chopper -i {} -q {} --minlength {} --headcrop {} --tailcrop {} --threads {} | gzip > {}".format(output_samples, str(args.min_read_quality), str(args.length), str(args.headcrop), str(args.tailcrop), str(args.threads), filtered_samples)

    # print(cmd_filtering)
    execute_subprocess(cmd_filtering, isShell=True)


def ONT_quality(output_samples, out_qc, threads=30):

    # --fastq_rich: Data is in one or more fastq file(s) generated by albacore, MinKNOW or guppy with additional information concerning with channel and time
    # --N50: Show the N50 mark in the read length histogram

    cmd_QC = ["NanoPlot", "--fastq", output_samples,
              "--N50", "-o", out_qc, "-t", str(threads)]

    # print(cmd_QC)
    execute_subprocess(cmd_QC, isShell=False)


### Functions from Varcal_virion.py ###

def minimap2_mapping(filename, filename_bam_out, reference):
    """
    https://github.com/lh3/minimap2
        # Oxford Nanopore genomic reads
        minimap2 -ax map-ont ref.fa ont.fq.gz > aln.sam
    http://www.htslib.org/doc/samtools.html
    """

    # -a: Output in the SAM format
    # -x: Preset (always applied before other options; see minimap2.1 for details) []
    #    - map-pb/map-ont - PacBio CLR/Nanopore vs reference mapping
    #    - map-hifi - PacBio HiFi reads vs reference mapping
    #    - ava-pb/ava-ont - PacBio/Nanopore read overlap
    #    - asm5/asm10/asm20 - asm-to-ref mapping, for ~0.1/1/5% sequence divergence
    #    - splice/splice:hq - long-read/Pacbio-CCS spliced alignment
    #    - sr - genomic short-read mapping
    # -t: Number of threads

    # -b: Output BAM
    # -S: Ignored (input format is auto-detected)
    # -F: Only include reads with none of the FLAGS in INT present
    # --threads: Number of additional threads to use

    cmd_minimap2 = "minimap2 -ax map-ont {} {} | samtools view -bS -F 4 - | samtools sort -o {}".format(
        reference, filename, filename_bam_out)
    # print(cmd_minimap2)
    execute_subprocess(cmd_minimap2, isShell=True)

    cmd_indexing = "samtools", "index", filename_bam_out
    # print(cmd_indexing)
    execute_subprocess(cmd_indexing, isShell=False)


def ivar_variants(reference, input_bam, out_variant_dir, sample, annotation, min_quality=15, min_frequency_threshold=0.2, min_depth=20):
    """
    https://andersen-lab.github.io/ivar/html/manualpage.html
    Usage: samtools mpileup -aa -A -d 0 -B -Q 0 --reference [<reference-fasta] <input.bam> | ivar variants -p <prefix> [-q <min-quality>] [-t <min-frequency-threshold>] [-m <minimum depth>] [-r <reference-fasta>] [-g GFF file]

    Note : samtools mpileup output must be piped into ivar variants
    """

    # -aa: Output absolutely all positions, including unused reference sequences. Note that when used in conjunction with a BED file the -a option may sometimes operate as if -aa was specified if the reference sequence has coverage outside of the region specified in the BED file.
    # -A: Do not skip anomalous read pairs in variant calling. Anomalous read pairs are those marked in the FLAG field as paired in sequencing but without the properly-paired flag set.
    # -d: At a position, read maximally INT reads per input file. Setting this limit reduces the amount of memory and time needed to process regions with very high coverage. Passing zero for this option sets it to the highest possible value, effectively removing the depth limit.
    # -B: Disable base alignment quality (BAQ) computation.
    # -Q: Minimum base quality for a base to be considered.

    # -p: (Required) Prefix for the output tsv variant file.
    # -q: Minimum quality score threshold to count base (Default: 20)
    # -t: Minimum frequency threshold(0 - 1) to call variants (Default: 0.03)
    # -m: Minimum read depth to call variants (Default: 0)
    # -r: Reference file used for alignment. This is used to translate the nucleotide sequences and identify intra host single nucleotide variants.
    # -g: A GFF file in the GFF3 format can be supplied to specify coordinates of open reading frames (ORFs). In absence of GFF file, amino acid translation will not be done.

    ivar_raw = os.path.join(out_variant_dir, "ivar_raw")
    check_create_dir(ivar_raw)
    prefix = ivar_raw + '/' + sample

    cmd_ivar = "samtools mpileup -aa -A -d 0 -B -Q 0 --reference {} {} | ivar variants -p {} -q {} t {} -m {} -r {} -g {}".format(
        reference, input_bam, prefix, str(min_quality), str(min_frequency_threshold), str(min_depth), reference, annotation)

    # print(cmd_ivar)
    execute_subprocess(cmd_ivar, isShell=True)


def filter_tsv_variants(tsv_file, output_filtered, min_frequency=0.6, min_total_depth=10, min_alt_dp=4, is_pass=True, only_snp=False):

    input_file_name = os.path.basename(tsv_file)
    input_file = os.path.abspath(tsv_file)
    output_file = os.path.join(output_filtered, input_file_name)

    df = pd.read_csv(input_file, sep='\t')
    df = df.drop_duplicates(subset=['POS', 'REF', 'ALT'], keep="first")
    filtered_df = df[(df.PASS == is_pass) &
                     (df.TOTAL_DP >= min_total_depth) &
                     (df.ALT_DP >= min_alt_dp) &
                     (df.ALT_FREQ >= min_frequency)]

    if only_snp == True:
        final_df = filtered_df[~(filtered_df.ALT.str.startswith(
            '+') | filtered_df.ALT.str.startswith('-'))]
        final_df.to_csv(output_file, sep='\t', index=False)
    else:
        filtered_df.to_csv(output_file, sep='\t', index=False)


def ivar_consensus(input_bam, output_consensus, sample, min_quality=15, min_frequency_threshold=0.6, min_depth=20, uncovered_character='N'):
    """
    https://andersen-lab.github.io/ivar/html/manualpage.html
    Usage: samtools mpileup -aa -A -d 0 -Q 0 <input.bam> | ivar consensus -p <prefix> 
    Note : samtools mpileup output must be piped into ivar consensus
    """

    # -aa: Output absolutely all positions, including unused reference sequences. Note that when used in conjunction with a BED file the -a option may sometimes operate as if -aa was specified if the reference sequence has coverage outside of the region specified in the BED file.
    # -A: Do not skip anomalous read pairs in variant calling. Anomalous read pairs are those marked in the FLAG field as paired in sequencing but without the properly-paired flag set.
    # -d: At a position, read maximally INT reads per input file. Setting this limit reduces the amount of memory and time needed to process regions with very high coverage. Passing zero for this option sets it to the highest possible value, effectively removing the depth limit.
    # -B: Disable base alignment quality (BAQ) computation.
    # -Q: Minimum base quality for a base to be considered.

    # -p: (Required) Prefix for the output tsv variant file.
    # -q: Minimum quality score threshold to count base (Default: 20)
    # -t: Minimum frequency threshold(0 - 1) to call variants (Default: 0.03)
    # -m: Minimum read depth to call variants (Default: 0)
    # -n: Character to print in regions with less than minimum coverage(Default: N)

    prefix = output_consensus + '/' + sample

    cmd_consensus = "samtools mpileup -aa -A -d 0 -B -Q 0 {} |  ivar consensus -p {} -q {} -t {} -m {} -n {}".format(
        input_bam, prefix, min_quality, min_frequency_threshold, min_depth, uncovered_character)

    # print(cmd_consensus)
    execute_subprocess(cmd_consensus, isShell=True)


def replace_consensus_header(input_fasta):

    with open(input_fasta, 'r+') as f:
        content = f.read()
        header = content.split('\n')[0].strip('>')
        new_header = header.split('_')[1].strip()
        content = content.replace(header, new_header)
        f.seek(0)
        f.write(content)
        f.truncate()


def create_bamstat(input_bam, output_file, threads=36):

    cmd_bamstat = "samtools flagstat --threads {} {} > {}".format(
        str(threads), input_bam, output_file)

    # print(cmd_bamstat)
    execute_subprocess(cmd_bamstat, isShell=True)


def create_coverage(input_bam, output_file):

    cmd_coverage = "samtools depth -aa {} > {}".format(input_bam, output_file)
    # print(cmd_coverage)
    execute_subprocess(cmd_coverage, isShell=True)


######################################################################
########################### START PIPELINE ###########################
######################################################################


if __name__ == "__main__":

    args = get_arguments()

    input_dir = os.path.abspath(args.input_dir)

    output_dir = os.path.abspath(args.output)
    group_name = output_dir.split('/')[-1]
    check_create_dir(output_dir)

    reference = os.path.abspath(args.reference)
    annotation = os.path.abspath(args.annotation)

    # Declare folders created in pipeline and key files

    out_basecalling_dir = os.path.join(input_dir, 'Basecalling')
    check_create_dir(out_basecalling_dir)
    basecalling_summary = os.path.join(
        out_basecalling_dir, 'sequencing_summary.txt')
    out_barcoding_dir = os.path.join(input_dir, 'Barcoding')
    check_create_dir(out_barcoding_dir)
    barcoding_summary = os.path.join(
        out_barcoding_dir, 'barcoding_summary.txt')
    out_samples_dir = os.path.join(input_dir, 'Samples_Fastq')
    check_create_dir(out_samples_dir)
    fastq_basecalled = os.path.join(out_basecalling_dir, 'Dorado_basecalled.fastq')
    # out_samples_filtered_dir = os.path.join(out_samples_dir, "Filtered_Fastq")
    # check_create_dir(out_samples_filtered_dir)
    out_qc_dir = os.path.join(input_dir, 'Quality')
    check_create_dir(out_qc_dir)

    out_species_dir = os.path.join(output_dir, "Species")
    check_create_dir(out_species_dir)

    out_bam_dir = os.path.join(output_dir, "Bam")
    check_create_dir(out_bam_dir)

    out_variant_dir = os.path.join(output_dir, "Variants")
    check_create_dir(out_variant_dir)
    out_variant_ivar_dir = os.path.join(
        out_variant_dir, "ivar_raw")  # subfolder
    check_create_dir(out_variant_ivar_dir)
    out_filtered_ivar_dir = os.path.join(
        out_variant_dir, "ivar_filtered")  # subfolder
    check_create_dir(out_filtered_ivar_dir)

    out_consensus_dir = os.path.join(output_dir, "Consensus")
    check_create_dir(out_consensus_dir)
    out_consensus_ivar_dir = os.path.join(
        out_consensus_dir, "ivar")  # subfolder
    check_create_dir(out_consensus_ivar_dir)

    out_stats_dir = os.path.join(output_dir, "Stats")
    check_create_dir(out_stats_dir)
    out_stats_bamstats_dir = os.path.join(
        out_stats_dir, "Bamstats")  # subfolder
    check_create_dir(out_stats_bamstats_dir)
    out_stats_coverage_dir = os.path.join(
        out_stats_dir, "Coverage")  # subfolder
    check_create_dir(out_stats_coverage_dir)

    out_compare_dir = os.path.join(output_dir, "Compare")
    check_create_dir(out_compare_dir)

    out_annot_dir = os.path.join(output_dir, "Annotation")
    check_create_dir(out_annot_dir)
    out_annot_snpeff_dir = os.path.join(out_annot_dir, "snpeff")  # subfolder
    check_create_dir(out_annot_snpeff_dir)
    out_annot_user_dir = os.path.join(out_annot_dir, "user")  # subfolder
    check_create_dir(out_annot_user_dir)
    out_annot_user_aa_dir = os.path.join(out_annot_dir, "user_aa")  # subfolder
    check_create_dir(out_annot_user_aa_dir)

    # Logging
    # Create log file with date and time

    today = str(datetime.date.today())
    right_now = str(datetime.datetime.now())
    right_now_full = '_'.join(right_now.split(' '))
    log_filename = group_name + '_' + right_now_full + '.log'
    log_folder = os.path.join(output_dir, 'Logs')
    check_create_dir(log_folder)
    log_full_path = os.path.join(log_folder, log_filename)

    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)

    formatter = logging.Formatter('%(asctime)s:%(message)s')

    file_handler = logging.FileHandler(log_full_path)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)

    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)
    # stream_handler.setFormatter(formatter)

    logger.addHandler(stream_handler)
    logger.addHandler(file_handler)

    logger.info(
        "\n" + BLUE + '############### START PROCESSING FAST5 FILES ###############' + END_FORMATTING + "\n")
    logger.info(args)

    # Obtain all fast5 files from folder

    fast5 = extract_read_list(args.input_dir)

    # Check how many files will be analysed

    sample_list = []

    for sample in fast5:
        sample = extract_sample_list(sample)
        sample_list.append(sample)

    logger.info("\n" + CYAN + "{} Samples will be analysed: {}".format(
        len(sample_list), ",".join(sample_list)) + END_FORMATTING)


    ############### START PIPELINE ###############

    # Dorado

    if args.basecall == 'dorado':

        prior = datetime.datetime.now()

        logger.info("\n" + GREEN + "STARTING BASECALLING WITH DORADO" + END_FORMATTING + "\n")

        if os.path.exists(fastq_basecalled):

            logger.info("\n" + YELLOW + BOLD +
                        "Ommiting BASECALLING" + END_FORMATTING + "\n")

        else:

            if any(file.endswith('.fast5') for file in os.listdir(args.input_dir)):

                out_pod5_dir = os.path.join(input_dir, "Pod5")
                check_create_dir(out_pod5_dir)

                pod5_conversion(input_dir, out_pod5_dir)

                basecalling_dorado(out_pod5_dir, out_basecalling_dir)

            elif any(file.endswith('.pod5') for file in os.listdir(args.input_dir)):

                    basecalling_dorado(input_dir, out_basecalling_dir)

            else:
                sys.exit(f"Error: The files in {input_dir} are not in .fast5 or pod5 format")

            after = datetime.datetime.now()
            print(("Done with function basecalling_dorado in: %s" % (after - prior) + "\n"))

    # Guppy

    elif args.basecall == 'guppy':

        prior = datetime.datetime.now()

        logger.info("\n" + GREEN + "STARTING BASECALLING WITH GUPPY" + END_FORMATTING + "\n")

        if any(file.endswith('.fast5') for file in os.listdir(args.input_dir)):

            if os.path.isfile(basecalling_summary):
                logger.info("\n" + YELLOW + BOLD + "Ommiting BASECALLING" + END_FORMATTING + "\n")

            else:
                basecalling_ion(input_dir, out_basecalling_dir, config=args.model, records=args.records_per_fastq)

            for root, _, files in os.walk(out_basecalling_dir):
                for name in files:
                    if name.startswith('guppy_basecaller_log'):
                        log_file = os.path.join(out_basecalling_dir, name)
                        os.remove(log_file)

        else:
            sys.exit(f"Error: The files in {input_dir} are not in .fast5 format")

        after = datetime.datetime.now()
        print(("Done with function basecalling_ion in: %s" % (after - prior) + "\n"))


    # Barcoding

    ## Dorado

    if args.basecall == 'dorado':

        prior = datetime.datetime.now()

        logger.info("\n" + GREEN + BOLD + "STARTING BARCODING" + END_FORMATTING)

        if any(f.endswith('unclassified.fastq') for f in os.listdir(out_barcoding_dir)):
            logger.info("\n" + YELLOW + BOLD +
                        "Ommiting BARCODING/DEMULTIPLEX" + END_FORMATTING + "\n")
        else:
            logger.info("\n" + GREEN +
                        "STARTING BARCODING/DEMULTIPLEX" + END_FORMATTING + "\n")
            demux_dorado(fastq_basecalled, out_barcoding_dir, barcode_kit=args.barcode_kit, require_barcodes_both_ends=args.require_barcodes_both_ends)
            reorganize_demux(out_barcoding_dir)

        after = datetime.datetime.now()
        print(("Done with function demux_dorado & reorganize_demux in: %s" % (after - prior) + "\n"))

    ## Guppy

    elif args.basecall == 'guppy':

        prior = datetime.datetime.now()

        logger.info("\n" + GREEN + BOLD + "STARTING BARCODING" + END_FORMATTING)

        if os.path.isfile(barcoding_summary):
            logger.info("\n" + YELLOW + BOLD +
                        "Ommiting BARCODING/DEMULTIPLEX" + END_FORMATTING + "\n")
        else:
            logger.info("\n" + GREEN +
                        "STARTING BARCODING/DEMULTIPLEX" + END_FORMATTING + "\n")
            barcoding_ion(out_basecalling_dir, out_barcoding_dir, barcode_kit=args.barcode_kit,
                        threads=args.threads, require_barcodes_both_ends=args.require_barcodes_both_ends)

        after = datetime.datetime.now()
        print(("Done with function barcoding_ion in: %s" % (after - prior) + "\n"))


    # Read Filtering

    prior = datetime.datetime.now()

    logger.info("\n" + GREEN + BOLD +
                "STARTING SAMPLE FILTERING" + END_FORMATTING)

    if args.samples == None:
        logger.info(
            '\n' + GREEN + 'Filtering samples' + END_FORMATTING)
        for root, _, files in os.walk(out_barcoding_dir):
            for subdirectory in _:
                if subdirectory.startswith('barcode'):
                    barcode_dir = os.path.join(root, subdirectory)
                    for root2, _, files2 in os.walk(barcode_dir):
                        if len(files2) > 1:
                            barcode_path = root2
                            sample = barcode_path.split('/')[-1]
                            # print(sample)
                            output_samples = os.path.join(
                                out_samples_dir, sample + '.fastq')
                            # print(output_samples)
                            filtered_samples = os.path.join(
                                output_dir, sample + '.fastq.gz')
                            # print(filtered_samples)

                            logger.info('\n' + BLUE + BOLD +
                                        sample + END_FORMATTING)

                            sum_files = []
                            for name in files2:
                                filename = os.path.join(barcode_path, name)
                                # print(filename)
                                sum_files.append(filename)
                            logger.info(MAGENTA + BOLD + "Processing {} files in {}".format(
                                len(sum_files), sample) + END_FORMATTING)

                            if os.path.isfile(filtered_samples):
                                logger.info(
                                    YELLOW + sample + ' sample already renamed and filtered' + END_FORMATTING)
                            else:
                                logger.info(
                                    GREEN + 'Renaming sample ' + sample + END_FORMATTING)
                                rename_files(output_samples)

                                logger.info(
                                    GREEN + 'Filtering sample ' + sample + END_FORMATTING)
                                ONT_QC_filtering(
                                    output_samples, filtered_samples)
                        else:
                            None

    else:
        logger.info(
            '\n' + GREEN + 'Filtering & Renaming' + END_FORMATTING)
        with open(args.samples, 'r') as f:
            for line in f:
                barcode, sample = line.split('\t')
                # print(barcode,sample)
                barcode_path = os.path.join(out_barcoding_dir, barcode)
                # print(barcode_path)
                output_samples = os.path.join(
                    out_samples_dir, sample.strip() + '.fastq')
                # print(output_samples)
                filtered_samples = os.path.join(
                    output_dir, sample.strip() + '.fastq.gz')
                # print(filtered_samples)

                logger.info('\n' + BLUE + BOLD + sample + END_FORMATTING)

                sum_files = []
                for root, _, files in os.walk(barcode_path):
                    for name in files:
                        filename = os.path.join(barcode_path, name)
                        # print(filename)
                        sum_files.append(filename)
                    logger.info(MAGENTA + BOLD + "Processing {} files in {}".format(
                        len(sum_files), sample) + END_FORMATTING)

                    if os.path.isfile(filtered_samples):
                        logger.info(
                            YELLOW + sample + ' sample already renamed and filtered' + END_FORMATTING)
                    else:
                        logger.info(GREEN + 'Renaming sample ' +
                                    sample + END_FORMATTING)
                        rename_files(output_samples)

                        logger.info(GREEN + 'Filtering sample ' +
                                    sample + END_FORMATTING)
                        ONT_QC_filtering(output_samples, filtered_samples)

    after = datetime.datetime.now()
    print(('\n' + "Done with function rename_files & ONT_QC_filtering in: %s" %
           (after - prior) + "\n"))

    # Quality Check

    prior = datetime.datetime.now()

    logger.info("\n" + GREEN + BOLD +
                "QUALITY CHECK IN RAW" + END_FORMATTING + '\n')

    for root, _, files in os.walk(output_dir):
        for name in files:
            if name.endswith('.fastq.gz'):
                filtered_sample = os.path.join(root, name)
                # print(filtered_sample)
                out_qc = os.path.join(
                    out_qc_dir, os.path.basename(filtered_sample.split(".")[0]))
                check_create_dir(out_qc)
                # print(out_qc)
                report = [x for x in os.listdir(
                    out_qc) if "NanoPlot-report" in x]
                report_file = os.path.join(out_qc, "".join(report))
                # print(report)

                if os.path.isfile(report_file):
                    logger.info(YELLOW + report_file +
                                " EXIST\n" + BOLD + "Ommiting QC for sample " + name + END_FORMATTING)
                else:
                    logger.info(GREEN + "Checking quality in sample " +
                                name + END_FORMATTING)
                    ONT_quality(filtered_sample, out_qc, threads=args.threads)

    after = datetime.datetime.now()
    print(("\n" + "Done with function ONT_quality in: %s" % (after - prior) + "\n"))

    logger.info("\n" + MAGENTA + BOLD +
                "#####END OF ONT DATA PROCESSING #####" + END_FORMATTING + "\n")

    logger.info(
        "\n" + BLUE + "############### START VARIANT CALLING ###############" + END_FORMATTING + "\n")

    logger.info(args)

    # Obtain all fastq files from folder

    fastq = extract_read_list(output_dir)

    # Check how many files will be analysed

    sample_list = []

    for sample in fastq:
        sample = extract_sample_list(sample)
        # sample = sample.split('_')[1]
        sample_list.append(sample)

    # Check if there are samples to filter out

    sample_list_F = []
    if args.sample_list == None:
        logger.info("\n" + "No samples to filter" + "\n")
        for sample in fastq:
            sample = extract_sample_list(sample)
            sample_list_F.append(sample)
    else:
        logger.info("Samples will be filtered")
        sample_list_F = file_to_list(args.sample_list)

    new_samples = check_reanalysis(args.output, sample_list_F)

    logger.info(CYAN + "\n%d samples will be analysed: %s" %
                (len(sample_list_F), ",".join(sample_list_F)) + END_FORMATTING + '\n')

    logger.info(CYAN + "\n%d NEW samples will be analysed: %s" %
                (len(new_samples), ",".join(new_samples)) + END_FORMATTING + '\n')

    new_sample_number = 0

    for sample in fastq:
        # Extract sample name
        sample = extract_sample_list(sample)
        args.sample = sample

        if sample in sample_list_F:
            sample_number = str(sample_list_F.index(sample) + 1)
            sample_total = str(len(sample_list_F))

            if sample in new_samples:
                new_sample_number = str(int(new_sample_number) + 1)
                new_sample_total = str(len(new_samples))
                logger.info("\n" + WHITE_BG + "STARTING SAMPLE: " + sample + " (" + sample_number + "/" +
                            sample_total + ")" + " (" + new_sample_number + "/" + new_sample_total + ")" + END_FORMATTING + '\n')
            else:
                logger.info("\n" + WHITE_BG + "STARTING SAMPLE: " + sample +
                            " (" + sample_number + "/" + sample_total + ")" + END_FORMATTING + '\n')

            HQ_filename = os.path.join(output_dir, sample + ".fastq.gz")
            # print(HQ_filename)
            # filename_out = sample.split('.')[0].split('_')[1]
            filename_out = sample
            # print(filename_out)

            ##### SPECIES DETERMINATION #####

            prior = datetime.datetime.now()

            # Species determination with kraken2 and its standard database and visualization with ktImportTaxonomy from kronatools kit

            sample_species_dir = os.path.join(out_species_dir, sample)
            # print(sample_species_dir)
            check_create_dir(sample_species_dir)
            report = os.path.join(sample_species_dir, sample)
            krona_html = os.path.join(report + ".html")
            mash_output = os.path.join(report + ".screen.tab")

            if args.kraken2_db != False:
                if os.path.isfile(krona_html):
                    logger.info(
                        YELLOW + krona_html + " EXIST\nOmmiting species determination with Kraken2 for " + sample + END_FORMATTING)
                else:
                    logger.info(
                        GREEN + "Species determination with Kraken2 for sample " + sample + END_FORMATTING)
                    kraken(HQ_filename, report, args.kraken2_db,
                           krona_html, threads=args.threads)
            else:
                logger.info(
                    YELLOW + BOLD + "No Kraken database suplied, skipping specie assignation in group " + group_name + END_FORMATTING)

            # Species determination with mash and its bacterial database

            if args.mash_db != False:
                if os.path.isfile(mash_output):
                    logger.info(
                        YELLOW + mash_output + " EXIST\nOmmiting species determination with Mash screen for " + sample + END_FORMATTING)
                else:
                    logger.info(
                        GREEN + "Species determination with Mash for sample " + sample + END_FORMATTING)

                    mash_screen(HQ_filename, mash_output,
                                args.mash_db, winner=True, threads=args.threads)

                    # Name the columns of the mash output and sort them in descending order by identity
                    output_sort_species = pd.read_csv(mash_output, sep='\t', header=None, names=[
                                                      'Identity', 'Share-hashes', 'Median-multiplicity', 'p-value', 'ID accession', 'Organism']).sort_values(by=['Identity'], ascending=False)
                    output_sort_species.to_csv(
                        mash_output, sep='\t', index=None)
            else:
                logger.info(
                    YELLOW + BOLD + "No MASH database suplied, skipping specie assignation in group " + group_name + END_FORMATTING)

            after = datetime.datetime.now()
            print(("Done with function kraken & mash_screen in: %s" %
                   (after - prior) + "\n"))

            ##### MAPPING #####

            # Mapping with minimap2, sorting Bam and indexing it (also can be made with bwa index & bwa mem -x ont2d)

            filename_bam_out = os.path.join(
                out_bam_dir, filename_out + '.sort.bam')
            filename_bam_trimmed = os.path.join(
                out_bam_dir, filename_out + '.trimmed.sort.bam')
            filename_bai_out = os.path.join(
                out_bam_dir, filename_out + '.sort.bam.bai')
            filename_bai_trimmed = os.path.join(
                out_bam_dir, filename_out + '.trimmed.sort.bam.bai')
            # print(filename_bai_out)

            logger.info(GREEN + BOLD +
                        'STARTING ANALYSIS FOR SAMPLE ' + filename_out + END_FORMATTING + '\n')

            prior = datetime.datetime.now()

            if os.path.isfile(filename_bai_out):
                logger.info(YELLOW + filename_bam_out +
                            ' EXIST\nOmmiting mapping for ' + filename_out + END_FORMATTING)
            else:
                logger.info(GREEN + 'Mapping sample ' +
                            filename_out + END_FORMATTING)
                minimap2_mapping(HQ_filename, filename_bam_out,
                                 reference=args.reference)

            after = datetime.datetime.now()
            print(('Done with function minimap2_mapping in: %s' %
                   (after - prior) + '\n'))

            ##### VARIANT CALLING #####

            # Variant calling with samtools mpileup & ivar variants (also can be made with nanopolish, we should use nanopolish index & nanopolish variants)

            prior = datetime.datetime.now()

            out_ivar_variant_name = filename_out + '.tsv'
            out_ivar_variant_file = os.path.join(
                out_variant_ivar_dir, out_ivar_variant_name)

            if os.path.isfile(out_ivar_variant_file):
                logger.info(YELLOW + out_ivar_variant_file +
                            " EXIST\nOmmiting variant call for  sample " + sample + END_FORMATTING)
            else:
                logger.info(
                    GREEN + "Calling variants with ivar in sample " + sample + END_FORMATTING)
                ivar_variants(reference, filename_bam_out, out_variant_dir, sample,
                              annotation, min_quality=args.min_quality, min_frequency_threshold=args.min_allele, min_depth=1)

            after = datetime.datetime.now()
            print(('Done with function ivar_variants in: %s' %
                   (after - prior) + '\n'))

            # Variant filtering by a frequency threshold

            prior = datetime.datetime.now()

            out_ivar_filtered_file = os.path.join(
                out_filtered_ivar_dir, out_ivar_variant_name)

            if os.path.isfile(out_ivar_filtered_file):
                logger.info(YELLOW + out_ivar_filtered_file +
                            " EXIST\nOmmiting variant filtering for  sample " + sample + END_FORMATTING)
            else:
                logger.info(GREEN + 'Filtering variants in sample ' +
                            sample + END_FORMATTING)
                filter_tsv_variants(out_ivar_variant_file, out_filtered_ivar_dir, min_frequency=args.min_frequency,
                                    min_total_depth=args.min_depth, min_alt_dp=4, is_pass=True, only_snp=False)

            after = datetime.datetime.now()
            print(('Done with function filter_tsv_variants in: %s' %
                   (after - prior) + '\n'))

            ##### CONSENSUS #####

            prior = datetime.datetime.now()

            out_ivar_consensus_name = sample + ".fa"
            out_ivar_consensus_file = os.path.join(
                out_consensus_ivar_dir, out_ivar_consensus_name)

            if os.path.isfile(out_ivar_consensus_file):
                logger.info(YELLOW + out_ivar_consensus_file +
                            " EXIST\nOmmiting consensus for sample " + sample + END_FORMATTING)
            else:
                logger.info(
                    GREEN + "Creating consensus with ivar in sample " + sample + END_FORMATTING)
                ivar_consensus(filename_bam_out, out_consensus_ivar_dir, sample, min_quality=args.min_quality,
                               min_frequency_threshold=args.min_frequency, min_depth=args.min_depth, uncovered_character='N')

                logger.info(
                    GREEN + "Replacing consensus header in " + sample + END_FORMATTING)
                replace_consensus_header(out_ivar_consensus_file)

            for root, _, files in os.walk(out_consensus_ivar_dir):
                for name in files:
                    if name.endswith('qual.txt'):
                        qual_consensus = os.path.join(
                            out_consensus_ivar_dir, name)
                        os.remove(qual_consensus)

            after = datetime.datetime.now()
            print(('Done with function ivar_consensus & replace_consensus_header in: %s' %
                   (after - prior) + '\n'))

        ##### CREATE STATS AND QUALITY FILTERS #####

        # Create Bamstats

        prior = datetime.datetime.now()

        out_bamstats_name = sample + ".bamstats"
        out_bamstats_file = os.path.join(
            out_stats_bamstats_dir, out_bamstats_name)

        if os.path.isfile(out_bamstats_file):
            logger.info(YELLOW + out_bamstats_file +
                        " EXIST\nOmmiting Bamstats for sample " + sample + END_FORMATTING)
        else:
            logger.info(GREEN + "Creating bamstats in sample " +
                        sample + END_FORMATTING)
            create_bamstat(filename_bam_out, out_bamstats_file,
                           threads=args.threads)

        after = datetime.datetime.now()
        print(("Done with function create_bamstat in: %s" %
               (after - prior) + "\n"))

        # Create Coverage

        prior = datetime.datetime.now()

        out_coverage_name = sample + ".cov"
        out_coverage_file = os.path.join(
            out_stats_coverage_dir, out_coverage_name)

        if os.path.isfile(out_coverage_file):
            logger.info(YELLOW + out_coverage_file +
                        " EXIST\nOmmiting Coverage for " + filename_out + END_FORMATTING)
        else:
            logger.info(GREEN + "Creating Coverage in sample " +
                        filename_out + END_FORMATTING)
            create_coverage(filename_bam_out, out_coverage_file)

        after = datetime.datetime.now()
        print(("Done with function create_coverage in: %s" %
               (after - prior) + "\n"))

    # Coverage output summary

    prior = datetime.datetime.now()

    logger.info(GREEN + BOLD + "Creating summary report for coverage results in group " +
                group_name + END_FORMATTING)

    obtain_group_cov_stats(out_stats_dir, group_name)

    # Reads and Variants output summary

    logger.info(GREEN + BOLD + "Creating overal summary report in group " +
                group_name + END_FORMATTING)

    obtain_overal_stats(out_stats_dir, output_dir, group_name)

    after = datetime.datetime.now()
    print(("Done with function obtain_group_cov_stats & obtain_overal_stats in: %s" % (
        after - prior) + "\n"))

    ##### ANNOTATION #####

    logger.info('\n' + GREEN + BOLD + 'STARTING ANNOTATION IN GROUP: ' +
                group_name + END_FORMATTING + '\n')

    # Annotation with SnpEFF

    prior = datetime.datetime.now()

    if args.snpeff_database != False:
        # Change for raw/filtered annotation
        for root, _, files in os.walk(out_filtered_ivar_dir):
            if root == out_filtered_ivar_dir:  # Change for raw/filtered annotation
                for name in files:
                    if name.endswith('.tsv'):
                        sample = name.split('.')[0]
                        filename = os.path.join(root, name)
                        out_annot_file = os.path.join(
                            out_annot_snpeff_dir, sample + '.annot')

                        if os.path.isfile(out_annot_file):
                            logger.info(
                                YELLOW + out_annot_file + ' EXIST\nOmmiting snpEff annotation for sample ' + sample + END_FORMATTING)
                        else:
                            logger.info(
                                GREEN + 'Annotation sample with snpEff: ' + sample + END_FORMATTING)
                            output_vcf = os.path.join(
                                out_annot_snpeff_dir, sample + '.vcf')
                            annotate_snpeff(
                                filename, output_vcf, out_annot_file, database=args.snpeff_database)

    else:
        logger.info(YELLOW + BOLD + 'No SnpEFF database suplied, skipping annotation in group ' +
                    group_name + END_FORMATTING)

    after = datetime.datetime.now()
    print(("Done with function annotate_snpeff in: %s" % (after - prior) + "\n"))

    # Annotation for user defined (bed & vcf annot)

    prior = datetime.datetime.now()

    if not args.annot_bed and not args.annot_vcf:
        logger.info(
            YELLOW + BOLD + 'Ommiting user annotation, no BED or VCF files supplied' + END_FORMATTING)
    else:
        # Change for raw/filtered annotation
        for root, _, files in os.walk(out_filtered_ivar_dir):
            if root == out_filtered_ivar_dir:  # Change for raw/filtered annotation
                for name in files:
                    if name.endswith('.tsv'):
                        sample = name.split('.')[0]
                        logger.info(
                            GREEN + 'User bed/vcf annotation in sample {}'.format(sample) + END_FORMATTING)
                        filename = os.path.join(root, name)
                        out_annot_file = os.path.join(
                            out_annot_user_dir, sample + '.tsv')
                        user_annotation(
                            filename, out_annot_file, vcf_files=args.annot_vcf, bed_files=args.annot_bed)

    after = datetime.datetime.now()
    print(("Done with function user_annotation in: %s" % (after - prior) + "\n"))

    # Annotation for user aa defined (aminoacid annot)

    prior = datetime.datetime.now()

    if not args.annot_aa:
        logger.info(
            YELLOW + BOLD + 'Ommiting user aa annotation, no AA files supplied' + END_FORMATTING)
    else:
        for root, _, files in os.walk(out_annot_snpeff_dir):
            if root == out_annot_snpeff_dir:
                for name in files:
                    if name.endswith('.annot'):
                        sample = name.split('.')[0]
                        logger.info(
                            GREEN + '\n' + 'User aa annotation in sample {}'.format(sample) + END_FORMATTING)
                        filename = os.path.join(root, name)
                        out_annot_aa_file = os.path.join(
                            out_annot_user_aa_dir, sample + '.tsv')
                        if os.path.isfile(out_annot_aa_file):
                            user_annotation_aa(
                                out_annot_aa_file, out_annot_aa_file, aa_files=args.annot_aa)
                        else:
                            user_annotation_aa(
                                filename, out_annot_aa_file, aa_files=args.annot_aa)

    after = datetime.datetime.now()
    print(("Done with function user_annotation_aa in: %s" % (after - prior) + "\n"))

    # Pangolin

    prior = datetime.datetime.now()

    with concurrent.futures.ThreadPoolExecutor(max_workers=args.threads) as executor:
        futures_pangolin = []

        if args.pangolin != False:

            out_annot_pangolin_dir = os.path.join(
                out_annot_dir, "pangolin")  # subfolder
            check_create_dir(out_annot_pangolin_dir)

            for root, _, files in os.walk(out_consensus_ivar_dir):
                if root == out_consensus_ivar_dir:
                    for name in files:
                        if name.endswith('.fa'):
                            sample = name.split('.')[0]
                            filename = os.path.join(root, name)
                            out_pangolin_filename = sample + '.lineage.csv'
                            out_pangolin_file = os.path.join(
                                out_annot_pangolin_dir, out_pangolin_filename)

                            if os.path.isfile(out_pangolin_file):
                                logger.info(
                                    YELLOW + out_pangolin_file + ' EXIST\nOmmiting lineage for sample ' + sample + END_FORMATTING)
                            else:
                                logger.info(
                                    GREEN + 'Obtaining lineage in sample ' + sample + END_FORMATTING)
                                future = executor.submit(
                                    annotate_pangolin, filename, out_annot_pangolin_dir, out_pangolin_filename, threads=args.threads, max_ambig=0.6)
                                futures_pangolin.append(future)

                    for future in concurrent.futures.as_completed(futures_pangolin):
                        logger.info(future.result())

        else:
            logger.info(YELLOW + BOLD + 'No Pangolin flag added, skipping lineage assignation in group ' +
                        group_name + END_FORMATTING)

    after = datetime.datetime.now()
    print(("Done with function annotate_pangolin in: %s" % (after - prior) + "\n"))

    # # User AA to html

    # prior = datetime.datetime.now()

    # annotated_samples = []

    # logger.info(
    #     GREEN + 'Adapting annotation to html in {}'.format(group_name) + END_FORMATTING)

    # for root, _, files in os.walk(out_annot_user_aa_dir):
    #     if root == out_annot_user_aa_dir:
    #         for name in files:
    #             if name.endswith('.tsv'):
    #                 sample = name.split('.')[0]
    #                 annotated_samples.append(sample)
    #                 filename = os.path.join(root, name)
    #                 annotation_to_html(filename, sample)

    # annotated_samples = [str(x) for x in annotated_samples]
    # report_samples_html_all = report_samples_html.replace(
    #     'ALLSAMPLES', ('","').join(annotated_samples))  # NEW

    # with open(os.path.join(out_annot_user_aa_dir, 'All_samples.html'), 'w+') as f:
    #     f.write(report_samples_html_all)

    # after = datetime.datetime.now()
    # print(("Done with function annotation_to_html in: %s" % (after - prior) + "\n"))

    ##### COMPARISON #####

    # SNPs comparison using tsv variant files

    prior = datetime.datetime.now()

    logger.info('\n' + GREEN + BOLD + 'STARTING COMPARISON IN GROUP: ' +
                group_name + END_FORMATTING + '\n')

    folder_compare = today + '_' + group_name
    path_compare = os.path.join(out_compare_dir, folder_compare)
    check_create_dir(path_compare)
    full_path_compare = os.path.join(path_compare, group_name)

    compare_snp_matrix_recal = full_path_compare + '.revised.final.tsv'
    compare_snp_matrix_INDEL = full_path_compare + ".revised_INDEL.final.tsv"
    compare_snp_matrix_recal_intermediate = full_path_compare + ".revised_intermediate.tsv"
    compare_snp_matrix_INDEL_intermediate = full_path_compare + \
        ".revised_INDEL_intermediate.tsv"

    recalibrated_snp_matrix_intermediate = ddbb_create_intermediate(
        out_variant_ivar_dir, out_stats_coverage_dir, min_freq_discard=args.min_allele, min_alt_dp=10, only_snp=args.only_snp)
    recalibrated_snp_matrix_intermediate.to_csv(
        compare_snp_matrix_recal_intermediate, sep="\t", index=False)

    after = datetime.datetime.now()
    print(("Done with function ddbb_create_intermediate in: %s" %
           (after - prior) + "\n"))

    prior = datetime.datetime.now()

    compare_snp_matrix_INDEL_intermediate_df = remove_position_range(
        recalibrated_snp_matrix_intermediate)
    compare_snp_matrix_INDEL_intermediate_df.to_csv(
        compare_snp_matrix_INDEL_intermediate, sep="\t", index=False)

    after = datetime.datetime.now()
    print(("Done with function remove_position_range in: %s" %
           (after - prior) + "\n"))

    prior = datetime.datetime.now()

    # Extract all low coverage o not covered positions
    symbol_file = full_path_compare + '_symbol_lowcov.tsv'
    symbol_lowcov = extract_lowcov(compare_snp_matrix_recal_intermediate) # It is made by the INDEL_intermediate.tsv, taking 0 and 1 into account, can also be made with intermediate.highfreq.tsv
    symbol_lowcov.to_csv(symbol_file, sep='\t', index=False)

    recalibrated_revised_df = revised_df(recalibrated_snp_matrix_intermediate, path_compare, min_freq_include=args.min_frequency,
                                         min_threshold_discard_sample=args.min_threshold_discard_sample, min_threshold_discard_position=args.min_threshold_discard_position, remove_faulty=True, drop_samples=True, drop_positions=True)
    # recalibrated_revised_df.to_csv(compare_snp_matrix_recal, sep="\t", index=False)

    recalibrated_revised_df = recalibrated_revised_df[~recalibrated_revised_df['Position'].isin(symbol_lowcov['Position'])]
    recalibrated_revised_df.to_csv(compare_snp_matrix_recal, sep='\t', index=False)

    # Extract all low coverage o not covered positions
    symbol_INDEL_file = full_path_compare + '_symbol_INDEL_lowcov.tsv'
    symbol_INDEL_lowcov = extract_lowcov(compare_snp_matrix_INDEL_intermediate) # It is made by the INDEL_intermediate.tsv, taking 0 and 1 into account, can also be made with intermediate.highfreq.tsv
    symbol_INDEL_lowcov.to_csv(symbol_INDEL_file, sep='\t', index=False)

    recalibrated_revised_INDEL_df = revised_df(compare_snp_matrix_INDEL_intermediate_df, path_compare, min_freq_include=args.min_frequency,
                                               min_threshold_discard_sample=args.min_threshold_discard_sample, min_threshold_discard_position=args.min_threshold_discard_position, remove_faulty=True, drop_samples=True, drop_positions=True)
    # recalibrated_revised_INDEL_df.to_csv(compare_snp_matrix_INDEL, sep="\t", index=False)

    recalibrated_revised_INDEL_df = recalibrated_revised_INDEL_df[~recalibrated_revised_INDEL_df['Position'].isin(symbol_INDEL_lowcov['Position'])]
    recalibrated_revised_INDEL_df.to_csv(compare_snp_matrix_INDEL, sep='\t', index=False)

    after = datetime.datetime.now()
    print(("Done with function revised_df in: %s" % (after - prior) + "\n"))

    prior = datetime.datetime.now()

    ddtb_compare(compare_snp_matrix_recal, distance=args.distance)
    ddtb_compare(compare_snp_matrix_INDEL, distance=args.distance, indel=True)

    after = datetime.datetime.now()
    print(("Done with function ddtb_compare in: %s" % (after - prior) + "\n"))

    logger.info(MAGENTA + BOLD + 'COMPARISON FINISHED IN GROUP: ' +
                group_name + END_FORMATTING + '\n')

    ##### REFINING CONSENSUS #####

    prior = datetime.datetime.now()

    logger.info('\n' + GREEN + BOLD +
                'CREATING REFINED CONSENSUS' + END_FORMATTING)

    create_consensus(reference, compare_snp_matrix_recal,
                     out_stats_coverage_dir, out_consensus_dir)

    after = datetime.datetime.now()
    print(("Done with function create_consensus in: %s" % (after - prior) + "\n"))

    logger.info("\n" + MAGENTA + BOLD +
                "##### END OF ONT VARIANT CALLING PIPELINE #####" + "\n" + END_FORMATTING)
