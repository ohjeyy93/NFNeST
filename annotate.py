#!/usr/bin/python3
import os
import sys
import glob
import time
import shutil
import logging
import argparse
import subprocess
import pandas as pd
from pathlib import Path
from itertools import repeat
from multiprocessing import Pool
from nest.bbduk import QualCheck
from nest.alignment import Bwa
from nest.alignment import Bowtie
from nest.alignment import BBMap
from nest.alignment import Snap
from nest.samtools import Samtools
from nest.gatk import GenAnTK
from nest.gatk import Picard
from nest.gatk import FreeBayes
from nest.kestrel import KestrelVar
#from nest.annotater import Annotate
from nest.kestrel import kes_runner
from nest.summarize import Summary
from nest.prepinputs import Prepper
from nest.parsers.vcfReader import Reader 
from nest.parsers.vcfmerge import Merge
from nest.parsers.vcfannotate import Annotate
from nest.parsers.vcfwriter import Writer
def main(arguments):
    ref_path=arguements[0]
    bed_path=arguements[1]
    out_path=arguements[2]
    vcf_path=arguements[3]
    bam_path=arguements[4]
    voi_path=arguements[5]

    #Filer  and annotate variant calls
    main_logger.debug('Annotating variants')
    annotate = Annotate()
    gvcf_path = annotate.getAnnotation(bed_path, gvcf_path, ref_path, out_path, bam_path)
    vcf_path = annotate.getAnnotation(bed_path, vcf_path, ref_path, out_path, bam_path)
    fvcf_path = annotate.getAnnotation(bed_path, fvcf_path, ref_path, out_path, bam_path)
    vcf_dict = {gvcf_path: 'GATK', vcf_path: 'Samtools', fvcf_path: 'Freebayes'}
    merger = Merge(out_path, vcf_dict, ref_path)
    merged_vcf = merger.splitter(list(vcf_dict.keys()))[0]
    final_vcf= '{0}/{1}_variants_merged_annotated.vcf'.format(out_path, sam_name)
    os.rename(merged_vcf, final_vcf)
    #final_path = annotate.getAnnotation(bed_path, final_vcf, ref_path, out_path, bam_path)
    main_logger.debug('Filetering low quality variants and merging GATK and Samtools calls')
    #merged_vcf = Vcf.Merge(gvcf_file, svcf_file, out_path).merge()
    summary = Summary(ref_path, bed_path, voi_path, out_dir)
    var_sum = summary.getVarStats(final_vcf)
    main_logger.info('Total variants : {0}; Verified calls : {1}; Exonic : {2}; Intronic : {3}; Synonymous : {4}; Non Synonymous : {5}; Transition : {6}; Transversion : {7}'.format(
                        var_sum[0], var_sum[1], var_sum[2], var_sum[3], var_sum[4], var_sum[5], var_sum[6], var_sum[7]))
    return(final_vcf, 0)

if __name__ == '__main__':
    #Define deffault paths and aligner informations
    def_path = "{0}/lib".format(os.path.abspath(os.path.dirname(os.path.realpath(__file__))))
    ref_def_path = "{0}/ref".format(os.path.abspath(os.path.dirname(os.path.realpath(__file__))))
    bbduk_def = 'bbduk.sh' #"{0}/bbmap/bbduk.sh".format(def_path)
    bbmap_def = 'bbmap.sh' #"{0}/bbmap/bbmap.sh".format(def_path)
    bwa_def = 'bwa' #"{0}/bwa/bwa".format(def_path)
    bowtie_def = 'bowtie2' #"{0}/bowtie2/bowtie2".format(def_path)
    snap_def = 'snap-alinger' #"{0}/snap/snap-aligner".format(def_path)
    smt_def = 'samtools' #"{0}/samtools/samtools".format(def_path)
    bft_def = 'bcftools' #"{0}/bcftools/bcftools".format(def_path)
    gatk_def = 'gatk' #"{0}/GenomeAnalysisTK.jar".format(def_path)
    pic_def = 'picard' #"{0}/picard.jar".format(def_path)
    sra_def = 'fastq-dump' #'{0}/sratoolkit/bin/fastq-dump'.format(def_path)
    voi_def = None #'{0}/Reportable_SNPs.csv'.format(ref_def_path)
    #if 'java version "1.8.' in str(subprocess.check_output(["java", "-version"], stderr=subprocess.STDOUT).decode('UTF-8').split('\n')[0]):
    java_def = 'java'
    #else:
    #    java_def = "{0}/jdk/bin/java".format(def_path)
    aligner_def = {'bwa' : bwa_def, 'snap' : snap_def, 'bowtie2': bowtie_def, 'bbmap': bbmap_def}
    #Get arguments
    parser = argparse.ArgumentParser(prog='NeST')
    parser.add_argument('-r', '--ref', dest='ref_path', type=str,
                        help='Path to Reference fasta file', required=True)
    parser.add_argument('-b', '--bed', dest='bed_path', type=str,
                        help='Path to Bed file for MDR regions', required=True)
    parser.add_argument('-o', '--outpath', dest='out_path', type=str,
                        help='Path where all outputs will be stored', required=True)
    parser.add_argument('-v', '--vcf', dest='vcf_path', type=str,
                        help='Sample name', required=True)
    parser.add_argument('-m', '--bam', dest='bam_path', type=str, required=True,
                         help='The aligner to used by MARs')    
    parser.add_argument('-voi', '--voi', dest='voi_path', type=str, required=True,
                         help='The aligner to used by MARs')                                     
    args = parser.parse_args()