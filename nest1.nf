#!/usr/bin/env nextflow
params.reads = "$params.input.fastq_path/*/*_insilico{0,1,2,3,4,5,6,7,8,9}{0,1,2,3,4,5,6,7,8,9}{0,1,2,3,4,5,6,7,8,9}{0,1,2,3,4,5,6,7,8,9}{0,1,2,3,4,5,6,7,8,9}{0,1,2,3,4,5,6,7,8,9}{0,1,2,3,4,5,6,7,8,9}_r{1,2}.fastq.gz"
fastq_path = Channel.fromFilePairs(params.reads, size: -1)
fasta_path = Channel.fromPath(params.input.fasta_path)
out_path = Channel.fromPath(params.output.folder)
bed_path = Channel.fromPath(params.input.bed_path)
variant_path = Channel.fromPath(params.input.variant_path)
params.method= "$params.input.method"

process preinputs {
    echo false
    #! /usr/bin/python3
    import os
    #Get FASTQs
    prepper =  Prepper(fastq_path, out_dir, sra_path)
    fastq_path = prepper.sra(sam_name, sra_list, file_list)
    ##Note: Generalize this, right now it will only work with SRA. This is a fix for NEJM
    rone_path = file_list[0]
    rtwo_path = file_list[1]
  
    if not os.path.exists(rone_path):
        raise FileNotFoundException('Forward read not found; Exiting MARs')
        sys.exit()

    if not os.path.exists(rtwo_path):
        raise FileNotFoundException('Reverse read not found; Exiting MARs')
        sys.exit()

    if not os.path.exists(ref_path):
        raise FileNotFoundException('Reference fasta file not found; Exiting MARs')
        sys.exit()

    if not os.path.exists(adp_path):
        raise FileNotFoundException('Adpater sequence not found; Exiting MARs')
        os.mkdir(out_path)

    #Create completion folder
    completion_path = '{0}/completion'.format(out_path)
    if not os.path.exists(completion_path):
        os.mkdir(completion_path)
}

process trimData {
    echo false
    #! /usr/bin/python3
    import os
    #Call Bbduk
    main_logger.debug('Running BBDuk')
    if os.path.exists('{0}/bbduk.rt'.format(completion_path)):
        brone = os.path.splitext(os.path.basename(rone_path))[0]
        brtwo = os.path.splitext(os.path.basename(rtwo_path))[0]
        rone_path = '{0}/{1}/{2}_cleaned.fq'.format(out_path, 'CleanedFastq', brone)
        rtwo_path = '{0}/{1}/{2}_cleaned.fq'.format(out_path, 'CleanedFastq', brtwo)
        main_logger.debug('Skipping BBDuk')
        bret = 0
    else:
        bbduk = QualCheck(bbduk_path, adp_path, out_path, java_path)
        rone_path, rtwo_path, bret = bbduk.bbduk(rone_path, rtwo_path)
        if bret == 0:
            Path('{0}/bbduk.rt'.format(completion_path)).touch()
    if bret != 0:
        raise RuntimeError('BBDuk failed to complete; Exiting MARs')
    else:
        main_logger.debug('BBDuk completed')
}

process alignSequence {
    echo false
    #! /usr/bin/python3
    import os
    if aligner == 'bwa':
        #Call BWA
        main_logger.debug('Running BWA')
        if os.path.exists('{0}/align.rt'.format(completion_path)):
            sam_path = '{0}/alignments/output.sam'.format(out_path)
            mret = 0
            main_logger.debug('Skipping BWA')
        else:
            bwa = Bwa(alinger_path, out_path, ref_path)
            sam_path, mret = bwa.bwamem(rone_path, rtwo_path)
            if mret == 0:
                Path('{0}/align.rt'.format(completion_path)).touch()
        if mret != 0:
            raise RuntimeError('Bwa mem failed to complete; Exiting MARs')
        else:
            main_logger.debug('BWA completed')

    elif aligner == 'bowtie2':
        #Call Bowtie2
        main_logger.debug('Running Bowtie2')
        if os.path.exists('{0}/aling.rt'.format(completion_path)):
            sam_path = '{0}/alignments/output.sam'.format(out_path)
            mret = 0
            main_logger.debug('Skipping Bowtie2')
        else:
            bowtie = Bowtie(alinger_path, out_path, ref_path)
            sam_path, mret = bowtie.bowtie(rone_path, rtwo_path)
            if mret == 0:
                Path('{0}/align.rt'.format(completion_path)).touch()
        if mret != 0:
            raise RuntimeError('Bowtie2 failed to complete; Exiting MARs')
        else:
            main_logger.debug('Bowtie2 completed')

    elif aligner == 'snap':
        #Call Snap
        main_logger.debug('Running Snap')
        snap = Snap(alinger_path, out_path, ref_path)
        sam_path, mret = snap.snap(rone_path, rtwo_path)
        if mret != 0:
            raise RuntimeError('Snap failed to complete; Exiting MARs')
        else:
            main_logger.debug('Snap completed')

    elif aligner == 'bbmap':
        #Call Bbmap
        main_logger.debug('Running BBMap')
        if os.path.exists('{0}/aling.rt'.format(completion_path)):
            sam_path = '{0}/alignments/output.sam'.format(out_path)
            mret = 0
        else:
            bbmap = BBMap(alinger_path, out_path, ref_path)
            sam_path, mret = bbmap.bbmap(rone_path, rtwo_path)
            if mret == 0:
                Path('{0}/align.rt'.format(completion_path)).touch()
        if mret != 0:
            raise RuntimeError('BBMap failed to complete; Exitinign MARs')
        else:
            main_logger.debug('BBMap completed')
}

process fixmate {
    #! /usr/bin/python3
    varengine = Samtools(smt_path, bft_path, out_path)
    if os.path.exists('{0}/fixmate.rt'.format(completion_path)):
        base = os.path.splitext(os.path.basename(sam_path))[0]
        bam_path = '{0}/alignments/{1}_FM.bam'.format(out_path, base)
        fret = 0
        main_logger.debug('Skipping fixmate')
    else:
        bam_path, fret = varengine.fixmate(sam_path)
        if fret == 0:
            Path('{0}/fixmate.rt'.format(completion_path)).touch()
    main_logger.debug('Running Samtools fixmate')
    if fret != 0:
        raise RuntimeError('Samtools fixmate failed to complete; Exiting MARs')
    else:
        main_logger.debug('Samtools fixmate completed')

    if os.path.exists('{0}/sort.rt'.format(completion_path)):
        base = os.path.splitext(os.path.basename(bam_path))[0]
        bam_path = '{0}/alignments/{1}_SR.bam'.format(out_path, base)
        sret = 0
        main_logger.debug('Skipping sort')
    else:
        bam_path, sret = varengine.sort(bam_path)
        if sret == 0:
            Path('{0}/sort.rt'.format(completion_path)).touch()
    main_logger.debug('Running Samtools sort')
    if sret != 0:
        raise RuntimeError('Samtools sort failed to complete; Exiting MARs')
    else:
        main_logger.debug('Samtools sort completed')

    if os.path.exists('{0}/dedup.rt'.format(completion_path)):
        base = os.path.splitext(os.path.basename(bam_path))[0]
        bam_path = '{0}/alignments/{1}_DD.bam'.format(out_path, base)
        dret = 0
        main_logger.debug('Skipping Dedup')
    else:
        bam_path, dret = varengine.dedup(bam_path)
        if dret == 0:
            Path('{0}/dedup.rt'.format(completion_path)).touch()
    main_logger.debug('Running Samtools dedup')
    if sret != 0:
        raise RuntimeError('Samtools dedup failed to complete; Exiting MARs')
    else:
        main_logger.debug('Samtools dedup completed')

    rgadder = Picard(java_path, pic_path, out_path)
    if os.path.exists('{0}/readgroup.rt'.format(completion_path)):
        base = os.path.splitext(os.path.basename(bam_path))[0]
        bam_path = '{0}/alignments/{1}_RG.bam'.format(out_path, base)
        aret = 0
        main_logger.debug('Skipping add read group')
    else:
        bam_path, aret = rgadder.picard(bam_path, sam_name)
        main_logger.debug('Running Picard AddOrReplaceReadGroups')
        if aret == 0:
            Path('{0}/readgroup.rt'.format(completion_path)).touch()
    if aret != 0:
        raise RuntimeError('Picard AddOrReplaceReadGroups failed to complete; Exiting MARs')
    else:
        main_logger.debug('Picard AddOrReplaceReadGroups completed')
}

process generateVCF{
    #! /usr/bin/python3
    if os.path.exists('{0}/pileup.rt'.format(completion_path)):
        bcf_path = '{0}/{1}_variants.bcf'.format(out_path, sam_name)
        pret = 0
        main_logger.debug('Skipping Pileup')
    else:
        bcf_path, pret = varengine.pileup(ref_path, bam_path, sam_name)
        main_logger.debug('Running Samtools mpileup')
        if pret == 0:
            Path('{0}/pileup.rt'.format(completion_path)).touch()
    if pret != 0:
        raise RuntimeError('Samtools mpileup failed to complete; Exiting MARs')
    else:
        main_logger.debug('Samtools mpileup completed')

    if os.path.exists('{0}/bcfindex.rt'.format(completion_path)):
        bret = 0
        main_logger.debug('Skipping Bcfindex')
    else:
        bret = varengine.bcfindex(bcf_path)
        main_logger.debug('Running Bcftools index')
        if bret ==0 :
            Path('{0}/bcfindex.rt'.format(completion_path)).touch()
    if bret != 0:
        raise RuntimeError('Bcftools index failed to complete; Exiting MARs')
    else:
        main_logger.debug('Bcftools index completed')

    if os.path.exists('{0}/bcfcall.rt'.format(completion_path)):
        vcf_path = '{0}/{1}_variants_samtools.vcf'.format(out_path, sam_name)
        bret = 0
        main_logger.debug('Skipping bcfcall')
    else:
        vcf_path, bret = varengine.bcftools(bcf_path, sam_name)
        main_logger.debug('Running Bcftools call')
        if bret == 0:
            Path('{0}/bcfcall.rt'.format(completion_path)).touch()

    if bret != 0:
        raise RuntimeError('Bcftools call failed to complete; Exiting MARs')
    else:
        main_logger.debug('Bcftools call completed')

    #Call GATK HaplotypeCaller to generate VCF files
    varcaller = GenAnTK(gatk_path, out_path, java_path, pic_path)
    main_logger.debug('Running GATK HaplotypeCaller')
    if os.path.exists('{0}/gatk.rt'.format(completion_path)):
        gvcf_path = '{0}/{1}_variants_gatk.vcf'.format(out_path, sam_name)
        gret = 0
        main_logger.debug('Skipping GATK')
    else:
        gvcf_path, gret = varcaller.hapCaller(bam_path, ref_path, sam_name)
        if gret == 0:
            Path('{0}/gatk.rt'.format(completion_path)).touch()
    if gret != 0:
        raise RuntimeError('GATK HaplotypeCaller failed to complete; Exiting MARs')
    else:
        main_logger.debug('GATK HaplotypeCaller stats completed')

    #Call Freebayes to generate VCF files
    varcaller = FreeBayes('freebayes', out_path)
    main_logger.debug('Running Freebayes')
    if os.path.exists('{0}/freebayes.rt'.format(completion_path)):
        fvcf_path = '{0}/{1}_variants_freebayes.vcf'.format(out_path, sam_name)
        fret = 0
        main_logger.debug('Skipping Freebayes')
    else:
        fvcf_path, fret = varcaller.freeBayes(bam_path, ref_path, sam_name)
        if fret == 0:
            Path('{0}/freebayes.rt'.format(completion_path)).touch()
    if fret != 0:
        raise RuntimeError('Freebayes failed to complete; Exiting MARs')
    else:
        main_logger.debug('Freebayes stats completed')
}

process annotate {
    #! /usr/bin/python3
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
    if purge:
       shutil.rmtree('{0}/RawFastq'.format(out_path))
       shutil.rmtree('{0}/CleanedFastq'.format(out_path))
       alignments = glob.glob('{0}/alignments/*'.format(out_path))
       for files in alignments:
           if 'output_FM_SR_DD_RG.ba' in files:
               continue
           else:
               os.remove(files)
       vcffiles = glob.glob('{0}/*.bcf*'.format(out_path))
       for files in vcffiles:
           os.remove(files)
    return(final_vcf, 0)
}

process batch {
    #! /usr/bin/python3
    logger = logging.getLogger('NeST')
    logger.setLevel(logging.DEBUG)
    #Create output paths for the run
    if not os.path.exists(os.path.abspath(out_dir)):
        os.mkdir(os.path.abspath(out_dir))
    # Creating a file handler which logs even debug messages
    fh = logging.FileHandler('{0}/nest.log'.format(os.path.abspath(out_dir)))
    if verbose:
        fh.setLevel(logging.DEBUG)
    else:
        fh.setLevel(logging.INFO)
    # Creating a console handler to log info messages
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    # Create formatter and add it to the handlers
    formatter = logging.Formatter('{asctime} - {name} - {levelname} - {message}', style="{")
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    # Add the handlers to the logger
    logger.addHandler(fh)
    logger.addHandler(ch)
    #Create file and console handlers for MaRS
    logger.info('Gathering input information from input path.')
    prep = Prepper(inp_path, out_dir, sra_path).prepInputs()
    samples, sra_list, files = list(), list(), list()
    logger.info('Running MaRS on {0} experiments'.format(len(prep)))
    #summary = Summary(ref_path, bed_path, voi_path, out_dir)
    #samples = config.keys()
    pools = Pool(threads)
    for sample in prep:
        samples.append(prep[sample].sample)
        files.append(prep[sample].files)
        sra_list.append(prep[sample].sra)
    #rone_list = list()
    #rtwo_list = list()
    #name_list = list()
    #for samples in config:
    #    name_list.append(config[samples].sample)
    #    rone_list.append(config[samples].files[0])
    #    rtwo_list.append(config[samples].files[1])

    #sra_list = files
    vcf_list = pools.map(main, zip(repeat(bbduk_path), repeat(aligner_path),
                repeat(smt_path), repeat(bft_path), repeat(gatk_path),
                samples, files, repeat(ref_path), repeat(adp_path),
                repeat(bed_path), repeat(out_dir), repeat(aligner),
                repeat(pic_path), repeat(voi_path),
                repeat(java_path), repeat(sra_path), repeat(purge), sra_list))
    logger.info('Summarizing variant calls from all {0} experiments'.format(len(prep)))
    summary = Summary(ref_path, bed_path, voi_path, out_dir)
    #Sumarize variants of intrest
    summary.getSummary()
    return(0)
}

#! /usr/bin/python3
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