#!/usr/bin/env nextflow
params.reads = $params.input.fastq_path+"/NeST_insilico{\d+}{\d+}{\d+}{\d+}{\d+}{\d+}{\d+}_r{1,2}.fastq.gz"
fastq_path = Channel.fromFilePairs(params.reads, size: -1)
fasta_path = Channel.fromPath(params.input.fasta_path)
out_path = Channel.fromPath(params.output.folder)
bed_path = Channel.fromPath(params.input.bed_path)
variant_path = Channel.fromPath(params.input.variant_path)
params.method= "$params.input.method"
bbduk_path = "$params.input.bbduk_def"

process trimData {
    echo false
    publishDir "$params.output.folder/trimFastq/${sample}", pattern: "*.fastq", mode : "copy"
    label 'gizmo'
    scratch "$task.scratch"
    input:
        each samp from fastq_path

    output:
        tuple val(sample), path("${sample}_trimmed_R1.fastq"), path("${sample}_trimmed_R2.fastq"), path("${sample}_*.txt") into trim_out

    script:
        """
        #! /usr/bin/python3
        import os
        from nest.bbduk import QualCheck
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
        """
}