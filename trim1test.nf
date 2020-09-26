#!/usr/bin/env nextflow(
params.reads = "$params.input.fastq_path/NeST_insilico\d+_r{1,2}.fastq"
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

    script:
        """
        print(sample)
        """
}