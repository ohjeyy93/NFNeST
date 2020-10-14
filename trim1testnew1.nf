#!/usr/bin/env nextflow
BBDUK = file(params.input.bbduk)
params.reads = params.input.fastq_path+'/NeST_insilico*_r{1,2}.fastq'
fastq_path = Channel.fromFilePairs(params.reads, size: -1)
/*fastq_path.view()*/
fasta_path = Channel.fromPath(params.input.fasta_path)
out_path = Channel.fromPath(params.output.folder)
bed_path = Channel.fromPath(params.input.bed_path)
variant_path = Channel.fromPath(params.input.variant_path)
bbduk_path = "$params.input.bbduk_def"
params.genome = "$baseDir/ref/pfalciparum/mdr.fa"
params.fas = "$baseDir/ref/pfalciparum/adapters.fa"

process combineFastq {

    publishDir "$params.output.folder/trimFastq/${pair_id}", pattern: "*.fastq", mode : "copy"
    scratch true
    input:
        set pair_id, path(fastq_group) from fastq_path

    
    output:
        tuple val(pair_id), path("${pair_id}_R1.fastq"), path("${pair_id}_R2.fastq") into comb_out

    script:
        """
        cat *r1.fastq > ${pair_id}_R1.fastq
        cat *r2.fastq > ${pair_id}_R2.fastq
        """
}

comb_out.into{preqc_path; trim_path}

process preFastQC {
    publishDir "$params.output.folder/preFastqQC/${sample}", mode : "copy"
    scratch true
    input:
        set val(sample), path(read_one), path(read_two) from preqc_path
    
    output:
        tuple val(sample), path("*") into preqc_out
    
    script:
        """
        fastqc --extract -f fastq -o ./ -t $task.cpus ${read_one} ${read_two}
        """
}

process trimFastq {
    publishDir "$params.output.folder/trimFastq/${sample}", pattern: "*.fastq", mode : "copy"
    publishDir "$params.output.folder/trimFastq/${sample}/Stats", pattern: "*.txt", mode : "copy"
    scratch true
    input:
        set val(sample), path(read_one), path(read_two) from trim_path
    
    output:
        tuple val(sample), path("${sample}_trimmed_R1.fastq"), path("${sample}_trimmed_R2.fastq"), path("${sample}_*.txt") into trim_out

    script:
        """
        bbduk.sh in=${read_one} in2=${read_two} out=${sample}_trimmed_R1.fastq \\
        out2=${sample}_trimmed_R2.fastq stats=${sample}_stats.txt statscolumns=5 \\
        bhist=${sample}_bhist.txt qhist=${sample}_qhist.txt qchist=${sample}_qchist.txt \\
        aqhist=${sample}_aqhist.txt bqhist=${sample}_bqhist.txt lhist=${sample}_lhist.txt \\
        gchist=${sample}_gchist.txt k=13 t=$task.cpus  \\
        qtrim=rl -Xmx8g
        """
}

trim_out.into{align_path; postqc_path}

process postFastQC {
    echo false
    publishDir "$params.output.folder/postFastqQC/${sample}", mode : "copy"
    input:
        set val(sample), path(read_one), path(read_two), path(stats_path) from postqc_path
    
    output:
        tuple val(sample), path("*") into postqc_out
    
    script:
        """
        fastqc --extract -f fastq -o ./ -t $task.cpus ${read_one} ${read_two}
        """
}

process buildIndex {
    tag "$genome.baseName"
    publishDir "$params.output.folder/Bowtie2Index", mode : "copy"

    input:
    path genome from params.genome

    output:
    file 'genome.index*' into index_ch

    """
    bowtie2-build ${genome} genome.index
    """
}



process alignReads {
    echo false
    publishDir "$params.output.folder/alignedReads/${sample}", pattern: "*.sam", mode : "copy"
    publishDir "$params.output.folder/alignedReads/${sample}/Stats", pattern: "*.txt", mode : "copy"
    input:
        set val(sample), path(read_one), path(read_two), path(stats_path) from align_path
        file index from index_ch
    
    output:
        tuple val(sample), path("${sample}.sam"), path("${sample}_bmet.txt") into align_out

    script:
        index_base = index[0].toString() - ~/.rev.\d.bt2?/ - ~/.\d.bt2?/
        """
        bowtie2 --very-sensitive  --dovetail --met-file ${sample}_bmet.txt -p $task.cpus \\
        -x $baseDir/$params.output.folder/Bowtie2Index/$index_base \\
        -1 ${read_one} -2 ${read_two} -p 4 -S ${sample}.sam 
        """ 
}

process processAlignments {
    echo false
    publishDir "$params.output.folder/alignedReads/${sample}", pattern: "*.ba*", mode : "copy"
    publishDir "$params.output.folder/alignedReads/${sample}/Stats", pattern: "*.txt", mode : "copy"
    input:
        set val(sample), path(sam_path), path(align_met) from align_out
        path fas from params.fas
    
    output:
        tuple val(sample), path("${sample}.bam"), path("${sample}.bai"), path("${sample}_*.txt") into postal_out
        path "${sample}.bam" into config_out

    script:
        """
        java -jar $baseDir/picard.jar SortSam -I ${sam_path} -O ${sample}_SR.bam -SORT_ORDER coordinate --CREATE_INDEX true -R ${fas}
        java -jar $baseDir/picard.jar MarkDuplicates -I ${sample}_SR.bam -M ${sample}_dmetrics.txt -O ${sample}.bam --CREATE_INDEX true -R ${fas}
        java -jar $baseDir/picard.jar CollectAlignmentSummaryMetrics -I ${sample}.bam -O ${sample}_ametrics.txt
        """
    
}