#!/usr/bin/env nextflow
BBDUK = file(params.input.bbduk)
params.reads = params.input.fastq_path+'*{R1,R2}*.fastq.gz'
/*params.reads = params.input.fastq_path+'/*_{1,2}.fastq.gz'*/
/*params.reads = params.input.fastq_path+'/18HA*R{1,2}*.fastq.gz'*/
fastq_path = Channel.fromFilePairs(params.reads, size: -1)
/*fastq_path.view()*/
fasta_path = Channel.fromPath(params.input.fasta_path)
out_path = Channel.fromPath(params.output.folder)
bed_path = Channel.fromPath(params.input.bed_path)
variant_path = Channel.fromPath(params.input.variant_path)
bbduk_path = "$params.input.bbduk_def"
params.genome = "$baseDir/ref/pfalciparum/New_6_genes.fa"
params.bed = "$baseDir/ref⁩/pfalciparum⁩/New_6_genes.bed"
params.fas = "$baseDir/ref/pfalciparum/adapters.fa"
params.gatk= "$baseDir/gatk-4.1.9.0/gatk"
mode = "$params.input.mode"
vcfmode = "$params.input.vcfmode"
logpath="$baseDir/.nextflow.log"
workpath="$baseDir/work"
CTexcelpath="$baseDir/ANG_2019_TES_master.xlsx"
datasettype = "$params.input.datasettype"

process combineFastq {
    publishDir "$params.output.folder/trimFastq/${pair_id}", pattern: "*.fastq", mode : "copy"
    input:
        set val(pair_id), path(fastq_group) from fastq_path

    
    output:
        tuple val(pair_id), path("${pair_id}_1.fastq"), path("${pair_id}_2.fastq") into comb_out

    script:
        """
        zless *R1*.fastq.gz > ${pair_id}_1.fastq
        zless *R2*.fastq.gz > ${pair_id}_2.fastq
        """
}

comb_out.into{preqc_path; trim_path}

process trimFastq {
    publishDir "$params.output.folder/trimFastq/${sample}", pattern: "*.fastq", mode : "copy"
    publishDir "$params.output.folder/trimFastq/${sample}/Stats", pattern: "*.txt", mode : "copy"
    errorStrategy 'ignore'

    input:
        set val(sample), path(read_one), path(read_two) from trim_path
        path fas from params.fas
    
    output:
        tuple val(sample), path("${sample}_trimmed_R1.fastq"), path("${sample}_trimmed_R2.fastq"), path("${sample}_*.txt") into trim_out

    script:
        """          
        bbduk.sh -Xmx1g ktrimright=t k=27 hdist=1 edist=0 ref=${fas} \\
        qtrim=rl trimq=35 minlength=150 trimbyoverlap=t minoverlap=24 ordered=t qin=33 in=${read_one} in2=${read_two} \\
        out=${sample}_trimmed_R1.fastq out2=${sample}_trimmed_R2.fastq stats=${sample}_stats.txt 
        """
}

trim_out.into{align_path; postqc_path}

process buildIndex {
    tag "$genome.baseName"
    publishDir "$params.output.folder/Bowtie2Index", mode : "copy"

    input:
    path genome from params.genome

    output:
    file 'genome.index*' into index_ch

    script:
    if( mode == 'Bowtie')
    """
    $baseDir/bowtie2-build ${genome} genome.index
    """
}



process alignReads {
    publishDir "$params.output.folder/alignedReads/${sample}", pattern: "*.sam", mode : "copy"
    publishDir "$params.output.folder/alignedReads/${sample}/Stats", pattern: "*.txt", mode : "copy"
    input:
        set val(sample), path(read_one), path(read_two), path(stats_path) from align_path
        file index from index_ch
        path genome from params.genome
    
    output:
        tuple val(sample), path("${sample}.sam"), path("${sample}_bmet.txt") into align_out
        tuple val(sample), path("${sample}.sam") into sam_out

    script:
        index_base = index[0].toString() - ~/.rev.\d.bt2?/ - ~/.\d.bt2?/
        if( mode == 'Bowtie')
            """
            $baseDir/bowtie2 --very-sensitive  --dovetail --met-file ${sample}_bmet.txt -p $task.cpus \\
            -x $baseDir/$params.output.folder/Bowtie2Index/$index_base \\
            -1 ${read_one} -2 ${read_two} -p 4 -S ${sample}.sam

            """
        else if( mode == 'Bwa' )
            """
            bwa mem -t 4 ${genome} ${read_one} ${read_two} > ${sample}.sam 
            """
        else if( mode == 'BBMap' )
            """
            bbmap ref=${genome} in=${read_one} in2=${read_two} out=${sample}.sam 
            """
        else if( mode == 'Snap' )
            """
            snap paired ${genome} ${read_one} ${read_two} -t 4 -o -sam ${sample}.sam $task.cpus \\
            """
}

process processAlignments {
    publishDir "$params.output.folder/alignedReads/${sample}", pattern: "*.ba*", mode : "copy"
    errorStrategy 'retry'

    input:
        set val(sample), path(sam_path), path(align_met) from align_out
        path genome from params.genome
    
    output:
        tuple val(sample), path("${sample}_SR.bam") into postal_out
        tuple val(sample), path("${sample}_SR.bam") into postal_out2
        tuple val(sample), path("${sample}_SR.bam") into postal_out3
        tuple val(sample), path("${sample}_SR.bam") into postal_out4
        tuple val(sample), path("${sample}_SR.bam") into postal_out5

    script:
        """
        java -jar $baseDir/picard.jar AddOrReplaceReadGroups -I ${sam_path} -O ${sample}_SR.bam -SORT_ORDER coordinate --CREATE_INDEX true \\
        -LB ExomeSeq -DS ExomeSeq -PL Illumina -CN AtlantaGenomeCenter -DT 2016-08-24 -PI null -ID ${sample} \\
        -PG ${sample} -PM ${sample} -SM ${sample} -PU HiSeq2500 > file
        """
    
}

process bbmap {
    publishDir "$params.output.folder/bbmap/${sample}", pattern: "*.fastq", mode : "copy"

    input:
        set val(sample), path(bam_path) from postal_out
    
    output:
        tuple val(sample), path("${sample}_un.fastq") into bbmap_out

    script:
        """          
        reformat.sh in=${bam_path} out=${sample}_un.fastq unmappedonly primaryonly
        """
}

process blast {
    publishDir "$params.output.folder/blast/${sample}", pattern: "*.fastq", mode : "copy"
    errorStrategy 'ignore'

    input:
        set val(sample), path("${sample}_un.fastq") from bbmap_out

    script:
        """          
        $baseDir/ncbi-blast-2.12.0+/bin/blastn -query ${sample}_un.fastq -db nt -remote -out ${sample}.out
        """
}

