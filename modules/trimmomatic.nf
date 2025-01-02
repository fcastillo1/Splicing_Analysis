process trimmomatic {
    tag "$sample_id"
    publishDir "${params.outdir}/trimming", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads), val(condition)
    
    output:
    tuple val(sample_id), path("*_trimmed.fastq.gz"), val(condition), emit: trimmed_reads
    path "*_unpaired_trimmed.fastq.gz", optional: true, emit: unpaired_reads
    
    script:
    def single_end = reads instanceof Path
    if (single_end) {
        """
        # Ejecucion del procesamiento de trimmomatic para SE
        java -jar \$TRIMMOMATIC_HOME/trimmomatic-0.39.jar SE \\
            -threads ${task.cpus} \\
            ${reads} \\
            ${sample_id}_trimmed.fastq.gz \\
            ILLUMINACLIP:${params.adapters}:${params.trimmomatic_illuminaclip} \\
            LEADING:${params.trimmomatic_leading} \\
            TRAILING:${params.trimmomatic_trailing} \\
            SLIDINGWINDOW:${params.trimmomatic_slidingwindow} \\
            MINLEN:${params.trimmomatic_minlen}
        """
    } else {
        """
        # Ejecucion del procesamiento de trimmomatic para PE
        java -jar \$TRIMMOMATIC_HOME/trimmomatic-0.39.jar PE \\
            -threads ${task.cpus} \\
            ${reads[0]} ${reads[1]} \\
            ${sample_id}_1_paired_trimmed.fastq.gz ${sample_id}_1_unpaired_trimmed.fastq.gz \\
            ${sample_id}_2_paired_trimmed.fastq.gz ${sample_id}_2_unpaired_trimmed.fastq.gz \\
            ILLUMINACLIP:${params.adapters}:${params.trimmomatic_illuminaclip} \\
            LEADING:${params.trimmomatic_leading} \\
            TRAILING:${params.trimmomatic_trailing} \\
            SLIDINGWINDOW:${params.trimmomatic_slidingwindow} \\
            MINLEN:${params.trimmomatic_minlen}
        """
    }
}
