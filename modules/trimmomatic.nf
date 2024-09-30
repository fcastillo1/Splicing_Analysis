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
        java -jar /mnt/disco_3/FReyes/Trimmomatic-0.39/trimmomatic-0.39.jar SE \
        $reads \
        ${sample_id}_trimmed.fastq.gz \
        ILLUMINACLIP:${params.adapters}:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
        """
    } else {
        """
        java -jar /mnt/disco_3/FReyes/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
        ${reads[0]} ${reads[1]} \
        ${sample_id}_1_paired_trimmed.fastq.gz ${sample_id}_1_unpaired_trimmed.fastq.gz \
        ${sample_id}_2_paired_trimmed.fastq.gz ${sample_id}_2_unpaired_trimmed.fastq.gz \
        ILLUMINACLIP:${params.adapters}:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
        """
    }
}
