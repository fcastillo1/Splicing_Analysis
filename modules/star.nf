process star {
    tag "$sample_id"
    publishDir "${params.outdir}/star", mode: 'copy'
    label 'process_high'

    input:
    tuple val(sample_id), path(reads), val(condition)
    path genomeDir
    path gtf

    output:
    tuple val(sample_id), path("*Aligned.out.bam"), val(condition), emit: aligned_reads
    path "*Log.final.out", emit: log_final
    path "*Log.out", emit: log_out
    path "*Log.progress.out", emit: log_progress
    path "*SJ.out.tab", emit: sj_out
    path "*ReadsPerGene.out.tab", emit: reads_per_gene
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${sample_id}"
    def read1 = reads[0]
    def read2 = reads.size() > 1 ? reads[1] : ''
    """
    STAR --genomeDir $genomeDir \\
        --readFilesIn $read1 $read2 \\
        --readFilesCommand zcat \\
        --readMatesLengthsIn NotEqual \\
        --outFileNamePrefix ${prefix}_ \\
        --runThreadN ${task.cpus} \\
        --sjdbGTFfile $gtf \\
        --sjdbOverhang ${params.overhang} \\
        --alignSJDBoverhangMin ${params.sjdbOverhangMin} \\
        --outFilterScoreMinOverLread ${params.filterScore} \\
        --outFilterMatchNminOverLread ${params.filterScore} \\
        --outFilterMismatchNmax ${params.outFilterMismatchNmax} \\
        --outFilterMultimapNmax 20 \\
        --alignMatesGapMax 1000000 \\
        --outSAMattributes All \\
        --outSAMtype BAM Unsorted \\
        --outFilterType BySJout \\
        --twopassMode Basic \\
        --alignEndsType Local \\
        --alignIntronMax ${params.alignIntronMax} \\
        --quantMode GeneCounts

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
    END_VERSIONS
    """
}
