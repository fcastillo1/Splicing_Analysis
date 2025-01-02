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
    # Ejecucion del alineamiento
    STAR --genomeDir $genomeDir \\
        --readFilesIn $read1 $read2 \\
        --readFilesCommand ${params.star_readFilesCommand} \\
        --readMatesLengthsIn ${params.star_readMatesLengthsIn} \\
        --outFileNamePrefix ${prefix}_ \\
        --runThreadN ${task.cpus} \\
        --sjdbGTFfile $gtf \\
        --sjdbOverhang ${params.star_overhang} \\
        --alignSJoverhangMin ${params.star_alignSJoverhangMin} \\
        --alignSJDBoverhangMin ${params.star_alignSJDBoverhangMin} \\
        --outFilterScoreMinOverLread ${params.star_filterScore} \\
        --outFilterMatchNminOverLread ${params.star_filterScore} \\
        --outFilterMismatchNmax ${params.star_outFilterMismatchNmax} \\
        --outFilterMultimapNmax ${params.star_outFilterMultimapNmax} \\
        --alignMatesGapMax ${params.star_alignMatesGapMax} \\
        --alignIntronMin ${params.star_alignIntronMin} \\
        --alignIntronMax ${params.star_alignIntronMax} \\
        --outSAMattributes ${params.star_outSAMattributes} \\
        --outSAMtype ${params.star_outSAMtype} \\
        --outFilterType ${params.star_outFilterType} \\
        --twopassMode ${params.star_twopassMode} \\
        --alignEndsType ${params.star_alignEndsType} \\
        --quantMode ${params.star_quantMode} \\
        ${params.star_soft_clipping ? '' : '--alignEndsType EndToEnd'} \\
        ${params.star_save_unmapped ? '--outReadsUnmapped Fastx' : ''}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
    END_VERSIONS
    """
}
