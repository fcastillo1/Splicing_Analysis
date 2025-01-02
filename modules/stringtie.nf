process stringtie {
    tag "$sample_id"
    publishDir "${params.outdir}/stringtie/${sample_id}", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam), path(bai), val(condition)
    path gtf
    
    output:
    tuple val(sample_id), path("${sample_id}.gtf"), emit: stringtie_gtf
    tuple val(sample_id), path("${sample_id}_for_DGE.gtf"), emit: stringtie_dge_gtf
    tuple val(sample_id), path(bam), path(bai), val(condition), emit: bam_output
    
    script:
    def strandedness = params.libTypeMap[params.stranded]['stringtie']
    """
    stringtie $bam -G $gtf -o ${sample_id}.gtf $strandedness -a ${params.stringtie_min_anchor} -p $task.cpus
    stringtie $bam -G $gtf -o ${sample_id}_for_DGE.gtf $strandedness -a ${params.stringtie_min_anchor} -e -p $task.cpus
    """
}
