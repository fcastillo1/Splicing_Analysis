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
  def strandedness = params.stranded ? (params.stranded == 'first-strand' ? '--rf' : '--fr') : ''
  """
  stringtie $bam -G $gtf -o ${sample_id}.gtf $strandedness -a 8 -p $task.cpus
  stringtie $bam -G $gtf -o ${sample_id}_for_DGE.gtf $strandedness -a 8 -e -p $task.cpus
  """
}
