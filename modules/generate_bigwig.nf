process generate_bigwig {
  tag "$sample_id"
  publishDir "${params.outdir}/star/bigwig", mode: 'copy'
  
  input:
  tuple val(sample_id), path(bam), path(bai), val(condition)
  
  output:
  path "${sample_id}.bw", emit: bigwig
  
  script:
  """
  # Proceso que genera archivos bigwig
  bamCoverage -b $bam -o ${sample_id}.bw
  """
}
