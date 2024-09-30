process SALMON_QUANT {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}/salmon_quant", mode: 'copy'

    input:
    tuple val(meta), path(reads)
    path index
    path gtf
    path transcript_fasta

    output:
    tuple val(meta), path("${meta.id}"), emit: results
    path "versions.yml" , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = meta.id
    def strandedness = 'IU'
    def input_reads = meta.single_end ? "-r $reads" : "-1 ${reads[0]} -2 ${reads[1]}"

    """
    salmon quant \\
        --index $index \\
        --libType=$strandedness \\
        $input_reads \\
        --geneMap $gtf \\
        --threads $task.cpus \\
        --validateMappings \\
        --numBootstraps 100 \\
        --seqBias \\
        --gcBias \\
        $args \\
        -o $prefix

    # Procesar quant.sf
    awk 'BEGIN {OFS="\\t"}
    NR==1 {print \$0}
    NR>1 {split(\$1,a,"|"); print a[1],\$2,\$3,\$4,\$5}' ${prefix}/quant.sf > ${prefix}/quant_temp.sf
    mv ${prefix}/quant_temp.sf ${prefix}/quant.sf

    # Procesar quant.genes.sf
    awk 'BEGIN {OFS="\\t"}
    NR==1 {print \$0}
    NR>1 {split(\$1,a,"|"); print a[1],\$2,\$3,\$4,\$5}' ${prefix}/quant.genes.sf > ${prefix}/quant.genes_temp.sf
    mv ${prefix}/quant.genes_temp.sf ${prefix}/quant.genes.sf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        salmon: \$(salmon --version | sed -e "s/salmon //g")
    END_VERSIONS
    """
}
