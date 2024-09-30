process SALMON_INDEX {
    tag "$transcript_fasta"
    label 'process_high'
    publishDir "${params.outdir}/salmon_index", mode: 'copy'

    input:
    path genome_fasta
    path transcript_fasta

    output:
    path "salmon_index"   , emit: index
    path "versions.yml"   , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    grep "^>" <(gunzip -c $genome_fasta) | cut -d " " -f 1 > decoys.txt
    sed -i.bak -e 's/>//g' decoys.txt
    cat $transcript_fasta $genome_fasta > gentrome.fa.gz

    salmon index \
        --threads $task.cpus \
        -t gentrome.fa.gz \
        -d decoys.txt \
        -i salmon_index \
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        salmon: \$(salmon --version | sed -e "s/salmon //g")
    END_VERSIONS
    """
}
