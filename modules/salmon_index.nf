process SALMON_INDEX {
    tag "Salmon index"
    publishDir "${params.outdir}/salmon_index", mode: 'copy'

    input:
    path transcript_fasta
    path genome_fasta

    output:
    path "salmon_index", emit: index
    path "versions.yml", emit: versions

    script:
    """
    # Descomprimir el transcriptoma si está comprimido
    if [[ ${transcript_fasta} == *.gz ]]; then
        gunzip -c ${transcript_fasta} > transcripts.fa
        TRANSCRIPT_FILE=transcripts.fa
    else
        TRANSCRIPT_FILE=${transcript_fasta}
    fi

    # Obtener IDs de decoys
    grep '^>' ${genome_fasta} | cut -d ' ' -f 1 | cut -d \$'\\t' -f 1 > decoys.txt
    sed -i.bak -e 's/>//g' decoys.txt

    # Combinar transcriptoma descomprimido y genoma
    cat \$TRANSCRIPT_FILE ${genome_fasta} > gentrome.fa

    # Crear índice
    salmon index \\
        -t gentrome.fa \\
        -i salmon_index \\
        -d decoys.txt \\
        -p ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        salmon: \$(salmon --version | sed 's/salmon //g')
    END_VERSIONS
    """
}
