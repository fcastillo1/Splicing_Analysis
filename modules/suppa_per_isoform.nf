process suppa_per_isoform {
    //cache false
    publishDir "${params.outdir}/suppa/per_isoform", mode: 'copy'

    input:
    path gtf
    path tpm

    output:
    path "*.psi", emit: psi_isoform
    //path "versions.yml", emit: versions

    script:
    """
    suppa.py psiPerIsoform \
        -g $gtf \
        -e $tpm \
        -o per_isoform

    """
}
