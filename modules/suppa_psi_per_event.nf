process suppa2_psi_per_event {
    //cache false
    tag "psi_per_event"
    publishDir "${params.outdir}/suppa/psi_per_event", mode: 'copy'

    input:
    path ioe_file
    path tpm

    output:
    path "*.psi", emit: psi

    script:
    """
    suppa.py psiPerEvent -i $ioe_file -e $tpm -o ${ioe_file.baseName}
    """
}
