process suppa_per_isoform {
    publishDir "${params.outdir}/suppa/per_isoform", mode: 'copy'
    
    input:
    path gtf
    path tpm
    
    output:
    path "*.psi", emit: psi_isoform
    
    script:
    """
    #!/bin/bash
    export PATH=/opt/conda/bin:\$PATH
    # Ejecucion de suppa para obtener las isoformas
    suppa.py psiPerIsoform -g $gtf -e $tpm -o per_isoform
    """
}
