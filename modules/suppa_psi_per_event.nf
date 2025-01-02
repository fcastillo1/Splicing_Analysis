process suppa2_psi_per_event {
    tag { ioe_file.baseName }
    publishDir "${params.outdir}/suppa/psi_per_event", mode: 'copy'
    
    input:
    path ioe_file
    path tpm
    
    output:
    path "*.psi", emit: psi
    
    script:
    """
    #!/bin/bash
    export PATH=/opt/conda/bin:\$PATH
    # Ejecucion de suppa para obtener el valor de PSI
    suppa.py psiPerEvent -i $ioe_file -e $tpm -o ${ioe_file.baseName}
    """
}
