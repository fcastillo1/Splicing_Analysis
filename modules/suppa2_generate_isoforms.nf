process suppa2_generate_isoforms {
    publishDir "${params.outdir}/suppa/generate_isoforms", mode: 'copy'
    
    input:
    path gtf_file
    
    output:
    path "*.ioi", emit: isoforms
    
    script:
    """
    #!/bin/bash
    export PATH=/opt/conda/bin:\$PATH
    # Ejeccion de suppa para obtener los eventos IOI
    suppa.py generateEvents -i ${gtf_file} -o isoforms -f ioi
    """
}
