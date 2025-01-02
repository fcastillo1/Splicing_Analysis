process suppa2_generate_events {
    publishDir "${params.outdir}/suppa/generate_events", mode: 'copy'
    
    input:
    path gtf_file
    
    output:
    path "events_*_strict.ioe", emit: events
    
    script:
    """
    #!/bin/bash
    export PATH=/opt/conda/bin:\$PATH
    which suppa.py
    # Ejecucion de eventos IOE
    suppa.py generateEvents -i ${gtf_file} -o events -f ioe -e ${params.events}
    """
}
