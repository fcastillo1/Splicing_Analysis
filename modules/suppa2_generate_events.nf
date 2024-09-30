process suppa2_generate_events {
   // cache false
    publishDir "${params.outdir}/suppa/generate_events", mode: 'copy'

    input:
    path gtf_file

    output:
    path "events_*_strict.ioe", emit: events

    script:
    """
    suppa.py generateEvents -i ${gtf_file} -o events -f ioe -e ${params.events}
    """
}
