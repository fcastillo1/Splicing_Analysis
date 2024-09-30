process suppa2_generate_isoforms {
    //cache false
    publishDir "${params.outdir}/suppa/generate_isoforms", mode: 'copy'
    
    input:
    path gtf_file
    
    output:
    path "*.ioi", emit: isoforms
    
    script:
    """
    suppa.py generateEvents -i ${gtf_file} -o isoforms -f ioi
    """
}
