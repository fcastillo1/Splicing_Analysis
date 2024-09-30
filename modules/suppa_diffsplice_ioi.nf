process suppa_diffsplice_ioi {
    
    tag "IOI_diffsplice"
    publishDir "${params.outdir}/suppa/diffsplice_ioi", mode: 'copy'

    input:
    tuple path(ioi_file), path(psi_basalA), path(psi_basalB), path(tpm_basalA), path(tpm_basalB)

    output:
    path "*.dpsi"
    path "*.psivec"
    path "cluster_ioi*"

    script:
    """
    suppa.py diffSplice \
        --method empirical \
        --input $ioi_file \
        --psi ${psi_basalA} ${psi_basalB} \
        --tpm ${tpm_basalA} ${tpm_basalB} \
        --area ${params.area} \
        --lower-bound ${params.lower_bound} \
        --alpha ${params.alpha} \
        --tpm-threshold ${params.tpm_threshold} \
        --nan-threshold ${params.nan_threshold} \
        -pa \
        -gc \
        -o diffsplice_output

    suppa.py clusterEvents \
        --dpsi diffsplice_output.dpsi \
        --psivec diffsplice_output.psivec \
        --sig-threshold ${params.sig_threshold} \
        --dpsi-threshold ${params. dpsi_threshold} \
        --eps ${params.eps} \
        --metric ${params.metric} \
        --min-pts ${params.min_pts} \
        --groups ${params.groups} \
        --clustering ${params.clustering_method} \
        -o cluster_ioi
    """
}
