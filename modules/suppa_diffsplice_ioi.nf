process suppa_diffsplice_ioi {
    tag "IOI_diffsplice"
    publishDir "${params.outdir}/suppa/diffsplice_ioi", mode: 'copy'

    input:
    tuple path(ioi_file), path(psi_basalA), path(psi_basalB), path(tpm_basalA), path(tpm_basalB)

    output:
    path "diffsplice_output.dpsi", emit: dpsi
    path "diffsplice_output.psivec", emit: psivec
    path "cluster_ioi*", emit: clusters
    path "ioi_diffsplice.log", emit: log

    script:
    """
    #!/bin/bash
    export PATH=/opt/conda/bin:\$PATH
    
    # Iniciar log
    echo "Iniciando análisis SUPPA diffSplice IOI: \$(date)" > ioi_diffsplice.log
    echo "Verificando archivos de entrada:" >> ioi_diffsplice.log
    ls -l ${ioi_file} ${psi_basalA} ${psi_basalB} ${tpm_basalA} ${tpm_basalB} >> ioi_diffsplice.log

    echo "Ejecutando diffSplice:" >> ioi_diffsplice.log
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
        -o diffsplice_output \
        >> ioi_diffsplice.log 2>&1

    # Verificar la salida de diffSplice
    if [ -s diffsplice_output.dpsi ] && [ -s diffsplice_output.psivec ]; then
        echo "Archivos .dpsi y .psivec generados correctamente" >> ioi_diffsplice.log
        echo "Ejecutando clusterEvents:" >> ioi_diffsplice.log
        
        set +e  # Permitir continuar si hay error en clustering
        suppa.py clusterEvents \
            --dpsi diffsplice_output.dpsi \
            --psivec diffsplice_output.psivec \
            --sig-threshold ${params.sig_threshold} \
            --dpsi-threshold ${params.dpsi_threshold} \
            --eps ${params.eps} \
            --metric ${params.metric} \
            --min-pts ${params.min_pts} \
            --groups ${params.groups} \
            --clustering ${params.clustering_method} \
            -o cluster_ioi \
            >> ioi_diffsplice.log 2>&1
        
        if [ \$? -ne 0 ]; then
            echo "Advertencia: Error en clusterEvents" >> ioi_diffsplice.log
        fi
        set -e
    else
        echo "Error: Archivos de salida de diffSplice vacíos o no generados" >> ioi_diffsplice.log
        exit 1
    fi

    echo "Proceso completado: \$(date)" >> ioi_diffsplice.log
    echo "Archivos generados:" >> ioi_diffsplice.log
    ls -l >> ioi_diffsplice.log
    """
}
