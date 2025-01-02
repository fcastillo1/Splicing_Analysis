process suppa_diffsplice_ioe {
    tag "$event_type"
    publishDir "${params.outdir}/suppa/diffsplice_ioe", mode: 'copy'

    input:
    tuple val(event_type), path(ioe_file), path(psi_basalA), path(psi_basalB), path(tpm_basalA), path(tpm_basalB)

    output:
    tuple val(event_type), path("diffsplice_${event_type}.dpsi"), path("diffsplice_${event_type}.psivec"), emit: results
    path "${event_type}_diffsplice.log", emit: logs
    path "cluster${event_type}.clustvec", emit: clustered_events, optional: true
    path "cluster${event_type}_scores.log", emit: cluster_scores, optional: true

    script:
    """
    #!/bin/bash
    export PATH=/opt/conda/bin:\$PATH

    set -e

    echo "Verificando archivos de entrada:" > ${event_type}_diffsplice.log
    ls -l ${ioe_file} ${psi_basalA} ${psi_basalB} ${tpm_basalA} ${tpm_basalB} >> ${event_type}_diffsplice.log

    echo "Ejecutando diffSplice:" >> ${event_type}_diffsplice.log
    suppa.py diffSplice \
        -m empirical \
        -i ${ioe_file} \
        -p ${psi_basalA} ${psi_basalB} \
        -e ${tpm_basalA} ${tpm_basalB} \
        --area ${params.area} \
        --lower-bound ${params.lower_bound} \
        --alpha ${params.alpha} \
        --tpm-threshold ${params.tpm_threshold} \
        --nan-threshold ${params.nan_threshold} \
        -gc \
        -o diffsplice_${event_type} \
        >> ${event_type}_diffsplice.log 2>&1

    echo "Verificando salida de diffSplice:" >> ${event_type}_diffsplice.log
    if [ -s diffsplice_${event_type}.dpsi ] && [ -s diffsplice_${event_type}.psivec ]; then
        echo "Archivos .dpsi y .psivec generados correctamente" >> ${event_type}_diffsplice.log
        head -n 5 diffsplice_${event_type}.dpsi >> ${event_type}_diffsplice.log
        head -n 5 diffsplice_${event_type}.psivec >> ${event_type}_diffsplice.log

        event_count=\$(awk 'NR>1' diffsplice_${event_type}.dpsi | wc -l)
        echo "Número de eventos encontrados: \$event_count" >> ${event_type}_diffsplice.log

        echo "Analizando contenido de diffsplice_${event_type}.dpsi:" >> ${event_type}_diffsplice.log
        awk -F'\t' '{print NF}' diffsplice_${event_type}.dpsi | sort | uniq -c >> ${event_type}_diffsplice.log

        echo "Analizando contenido de diffsplice_${event_type}.psivec:" >> ${event_type}_diffsplice.log
        awk -F'\t' '{print NF}' diffsplice_${event_type}.psivec | sort | uniq -c >> ${event_type}_diffsplice.log

        if [ \$event_count -gt 1 ]; then
            echo "Ejecutando clusterEvents:" >> ${event_type}_diffsplice.log
            set +e  # Permitir que el script continúe si hay un error
            suppa.py clusterEvents \
                --dpsi diffsplice_${event_type}.dpsi \
                --psivec diffsplice_${event_type}.psivec \
                --sig-threshold ${params.sig_threshold} \
                --dpsi-threshold ${params.dpsi_threshold} \
                --eps ${params.eps} \
                --metric ${params.metric} \
                --min-pts ${params.min_pts} \
                --groups ${params.groups} \
                --clustering ${params.clustering_method} \
                -o cluster${event_type} \
                >> ${event_type}_diffsplice.log 2>&1
            
            if [ \$? -ne 0 ]; then
                echo "Error en clusterEvents. Analizando archivos de entrada:" >> ${event_type}_diffsplice.log
                echo "Contenido de diffsplice_${event_type}.dpsi:" >> ${event_type}_diffsplice.log
                head -n 10 diffsplice_${event_type}.dpsi >> ${event_type}_diffsplice.log
                echo "Contenido de diffsplice_${event_type}.psivec:" >> ${event_type}_diffsplice.log
                head -n 10 diffsplice_${event_type}.psivec >> ${event_type}_diffsplice.log
            fi
            set -e  # Reactivar la detención en caso de error
        else
            echo "No hay suficientes eventos para realizar el clustering" >> ${event_type}_diffsplice.log
        fi
    else
        echo "Error: Archivos de salida de diffSplice vacíos o no generados" >> ${event_type}_diffsplice.log
        exit 1
    fi

    echo "Archivos generados por SUPPA para ${event_type}:" >> ${event_type}_diffsplice.log
    ls -l >> ${event_type}_diffsplice.log
    """
}
