process RMATS_SASHIMIPLOT {
    publishDir "${params.outdir}/sashimiplot_results", mode: 'copy'

    input:
    path rmats_results_dir
    path bam1
    path bam2
    path merged_gtf  // Añadimos el GTF del merge como input

    output:
    path "sashimi_plots/**/*"
    path "sashimi_plot.log" // Añadimos esto para poder inspeccionar los archivos de entrada si es necesario

    script:
    """
    mkdir -p rmats_results
    
    # Mover archivos de rMATS, ignorando errores si algún archivo o directorio no existe
    for file in ${rmats_results_dir}/*; do
        if [ -e "\$file" ]; then
            mv "\$file" rmats_results/ 2>/dev/null || true
        fi
    done

    # Copiar el GTF del merge a rmats_results (no lo movemos por si se necesita en otros procesos)
    cp ${merged_gtf} rmats_results/merged.gtf

    Rscript ${projectDir}/modules/generate_sashimiplots.R \
        --rMATSdir rmats_results \
        --eventTypes ${params.sashimi_event_type ?: 'SE,A3SS,A5SS,MXE,RI'} \
        --bam1 ${bam1} \
        --bam2 ${bam2} \
        --outdir sashimi_plots \
        --label1 '${params.label1}' \
        --label2 '${params.label2}'
    """
}
