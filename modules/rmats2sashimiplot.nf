process RMATS_SASHIMIPLOT {
    tag "Generating Sashimi plots for ${event_type} in ${genes_of_interest}"
    publishDir "${params.outdir}/sashimi_plots", mode: 'copy'

    input:
    tuple val(event_type), path(rmats_output)
    path bam1
    path bam2
    path grouping_file
    val genes_of_interest
    val label1
    val label2

    output:
    path "*_sashimi_plot", emit: sashimi_plots, optional: true
    path "*.log", emit: logs

    script:
    """
    echo "Processing event type: ${event_type}" > process.log
    echo "Input file: ${rmats_output}" >> process.log
    echo "Genes of interest: ${genes_of_interest}" >> process.log

    # Filtrar por FDR (columna 20) y IncLevelDifference (columna 23)
    awk -F'\\t' '
        NR==1 { print \$0 > "filtered_significant.txt" }
        NR>1 && \$20 <= 0.05 && (\$23 >= 0.2 || \$23 <= -0.2) {
            print \$0 >> "filtered_significant.txt"
        }
    ' ${rmats_output}

    for gene in ${genes_of_interest.split(',').join(' ')}; do
        echo "Processing gene: \$gene" >> process.log
        
        # Extraer el encabezado y las líneas correspondientes al gen del archivo filtrado
        (head -n 1 filtered_significant.txt && grep "\$gene" filtered_significant.txt) > "\${gene}_${event_type}_events_with_header.txt"
        
        # Generar Sashimi plot si se encontraron eventos significativos
        if [ -s "\${gene}_${event_type}_events_with_header.txt" ] && [ \$(wc -l < "\${gene}_${event_type}_events_with_header.txt") -gt 1 ]; then
            echo "Significant events found for \$gene in ${event_type}, generating sashimi plot" >> process.log
            rmats2sashimiplot --b1 $bam1 \\
                --b2 $bam2 \\
                --event-type $event_type \\
                -e "\${gene}_${event_type}_events_with_header.txt" \\
                --exon_s 1 \\
                --intron_s 5 \\
                --l1 $label1 \\
                --l2 $label2 \\
                -o "\${gene}_${event_type}_sashimi_plot" \\
                --group-info $grouping_file \\
                --font-size 8 \\
                --fig-height 8 \\
                --fig-width 12 \\
                2> >(grep -v "Skipping read with CIGAR" >&2)
        else
            echo "Warning: No significant events found for gene \$gene in ${event_type} (FDR ≤ 0.05 and |IncLevelDifference| ≥ 0.2)" >> process.log
        fi
    done

    # Verificar si se generaron sashimi plots
    if ls *sashimiplot/*.pdf 1> /dev/null 2>&1; then
        echo "Sashimi plots generated successfully" >> process.log
    else
        echo "No sashimi plots were generated" >> process.log
    fi
    """
}
