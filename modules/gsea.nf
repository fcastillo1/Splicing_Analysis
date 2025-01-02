process GSEA {
    tag "$meta.id"
    publishDir "${params.outdir}/gsea", mode: 'copy'

    input:
    tuple val(meta), path(rnk)
    path gmt_files

    output:
    tuple val(meta), path("${meta.id}"), emit: results
    tuple val(meta), path("${meta.id}/**/*.html"), emit: html_reports 
    path "versions.yml", emit: versions

    when:
    params.gsea_run

    script:
    """
    mkdir -p ${meta.id}
    
    for gmt_file in $gmt_files; do
        db_name=\$(basename \$gmt_file .gmt)
        mkdir -p ${meta.id}/\$db_name
        
         # Ejecucion del comando de gsea
        gsea-cli.sh GSEAPreranked \
            -rnk $rnk \
            -gmx \$gmt_file \
            -out ${meta.id}/\$db_name \
            -nperm ${params.gsea_nperm} \
            -permute ${params.gsea_permute} \
            -scoring_scheme ${params.gsea_scoring_scheme} \
            -metric ${params.gsea_metric} \
            -sort ${params.gsea_sort} \
            -order ${params.gsea_order} \
            -set_max ${params.gsea_set_max} \
            -set_min ${params.gsea_set_min} \
            -norm ${params.gsea_norm} \
            -rnd_type ${params.gsea_rnd_type} \
            -make_sets ${params.gsea_make_sets} \
            -median ${params.gsea_median} \
            -num ${params.gsea_num} \
            -plot_top_x ${params.gsea_plot_top_x} \
            -rnd_seed ${params.gsea_rnd_seed} \
            -save_rnd_lists ${params.gsea_save_rnd_lists} \
            -zip_report ${params.gsea_zip_report}
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gsea: \$(gsea-cli.sh --version 2>&1 | sed 's/^GSEA version: //; s/ .*\$//')
    END_VERSIONS
    """
}
