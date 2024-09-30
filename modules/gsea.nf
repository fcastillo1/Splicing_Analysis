process GSEA {
    tag "$meta.id"
    publishDir "${params.outdir}/gsea", mode: 'copy'

    input:
    tuple val(meta), path(rnk), path(gmt_file)

    output:
    tuple val(meta), path("${meta.id}"), emit: results
    tuple val(meta), path("${meta.id}/**/*.html"), emit: html_reports
    path "versions.yml" , emit: versions

    script:
    def rnd_seed = params.gsea_rnd_seed == 'timestamp' ? "\$(date +%s)" : params.gsea_rnd_seed
    """
    export PATH=\$PATH:/home/francisca/GSEA_4.3.3
    /home/francisca/GSEA_4.3.3/gsea-cli.sh \\
        GSEAPreranked \\
        -rnk $rnk \\
        -gmx $gmt_file \\
        -out ${meta.id} \\
        -set_min ${params.gsea_set_min} \\
        -set_max ${params.gsea_set_max} \\
        -collapse ${params.gsea_collapse} \\
        -mode ${params.gsea_mode} \\
        -create_svgs ${params.gsea_create_svgs} \\
        -include_only_symbols ${params.gsea_include_only_symbols} \\
        -make_sets ${params.gsea_make_sets} \\
        -plot_top_x ${params.gsea_plot_top_x} \\
        -rnd_seed ${rnd_seed} \\
        -zip_report ${params.gsea_zip_report}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gsea: \$(gsea-cli --version 2>&1 | sed 's/^GSEA version: //; s/ .*\$//')
    END_VERSIONS
    """
}
