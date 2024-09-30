process DESEQ2_DIFFERENTIAL {
    tag "$meta.id"
    publishDir "${params.outdir}/deseq2", mode: 'copy'

    input:
    tuple val(meta), path(counts)
    path samplesheet
    val(contrast_variable)
    val(reference_level)
    val(target_level)
    path gtf_file

    output:
    tuple val(meta), path("*_results.csv"), emit: results
    tuple val(meta), path("*_ma_plot.pdf"), emit: ma_plot
    tuple val(meta), path("*_pca_plot.pdf"), emit: pca_plot
    tuple val(meta), path("*_heatmap.pdf"), emit: heatmap
    tuple val(meta), path("*_dispersion_plot.pdf"), emit: dispersion_plot
    tuple val(meta), path("*_volcano_plot_labeled.pdf"), emit: volcano_plot_labeled
    tuple val(meta), path("*_volcano_plot_no_labels.pdf"), emit: volcano_plot_no_labels
    tuple val(meta), path("*_dds.rds"), emit: rds
    tuple val(meta), path("*_sizefactors.csv"), emit: size_factors
    tuple val(meta), path("*_normalized_counts.csv"), emit: normalized_counts
    tuple val(meta), path("*_R_sessionInfo.log"), emit: log
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    Rscript $projectDir/modules/deseq_differential.R \\
        --count_file '${counts}' \\
        --sample_file '${samplesheet}' \\
        --contrast '${contrast_variable},${reference_level},${target_level}' \\
        --output_prefix '${prefix}' \\
        --gtf_file '${gtf_file}' \\
        --lfc_threshold ${params.lfc_threshold ?: 1} \\
        --alpha ${params.alpha ?: 0.05} \\
        --p_adjust_method '${params.p_adjust_method ?: 'BH'}' \\
        --shrink_lfc ${params.shrink_lfc ?: 'true'} \\
        --cores ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | grep "R version" | sed 's/R version //; s/ .*//')
        bioconductor-deseq2: \$(Rscript -e "library(DESeq2); cat(as.character(packageVersion('DESeq2')))")
    END_VERSIONS
    """
}
