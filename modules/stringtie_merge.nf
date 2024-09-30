process stringtie_merge {
    label 'mid_memory'
    publishDir "${params.outdir}/stringtie_merge", mode: 'copy'

    input:
    path gtfs
    path reference_gtf

    output:
    path "gffcmp.annotated.corrected.gtf", emit: merged_gtf
    path "gffcmp.*", emit: gffcmp

    script:
    """
    ls -1 *.gtf > assembly_gtf_list.txt
    stringtie --merge -G $reference_gtf -o stringtie_merged.gtf assembly_gtf_list.txt -p $task.cpus
    gffcompare -R -V -r $reference_gtf stringtie_merged.gtf
    Rscript $projectDir/modules/correct_gene_names.R
    gffread -E gffcmp.annotated.corrected.gff -T -o gffcmp.annotated.corrected.gtf
    """
}
