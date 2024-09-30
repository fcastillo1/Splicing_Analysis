#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Importar módulos necesarios
include { download_srr } from './modules/download_srr'
include { fastqc_raw } from './modules/fastqc_raw'
include { trimmomatic } from './modules/trimmomatic'
include { fastqc_trimmed } from './modules/fastqc_trimmed'
include { star } from './modules/star'
include { generate_bigwig } from './modules/generate_bigwig'
include { multiqc_run } from './modules/multiqc'
include { stringtie } from './modules/stringtie'
include { prep_de } from './modules/prep_de'
include { stringtie_merge } from './modules/stringtie_merge'
include { RMATS_ORIGINAL_GTF; RMATS_MERGED_GTF } from './modules/rmats'
include { RMATS_SASHIMIPLOT } from './modules/rmats2sashimiplot'
include { SALMON_INDEX } from './modules/salmon_index'
include { SALMON_QUANT } from './modules/salmon_quant'
include { GFFREAD_TX2GENE } from './modules/gffread_tx2gene'
include { SAMTOOLS_SORT_INDEX } from './modules/samtools'
include { TXIMPORT } from './modules/tximport'
include { suppa2_generate_events } from './modules/suppa2_generate_events'
include { suppa2_generate_isoforms } from './modules/suppa2_generate_isoforms'
include { suppa2_psi_per_event } from './modules/suppa2_psi_per_event'
include { suppa_per_isoform } from './modules/suppa_per_isoform'
include { split_psi_ioe } from './modules/split_psi_ioe'
include { split_psi_ioi } from './modules/split_psi_ioi'
include { split_tpm_file } from './modules/split_tpm_file'
include { suppa_diffsplice_ioe } from './modules/suppa_diffsplice_ioe'
include { suppa_diffsplice_ioi } from './modules/suppa_diffsplice_ioi'
include { DESEQ2_DIFFERENTIAL } from './modules/deseq2_differential'
include { GENERATE_INTERACTIVE_REPORT } from './modules/generate_interactive_report'
include { GSEA } from './modules/gsea'
include { GPROFILER2_GOST } from './modules/gprofiler'
include { GENERATE_BACKGROUND } from './modules/generate_background'

// Funciones auxiliares
def needsDownload(fastq, fastq2) {
    return (!fastq && !fastq2) || (fastq == '' && fastq2 == '')
}

def filterEventsByGenes(file, geneList) {
    def filteredEvents = [:]
    def header = file.readLines()[0]
    
    if (geneList.isEmpty()) {
        filteredEvents['all'] = file
    } else {
        geneList.each { gene ->
            def filteredLines = file.readLines().findAll { line ->
                line == header || line.split("\t")[1].trim() == gene
            }
            if (filteredLines.size() > 1) {
                def filteredFile = file("${params.outdir}/filtered_${gene}_${file.getName()}")
                filteredFile.text = filteredLines.join("\n")
                filteredEvents[gene] = filteredFile
            }
        }
    }
    
    return filteredEvents
}

def readGeneList(file) {
    return file.readLines().collect { it.trim() }.findAll { it }
}

def getLibType(tool, strandedness) {
    def libTypeMap = [
        'first-strand': ['rmats': 'fr-firststrand', 'salmon': 'ISR'],
        'second-strand': ['rmats': 'fr-secondstrand', 'salmon': 'ISF'],
        'unstranded': ['rmats': 'fr-unstranded', 'salmon': 'IU'],
        'unknown': ['rmats': 'fr-unstranded', 'salmon': 'A']
    ]
    return libTypeMap.get(strandedness, 'unknown')[tool]
}

workflow {
    log.info "Iniciando pipeline de análisis de splicing diferencial"

    // Inicializar multiqc_files como un canal vacío
    multiqc_files = Channel.empty()

    // Verificar parámetros obligatorios y existencia de archivos
    ['samplesheet', 'outdir', 'gtf', 'genome_fasta', 'transcript_fasta'].each { param ->
        if (!params[param]) error "El parámetro --${param} es obligatorio"
        if (!file(params[param]).exists()) error "El archivo/directorio ${params[param]} no existe"
    }

    // Crear canal de entrada desde el samplesheet
    input_samples = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            if (!row.sample_id || !row.condition) error "Las columnas 'sample_id' y 'condition' son obligatorias en el samplesheet"
            def needs_download = needsDownload(row.fastq, row.fastq2)
            tuple(row.sample_id, row.fastq, row.fastq2, row.condition, needs_download)
        }

    // Separar muestras que necesitan descarga de las que ya tienen archivos locales
    input_samples
        .branch {
            download: it[4]
            local: !it[4]
        }
        .set { branched_input }

    // Descargar SRR si es necesario
    download_ch = download_srr(branched_input.download.map { row -> tuple(row[0], row[3]) })

    // Preparar canal para muestras locales
    local_ch = branched_input.local.map { row ->
        def reads = params.single_end ? [file(row[1])] : [file(row[1]), file(row[2])]
        if (reads.any { !it.exists() }) {
            error "Archivo(s) de lecturas no encontrado(s) para la muestra ${row[0]}"
        }
        tuple(row[0], reads, row[3])
    }

    // Combinar muestras descargadas y locales
    read_pairs_ch = local_ch.mix(download_ch)

    // FastQC en lecturas crudas
    fastqc_raw_ch = fastqc_raw(read_pairs_ch)
    multiqc_files = multiqc_files.mix(fastqc_raw_ch.fastqc_output.map { it[1] }.flatten())

    // Trimmomatic
    trimmed_reads = trimmomatic(read_pairs_ch)

    // FastQC en lecturas recortadas
    fastqc_trimmed_ch = fastqc_trimmed(trimmed_reads.trimmed_reads)
    multiqc_files = multiqc_files.mix(fastqc_trimmed_ch.fastqc_output.map { it[1] }.flatten())
    
    // STAR alignment (si está habilitado)
    if (params.aligner_star || params.full_pipeline) {
        // Filtrar solo los archivos paired para STAR
        paired_reads = trimmed_reads.trimmed_reads.map { sample_id, reads, condition ->
            def paired = reads.findAll { it.name.contains("paired") && !it.name.contains("unpaired") }
            if (paired.size() != 2) {
                error "Se esperaban exactamente 2 archivos paired para la muestra $sample_id, pero se encontraron ${paired.size()}: ${paired*.name}"
            }
            tuple(sample_id, paired, condition)
        }

        star_results = star(paired_reads, file(params.genomeDir), file(params.gtf))

        // SAMTOOLS_SORT_INDEX
        SAMTOOLS_SORT_INDEX(star_results.bam)
        SAMTOOLS_SORT_INDEX.out.sorted_bam_and_index
            .ifEmpty { error "El proceso SAMTOOLS_SORT_INDEX no produjo resultados" }

        // Generar bigwig
        generate_bigwig(SAMTOOLS_SORT_INDEX.out.sorted_bam_and_index)

        // Ejecutar StringTie
        stringtie_results = stringtie(SAMTOOLS_SORT_INDEX.out.sorted_bam_and_index, file(params.gtf))

        // Recolectar GTFs para merge
        stringtie_gtfs = stringtie_results.stringtie_gtf
            .map { it[1] }
            .collect()
        
        stringtie_gtfs.ifEmpty { error "No se generaron archivos GTF en el proceso StringTie" }

        // Preparar matrices de conteo con prepDE.py
        stringtie_dge_gtfs = stringtie_results.stringtie_dge_gtf
            .map { it[1] }
            .collect()
        prep_de_results = prep_de(stringtie_dge_gtfs)

        // StringTie merge
        stringtie_merge_results = stringtie_merge(stringtie_gtfs, file(params.gtf))

        if (params.run_rmats) {
            log.info "Iniciando análisis rMATS"

            // Verificar que los archivos BAM existen
            bam1_file = file(params.bam1)
            bam2_file = file(params.bam2)
            if (!bam1_file.exists()) {
                error "El archivo BAM1 no existe: ${params.bam1}"
            }
            if (!bam2_file.exists()) {
                error "El archivo BAM2 no existe: ${params.bam2}"
            }

            // Verificar que el GTF original existe
            original_gtf = file(params.gtf)
            if (!original_gtf.exists()) {
                error "El archivo GTF original no existe: ${params.gtf}"
            }

            // Ejecutar rMATS con GTF original
            rmats_original_results = RMATS_ORIGINAL_GTF(original_gtf, bam1_file, bam2_file)

            // Verificar si el GTF fusionado existe
            if (stringtie_merge_results.merged_gtf) {
                merged_gtf = stringtie_merge_results.merged_gtf

                // Ejecutar rMATS con GTF fusionado
                rmats_merged_results = RMATS_MERGED_GTF(merged_gtf, bam1_file, bam2_file)

                // Asumir que necesitas el archivo con los resultados del GTF fusionado
                rmats_final_results = rmats_merged_results.results
            } else {
                log.warn "GTF fusionado no encontrado. Ejecutando rMATS solo con GTF original."
                rmats_final_results = rmats_original_results.results
            }

            // Ejecutar RMATS_SASHIMIPLOT
            def event_type = params.sashimi_event_type
            rmats_results = Channel
                .fromPath("${params.rmats_output_dir}/${event_type}.MATS.JC.txt")
                .map { file -> tuple(event_type, file) }

            RMATS_SASHIMIPLOT(
                rmats_results,
                bam1_file,
                bam2_file,
                file(params.sashimi_group_file),
                params.genes_of_interest,
                params.label1,
                params.label2
            )

            // Recoger y publicar los logs
            RMATS_SASHIMIPLOT.out.logs.collectFile(name: 'sashimi_plot_generation.log', storeDir: params.outdir)
        }
    }
    
    if (params.aligner_salmon || params.full_pipeline) {
        SALMON_INDEX(
            file(params.genome_fasta),
            file(params.transcript_fasta)
        )

        // Preparar el canal de entrada para Salmon
        trimmed_reads_for_salmon = trimmed_reads.trimmed_reads.map { sample_id, reads, condition ->
            def meta = [id: sample_id, single_end: false]
            tuple(meta, reads)
        }

        // Salmon Quantification
        SALMON_QUANT(
            trimmed_reads_for_salmon,
            SALMON_INDEX.out.index,
            file(params.gtf),
            file(params.transcript_fasta)
        )

        GFFREAD_TX2GENE(file(params.gtf))

        TXIMPORT(
            SALMON_QUANT.out.results.map { it[1] }.collect(),
            GFFREAD_TX2GENE.out.tx2gene
        )

        // Create a channel from tximport results for DESeq2
        deseq2_input = TXIMPORT.out.counts_gene
            .map { file ->
                if (file.exists()) {
                    return tuple([ id: 'deseq2' ], file)
                } else {
                    error "The Salmon-tximport count file was not generated correctly: $file"
                }
            }

        // Verificar si el canal está vacío
        deseq2_input.ifEmpty { error "tximport no generó el archivo de conteos de genes" }

        DESEQ2_DIFFERENTIAL(
            deseq2_input,
            file(params.samplesheet),
            "${params.contrast_variable},${params.reference_level},${params.target_level}",
            file(params.gtf)
        )

        DESEQ2_DIFFERENTIAL.out.results
            .ifEmpty { error "No se generaron resultados de DESeq2" }

        // Agregar resultados a multiqc_files
        multiqc_files = multiqc_files.mix(DESEQ2_DIFFERENTIAL.out.results.map { it[1] })
        multiqc_files = multiqc_files.mix(DESEQ2_DIFFERENTIAL.out.ma_plot.map { it[1] })
        multiqc_files = multiqc_files.mix(DESEQ2_DIFFERENTIAL.out.volcano_plot_labeled.map { it[1] })
        multiqc_files = multiqc_files.mix(DESEQ2_DIFFERENTIAL.out.volcano_plot_no_labels.map { it[1] })
        multiqc_files = multiqc_files.mix(DESEQ2_DIFFERENTIAL.out.pca_plot.map { it[1] })
        multiqc_files = multiqc_files.mix(DESEQ2_DIFFERENTIAL.out.heatmap.map { it[1] })
        multiqc_files = multiqc_files.mix(DESEQ2_DIFFERENTIAL.out.dispersion_plot.map { it[1] })
       
        // Agregar resultados de Salmon a multiqc_files
        multiqc_files = multiqc_files.mix(SALMON_QUANT.out.results.collect{it[1]})
        
        // Preparar input para GSEA
        deseq2_results_for_gsea = DESEQ2_DIFFERENTIAL.out.results
            .map { meta, results_file ->
                def gsea_dir = file("${params.outdir}/gsea")
                if (!gsea_dir.exists()) {
                    gsea_dir.mkdirs()
                }
                
                def rnk_file = gsea_dir.resolve("${meta.id}_for_gsea.rnk")
                
                def lines = results_file.readLines()
                if (lines.size() > 1) {
                    def header = lines[0].split(',')
                    
                    def geneNameIndex = header.findIndexOf { it.replaceAll('"', '').toLowerCase() == "gene_name" }
                    def log2FCIndex = header.findIndexOf { it.replaceAll('"', '').toLowerCase() == "log2foldchange" }
                    
                    if (geneNameIndex == -1 || log2FCIndex == -1) {
                        error "No se encontraron las columnas necesarias en el archivo de resultados de DESeq2"
                    }
                    
                    def rnkContent = lines.subList(1, lines.size()).collect { line ->
                        def fields = line.split(',').collect { it.replaceAll('"', '') }
                        if (fields.size() > geneNameIndex && fields.size() > log2FCIndex) {
                            "${fields[geneNameIndex]}\t${fields[log2FCIndex]}"
                        } else {
                            null
                        }
                    }.findAll { it != null }
                    
                    rnk_file.text = rnkContent.sort { -it.split('\t')[1].toDouble() }.join('\n')
                } else {
                    log.warn "El archivo de resultados de DESeq2 está vacío o solo contiene encabezado"
                    rnk_file.text = ""
                }
                
                [ meta, rnk_file ]
            }
        
        // Crear un canal para los archivos GMT
        gmt_files = Channel.fromPath(params.gsea_gmt_symbols)
        
        // Combinar los resultados de DESeq2 con cada archivo GMT
        gsea_input = deseq2_results_for_gsea
            .combine(gmt_files)
            .map { meta, rnk, gmt -> [ meta, rnk, gmt ] }
        
        // Ejecutar GSEA para cada combinación de resultados y archivo GMT
        GSEA(gsea_input)

        GENERATE_BACKGROUND(DESEQ2_DIFFERENTIAL.out.results)

        GPROFILER2_GOST(
            DESEQ2_DIFFERENTIAL.out.results,
            GENERATE_BACKGROUND.out.background
        )
                
        GENERATE_INTERACTIVE_REPORT(
            DESEQ2_DIFFERENTIAL.out.results,
            DESEQ2_DIFFERENTIAL.out.normalized_counts,
            file(params.gtf),
            file(params.samplesheet),
            GPROFILER2_GOST.out.all_enrich,
            GPROFILER2_GOST.out.plot_png,
            file("$projectDir/modules/DESeq2_Interactive_Report.Rmd")
        )

        // Agregar resultados de gProfiler a multiqc_files
        multiqc_files = multiqc_files.mix(GPROFILER2_GOST.out.all_enrich.map { it[1] })
    }

    if (params.run_suppa) {
        // Generar eventos (IOE)
        suppa2_generate_events_ch = suppa2_generate_events(file(params.gtf))

        // Calcular PSI por evento (IOE)
        suppa2_psi_per_event_ch = suppa2_psi_per_event(
            suppa2_generate_events_ch.events.flatten(),
            TXIMPORT.out.suppa_tpm
        )

        // Splitting de TPM
        split_tpm_ch = split_tpm_file(TXIMPORT.out.suppa_tpm, file(params.samplesheet))

        // Splitting de PSI para IOE
        split_psi_ioe_ch = split_psi_ioe(suppa2_psi_per_event_ch.psi.flatten(), file(params.samplesheet))

        // Preparar input para diffsplice IOE
        def event_types = ['A3', 'A5', 'AF', 'AL', 'MX', 'RI', 'SE']
        diffsplice_input_ioe_ch = Channel.fromList(event_types)
            .combine(split_psi_ioe_ch.collect())
            .combine(split_tpm_ch.split_tpm.collect())
            .map { it ->
                def (event_type, psi_files, tpm_files) = [it[0], it[1..-3], it[-2..-1]]
                def ioe_file = file("${params.outdir}/suppa/generate_events/events_${event_type}_strict.ioe")
                def psi_basalA = psi_files.find { it.name.contains("${event_type}_strict_basalA") }
                def psi_basalB = psi_files.find { it.name.contains("${event_type}_strict_basalB") }
                def tpm_basalA = tpm_files.find { it.name.contains('basalA') }
                def tpm_basalB = tpm_files.find { it.name.contains('basalB') }
                
                if (!ioe_file.exists()) {
                    log.warn "IOE file not found for event type ${event_type}: ${ioe_file}"
                    return null
                }
                if (!psi_basalA || !psi_basalB) {
                    log.warn "PSI files not found for event type ${event_type}: basalA=${psi_basalA}, basalB=${psi_basalB}"
                    return null
                }
                if (!tpm_basalA || !tpm_basalB) {
                    log.warn "TPM files not found for event type ${event_type}: basalA=${tpm_basalA}, basalB=${tpm_basalB}"
                    return null
                }
                
                tuple(event_type, ioe_file, psi_basalA, psi_basalB, tpm_basalA, tpm_basalB)
            }
            .filter { it != null }

        // Ejecutar diffsplice para IOE
        suppa_diffsplice_ioe(diffsplice_input_ioe_ch)

        // Generar isoformas (IOI)
        suppa2_generate_isoforms_ch = suppa2_generate_isoforms(file(params.gtf))

        // Calcular PSI por isoforma (IOI)
        suppa_per_isoform_ch = suppa_per_isoform(
            file(params.gtf),
            TXIMPORT.out.suppa_tpm
        )

        // Splitting de PSI para IOI
        split_psi_ioi_ch = split_psi_ioi(suppa_per_isoform_ch.psi_isoform, file(params.samplesheet))
        
        // Preparar input para diffsplice IOI
        diffsplice_input_ioi_ch = suppa2_generate_isoforms_ch.isoforms
            .combine(split_psi_ioi_ch.collect())
            .combine(split_tpm_ch.split_tpm.collect())
            .map { it ->
                def (ioi_file, psi_files, tpm_files) = [it[0], it[1..-3], it[-2..-1]]
                def psi_basalA = psi_files.find { it.name.contains("basalA") }
                def psi_basalB = psi_files.find { it.name.contains("basalB") }
                def tpm_basalA = tpm_files.find { it.name.contains('basalA') }
                def tpm_basalB = tpm_files.find { it.name.contains('basalB') }
                
                if (!ioi_file.exists()) {
                    log.warn "IOI file not found: ${ioi_file}"
                    return null
                }
                if (!psi_basalA || !psi_basalB) {
                    log.warn "PSI IOI files not found: basalA=${psi_basalA}, basalB=${psi_basalB}"
                    return null
                }
                if (!tpm_basalA || !tpm_basalB) {
                    log.warn "TPM files not found: basalA=${tpm_basalA}, basalB=${tpm_basalB}"
                    return null
                }
                
                tuple(ioi_file, psi_basalA, psi_basalB, tpm_basalA, tpm_basalB)
            }
            .filter { it != null }

        // Ejecutar diffsplice para IOI
        suppa_diffsplice_ioi(diffsplice_input_ioi_ch)
    }

    // Ejecutar MultiQC
    multiqc_run(multiqc_files.collect())
}

// Mensaje de finalización
workflow.onComplete {
    log.info "Pipeline completado en: $workflow.complete"
    log.info "Estado de ejecución: ${workflow.success ? 'OK' : 'fallido'}"
    log.info "Tiempo de ejecución: $workflow.duration"
    log.info "Directorio de salida: ${params.outdir}"
}

// Manejo de errores
workflow.onError {
    log.error "Pipeline fallido: $workflow.errorMessage"
}
