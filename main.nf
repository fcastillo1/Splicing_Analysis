#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// =======================================
// IMPORTACION DE MODULOS
// =======================================

// Módulos de control de calidad y preprocesamiento
include { download_srr } from './modules/download_srr'          
include { fastqc_raw } from './modules/fastqc_raw'             
include { trimmomatic } from './modules/trimmomatic'           
include { fastqc_trimmed } from './modules/fastqc_trimmed'     

// Módulos de alineamiento y cuantificación
include { star } from './modules/star'                         
include { generate_bigwig } from './modules/generate_bigwig'   
include { SALMON_INDEX } from './modules/salmon_index'         
include { SALMON_QUANT } from './modules/salmon_quant'         
include { SAMTOOLS_SORT_INDEX } from './modules/samtools'      

// Módulos de splicing alternativo
include { RMATS_ORIGINAL_GTF } from './modules/rmats'          
include { RMATS_SASHIMIPLOT } from './modules/rmats2sashimiplot'
include { GFFREAD_TX2GENE } from './modules/gffread_tx2gene'   
include { TXIMPORT } from './modules/tximport'                 

// Módulos SUPPA
include { suppa2_generate_events } from './modules/suppa2_generate_events'
include { suppa2_generate_isoforms } from './modules/suppa2_generate_isoforms'
include { suppa2_psi_per_event } from './modules/suppa2_psi_per_event'
include { suppa_per_isoform } from './modules/suppa_per_isoform'
include { split_psi_ioe } from './modules/split_psi_ioe'
include { split_psi_ioi } from './modules/split_psi_ioi'
include { split_tpm_file } from './modules/split_tpm_file'
include { suppa_diffsplice_ioe } from './modules/suppa_diffsplice_ioe'
include { suppa_diffsplice_ioi } from './modules/suppa_diffsplice_ioi'

// Módulos de análisis diferencial y enriquecimiento
include { DESEQ2_DIFFERENTIAL } from './modules/deseq2_differential'
include { GSEA } from './modules/gsea'
include { GPROFILER } from './modules/gprofiler'
include { GENERATE_BACKGROUND } from './modules/generate_background'
include { GENERATE_DESEQ2_REPORT } from './modules/generate_deseq2_report'
include { PREPARE_GSEA_INPUT } from './modules/prepare_gsea_input'

// Módulo de control de calidad global
include { multiqc_run } from './modules/multiqc'              

// Funciones auxiliares
def needsDownload(fastq, fastq2, trimmed_fastq, trimmed_fastq2) {
    if (params.use_trimmed) {
        return (!trimmed_fastq && !trimmed_fastq2) || (trimmed_fastq == '' && trimmed_fastq2 == '')
    }
    return (!fastq && !fastq2) || (fastq == '' && fastq2 == '')
}

// Funcion que filtra los eventos por genes
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

// Funcion que lee la lista de geens
def readGeneList(file) {
    return file.readLines().collect { it.trim() }.findAll { it }
}

// Funcion que permite obtener el tipo de libreria
def getLibType(tool, strandedness) {
    def libTypeMap = [
        'first-strand': ['rmats': 'fr-firststrand', 'salmon': 'ISR'],
        'second-strand': ['rmats': 'fr-secondstrand', 'salmon': 'ISF'],
        'unstranded': ['rmats': 'fr-unstranded', 'salmon': 'IU'],
        'unknown': ['rmats': 'fr-unstranded', 'salmon': 'A']
    ]
    return libTypeMap.get(strandedness, 'unknown')[tool]
}

// =======================================
// INICIO DE WORKFLOW DE TRABAJO
// =======================================
workflow {
    log.info "Iniciando pipeline de análisis de splicing diferencial y expresión diferencial"

// Inicializar multiqc_files como un canal vacío
    multiqc_files = Channel.empty()

    // Condiciones de parametros que tienen que ser obligatorios para el funcionamiento del pipelien
    ['samplesheet', 'outdir', 'gtf', 'genome_fasta', 'transcript_fasta'].each { param ->
        if (!params[param]) error "El parámetro --${param} es obligatorio"
        if (!file(params[param]).exists()) error "El archivo/directorio ${params[param]} no existe"
    }

    // Funcion que lee los inputs
    input_samples = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            if (!row.sample_id || !row.condition) error "Las columnas 'sample_id' y 'condition' son obligatorias"
            
            def needs_download = needsDownload(row.fastq, row.fastq2, row.trimmed_fastq, row.trimmed_fastq2)
            
            if (params.use_trimmed && row.trimmed_fastq && row.trimmed_fastq2) {
                def reads = [file(row.trimmed_fastq), file(row.trimmed_fastq2)]
                // Importante: cambiamos el orden de los elementos en el tuple
                tuple(row.sample_id, reads, row.condition, false)  // false porque tenemos los archivos trimados
            } else if (!params.use_trimmed && (row.fastq || row.fastq2)) {
                // Mantenemos la estructura para archivos crudos
                tuple(row.sample_id, row.fastq, row.fastq2, row.condition, needs_download)
            } else {
                error "No se encontraron archivos válidos para la muestra ${row.sample_id}"
            }
        }

    input_samples
        .branch {
            download: it[4]
            local: !it[4]
        }
        .set { branched_input }

    // Funcion que incia la descarga de ser necesario
    download_ch = download_srr(
        branched_input.download.map { row -> tuple(row[0], row[3]) }
    )

    // Funcion que recoge los datos locales
    local_ch = branched_input.local.map { row ->
        try {
            // Si estamos usando archivos trimados
            if (params.use_trimmed) {
                // Los archivos ya vienen como una lista en reads
                tuple(row[0], row[1], row[2])  // sample_id, reads, condition
            } else {
                // Si son archivos crudos
                if (params.single_end) {
                    def read = file(row[1])
                    if (!read.exists()) {
                        error "Archivo no encontrado: ${read}"
                    }
                    tuple(row[0], [read], row[3])
                } else {
                    def read1 = file(row[1])
                    def read2 = file(row[2])
                    if (!read1.exists() || !read2.exists()) {
                        error "Archivos no encontrados: ${read1} y/o ${read2}"
                    }
                    tuple(row[0], [read1, read2], row[3])
                }
            }
        } catch (Exception e) {
            error "Error procesando archivos para muestra ${row[0]}: ${e.message}"
        }
    }

    read_pairs_ch = local_ch.mix(download_ch)

    // Cuando se activa la condicion de usar datos trimmeados se ejecuta de esta manera el pipeline
    if (!params.use_trimmed) {
        // Flujo normal con trimming
        fastqc_raw_ch = fastqc_raw(read_pairs_ch)
        // MultiQC
        multiqc_files = multiqc_files.mix(fastqc_raw_ch.fastqc_output.map { it[1] }.flatten())
        // Trimmed
        trimmed_reads = trimmomatic(read_pairs_ch)
        
        processed_reads = trimmed_reads.trimmed_reads.map { sample_id, reads, condition ->
            def processed = params.single_end ? reads[0] : reads
            tuple(sample_id, processed, condition)
        }

        fastqc_trimmed_ch = fastqc_trimmed(processed_reads)
        multiqc_files = multiqc_files.mix(fastqc_trimmed_ch.fastqc_output.map { it[1] }.flatten())
    } else {
        processed_reads = read_pairs_ch
    }

// =================================================================
// PROCESAMIENTO DE LOS DATOS CON STAR (RAMA IZQUIERDA DEL PIPELINE)
// =================================================================
    if (params.aligner_star || params.full_pipeline) {
        paired_reads = processed_reads.map { sample_id, reads, condition ->
            if (!params.use_trimmed) {
                def paired = reads.findAll { it.name.contains("paired") && !it.name.contains("unpaired") }
                if (paired.size() != 2 && !params.single_end) {
                    error "Se esperaban 2 archivos paired para $sample_id"
                }
                tuple(sample_id, paired, condition)
            } else {
                tuple(sample_id, reads, condition)
            }
        }
        // Ejecucion de STAR
        star_results = star(paired_reads, file(params.genomeDir), file(params.gtf))
        // Conversion con Samtools
        bam_results = SAMTOOLS_SORT_INDEX(star_results.aligned_reads)
        // Generacion de bigwig
        generate_bigwig(bam_results.sorted_bam_and_index)

        // Si se activa rMATS se inicializa el proceso
        if (params.run_rmats) {
            log.info "Iniciando análisis rMATS"
            // Se cargan los bam y la anotacion
            bam1_file = file(params.bam1)
            bam2_file = file(params.bam2)
            original_gtf = file(params.gtf)

            // Se generan las condiciones cuando no existan los archivos
            if (!bam1_file.exists()) error "El archivo BAM1 no existe: ${params.bam1}"
            if (!bam2_file.exists()) error "El archivo BAM2 no existe: ${params.bam2}"
            if (!original_gtf.exists()) error "El archivo GTF original no existe: ${params.gtf}"

            rmats_trigger = bam_results.sorted_bam_and_index.collect()
            // Se ejecuta rMATS
            RMATS_ORIGINAL_GTF(original_gtf, bam1_file, bam2_file)
                .subscribe { log.info "rMATS completado" }
        }
    }

// =================================================================
// PROCESAMIENTO DE LOS DATOS CON SALMON (RAMA DERECHA DEL PIPELINE)
// =================================================================
    if (params.aligner_salmon || params.full_pipeline) {
        Channel.fromPath(params.gtf).set { gtf_ch }
        Channel.fromPath(params.samplesheet).set { samplesheet_ch }

        // Se realiza el indice de salmon
        SALMON_INDEX(file(params.transcript_fasta), file(params.genome_fasta))

        salmon_reads = processed_reads.map { sample_id, reads, condition ->
            if (params.single_end) {
                def read = params.use_trimmed ? reads[0] : reads.find { it.name.contains("trimmed.fastq.gz") }
                tuple(sample_id, read, condition)
            } else {
                def paired_reads = params.use_trimmed ? reads : reads.findAll {
                    it.name.contains("paired_trimmed.fastq.gz") && !it.name.contains("unpaired")
                }
                tuple(sample_id, paired_reads, condition)
            }
        }

        // Se genera la cuantificacion
        SALMON_QUANT(salmon_reads, SALMON_INDEX.out.index)
        GFFREAD_TX2GENE(gtf_ch)
        
        // Se realiza tximport 
        TXIMPORT(
            SALMON_QUANT.out.results.map { it[1] }.collect(),
            GFFREAD_TX2GENE.out.tx2gene,
            file("$projectDir/tximport.R")
        )

        // Se inicia el analisis de expresion diferencial 
        DESEQ2_DIFFERENTIAL(
            TXIMPORT.out.counts_gene_length_scaled,
            samplesheet_ch,
            params.contrast_variable,
            params.reference_level,
            params.target_level,
            gtf_ch
        )

        // Se genera el reporte de Deseq2
        GENERATE_DESEQ2_REPORT(
            DESEQ2_DIFFERENTIAL.out.results,
            DESEQ2_DIFFERENTIAL.out.normalized_counts,
            file(params.samplesheet),
            file("$projectDir/modules/DESeq2_Interactive_Report.Rmd")
        )

        // === ANÁLISIS DE ENRIQUECIMIENTO ===
        if (params.gsea_run) {
            ch_deseq_results = DESEQ2_DIFFERENTIAL.out.results
                .map { file -> 
                    def meta = [id: file.simpleName]
                    tuple(meta, file)
                }
            
            PREPARE_GSEA_INPUT(ch_deseq_results)
            ch_gmt_files = Channel.fromPath(params.gsea_gmt_symbols).collect()
            GSEA(PREPARE_GSEA_INPUT.out.rnk, ch_gmt_files)
        }
        
        // === ANÁLISIS DE ONTOLOGIA DE GENES ===
        if (params.gprofiler_run) {
            ch_deseq_results = DESEQ2_DIFFERENTIAL.out.results
                .map { file ->
                    def meta = [id: file.simpleName]
                    tuple(meta, file)
                }
            GPROFILER(ch_deseq_results)
        }
        
        // === ANÁLISIS SUPPA ===
        if (params.run_suppa) {
            // Generación de eventos
            suppa2_generate_events_ch = suppa2_generate_events(file(params.gtf))

            event_files_ch = suppa2_generate_events_ch.events
                .flatten()
                .map { file -> 
                    def event_type = file.name.toString().find(/events_(.+)_strict\.ioe/)[1]
                    tuple(event_type, file)
                }

            psi_input_ch = event_files_ch
                .combine(TXIMPORT.out.suppa_tpm)

            // Cálculo de PSI
            suppa2_psi_per_event_ch = suppa2_psi_per_event(
                psi_input_ch.map { event_type, ioe_file, tpm -> ioe_file },
                psi_input_ch.map { event_type, ioe_file, tpm -> tpm }
            )

            // Procesamiento de TPM y PSI
            split_tpm_ch = split_tpm_file(TXIMPORT.out.suppa_tpm, file(params.samplesheet))
            split_psi_ioe_ch = split_psi_ioe(
                suppa2_psi_per_event_ch.psi.flatten(),
                file(params.samplesheet)
            )

            // Preparación para diffsplice IOE
            def event_types = ['A3', 'A5', 'AF', 'AL', 'MX', 'RI', 'SE']
            diffsplice_input_ioe_ch = Channel.fromList(event_types)
                .combine(split_psi_ioe_ch.collect())
                .combine(split_tpm_ch.split_tpm.collect())
                .map { it ->
                    def (event_type, psi_files, tpm_files) = [it[0], it[1..-3], it[-2..-1]]
                    def ioe_file = file("${params.outdir}/suppa/generate_events/events_${event_type}_strict.ioe")
                    def psi_control = psi_files.find { it.name.contains("${event_type}_strict_control") }
                    def psi_treatment = psi_files.find { it.name.contains("${event_type}_strict_treatment") }
                    def tpm_control = tpm_files.find { it.name.contains('control') }
                    def tpm_treatment = tpm_files.find { it.name.contains('treatment') }
                    
                    // Verificación de archivos
                    if (!ioe_file.exists()) {
                        log.warn "IOE file not found for event type ${event_type}: ${ioe_file}"
                        return null
                    }
                    if (!psi_control || !psi_treatment) {
                        log.warn "PSI files not found for event type ${event_type}: control=${psi_control}, treatment=${psi_treatment}"
                        return null
                    }
                    if (!tpm_control || !tpm_treatment) {
                        log.warn "TPM files not found for event type ${event_type}: control=${tpm_control}, treatment=${tpm_treatment}"
                        return null
                    }
                    
                    tuple(event_type, ioe_file, psi_control, psi_treatment, tpm_control, tpm_treatment)
                }
                .filter { it != null }

            // Análisis diffsplice IOE
            suppa_diffsplice_ioe(diffsplice_input_ioe_ch)

            // Análisis de isoformas
            suppa2_generate_isoforms_ch = suppa2_generate_isoforms(file(params.gtf))
            suppa_per_isoform_ch = suppa_per_isoform(
                file(params.gtf),
                TXIMPORT.out.suppa_tpm
            )

            // Split PSI para IOI
            split_psi_ioi_ch = split_psi_ioi(
                suppa_per_isoform_ch.psi_isoform,
                file(params.samplesheet)
            )
            
            // Preparación para diffsplice IOI
            diffsplice_input_ioi_ch = suppa2_generate_isoforms_ch.isoforms
                .combine(split_psi_ioi_ch.collect())
                .combine(split_tpm_ch.split_tpm.collect())
                .map { it ->
                    def (ioi_file, psi_files, tpm_files) = [it[0], it[1..-3], it[-2..-1]]
                    def psi_control = psi_files.find { it.name.contains("control") }
                    def psi_treatment = psi_files.find { it.name.contains("treatment") }
                    def tpm_control = tpm_files.find { it.name.contains('control') }
                    def tpm_treatment = tpm_files.find { it.name.contains('treatment') }
                    
                    // Verificación de archivos
                    if (!ioi_file.exists()) {
                        log.warn "IOI file not found: ${ioi_file}"
                        return null
                    }
                    if (!psi_control || !psi_treatment) {
                        log.warn "PSI IOI files not found: control=${psi_control}, treatment=${psi_treatment}"
                        return null
                    }
                    if (!tpm_control || !tpm_treatment) {
                        log.warn "TPM files not found: control=${tpm_control}, treatment=${tpm_treatment}"
                        return null
                    }
                    
                    tuple(ioi_file, psi_control, psi_treatment, tpm_control, tpm_treatment)
                }
                .filter { it != null }

            // Análisis diffsplice IOI
            suppa_diffsplice_ioi(diffsplice_input_ioi_ch)
        }
    }

    // Generación de reporte MultiQC
    multiqc_run(multiqc_files.collect())
}

// Mensajes de finalización
workflow.onComplete {
    log.info "Pipeline completado en: $workflow.complete"
    log.info "Estado de ejecución: ${workflow.success ? 'OK' : 'fallido'}"
    log.info "Tiempo de ejecución: $workflow.duration"
    log.info "Directorio de salida: ${params.outdir}"
}

// Mensaje de error
workflow.onError {
    log.error "Pipeline fallido: $workflow.errorMessage"
}
