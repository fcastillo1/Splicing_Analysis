# Splicing_Analysis
# Pipeline de Análisis de Splicing Diferencial

## Descripción General

Este pipeline de Nextflow está diseñado para realizar un análisis completo de splicing diferencial en datos de RNA-seq. Integra varias herramientas bioinformáticas para proporcionar un análisis robusto y detallado del splicing alternativo.

## Características Principales

1. Control de calidad de lecturas antes y después del procesamiento (FastQC)
2. Recorte de adaptadores (Trimmomatic)
3. Alineamiento de lecturas (STAR)
4. Alineamiento y Cuantificación de transcritos (Salmon)
5. Reconstrucción del transcriptoma (StringTie)
6. Análisis de splicing diferencial (rMATS y SUPPA2)
7. Visualización de eventos de splicing (Sashimi plots)
8. Estimación de expresión a nivel de genes (Tximport)
9. Análisis de expresión diferencial (DESeq, GSEA y Gprofiler)
10. Informe de control de calidad integrado (MultiQC)

## Requisitos del Sistema

- Nextflow (v20.04.0 o superior)
- Java (v8 o superior)
- Docker
- Mínimo 64 GB de RAM
- Mínimo 128 GB de espacio en disco
- Sistema operativo: Linux (recomendado), macOS

## Instalación

1. Clona este repositorio:
   ```
   git clone https://github.com/fcastillo1/Splicing_Analysis.git
   cd Splicing_Analysis
   ```

2. Instala Nextflow:
   ```
   curl -s https://get.nextflow.io | bash
   ```

3. Si usas Docker, asegúrate de tener Docker instalado y en ejecución. (AGREGAR PASOS)


## Configuración

1. Archivo de configuración principal: `nextflow.config`
   - Aquí se definen los parámetros globales y las configuraciones de los procesos.

2. Archivos de configuración de Docker (PENDIENTE)

3. Módulos: Directorio `modules/`
   - Contiene los scripts de los procesos individuales (STAR, Salmon, rMATS, etc.)

## Uso

### Preparación de Datos

1. Prepara tu archivo de muestras (`samplesheet.csv`) con el siguiente formato:
   - En caso de que las muestras sean Single End
     ```
     sample_id,fastq,fastq2,condition
     muestra1,/ruta/a/muestra1_R1.fastq.gz,,control
     muestra2,/ruta/a/muestra2_R1.fastq.gz,,tratamiento
     ```
    - En caso de que las muestras sean Paired End
       ```
       sample_id,fastq,fastq2,condition
       muestra1,/ruta/a/muestra1_R1.fastq.gz,/ruta/a/muestra1_R2.fastq.gz,control
       muestra2,/ruta/a/muestra2_R1.fastq.gz,/ruta/a/muestra2_R2.fastq.gz,tratamiento
       ```
    - En caso de que las muestras necesiten ser descargadas
       ```
       sample_id,fastq,fastq2,condition
       SRR muestra1,,,control
       SRR muestra2,,,tratamiento
       ```

2. Asegúrate de tener los siguientes archivos de referencia:
   - Genoma de referencia (formato FASTA) puede ser comprimido.
   - Archivo de anotación (formato GTF) tiene que estar en formato gtf no .gz
   - Archivo del transcriptoma (formato FASTA) puede estar comprimido.
   - Índice de STAR en caso de no saber realizarlo puedes hacer:
        ```
      STAR --runThreadN X -runMode genomeGenerate \
           --genomeDir /path/star_index \
           --genomeFastaFiles /path/archivo.genome.fa \
           --sjdbGTFfile path/archivo.anotacion.gtf \
           --sjdbOverhang largo_pares_bases -1
       ```
3. Crear archivos necesarios para el funcionamiento del pipeline
   - Archivo grouping.nf el cual es necesario para hacer las comparaciones con rMATS2sashimiplot:
       ```
       condicion1: 1-3
       condicion2: 4-6
       ```
      Este archivo se genera con el número de las muestras (archivos BAM) separados por un guión. Este orden debe seguir a los archivos BAM.
     
   - Archivo constasts.csv para realizar las comparaciones con DESeq2. Por ejemplo:
       ```
      id,variable,reference,target,blocking
      basalB_vs_basalA,condition,basalA,basalB,
       ```
   - Generación de archivos BAM por cada condición para realizar el proceso de rMATS
       ```
      condicion 1: path/archivo1.bam,path/archivo2.bam,path/archivo3.bam
      condicion 2: path/archivo4.bam,path/archivo5.bam,path/archivo6.bam
       ```
       
### Ejecución del Pipeline

Ejecuta el pipeline con el siguiente comando: 

```
nextflow run main.nf -profile <docker/singularity/standard> \
    --samplesheet samplesheet.csv \
    --transcript_fasta \
    --genome /path/genoma.fa \
    --gtf /path/anotacion.gtf \
    --grouping /path/grouping.nf \
    --aligner_mode \
    --outdir resultados
```
Algunos parametros pueden ser modificados en nextflow.config

### Parámetros Principales

- `--samplesheet`: Ruta al archivo CSV de muestras (requerido)
- `--outdir`: Directorio de salida (por defecto: 'results')
- `--genome`: Ruta al archivo FASTA del genoma de referencia
- `--transcript_fasta`: Ruta al archivo FASTA de los transcritos
- `--gtf`: Ruta al archivo GTF de anotación
- `--star_index`: Ruta al índice de STAR
- `--single_end`: Usar si las lecturas son single-end (por defecto: false)
- `--stranded`: Especifica la orientación de la hebra ('first-strand', 'second-strand', o 'unstranded')
- Los modos de alineamientos corresponden a:
     - `aligner_star`: Realiza todo el proceso del pipeline del brazo de STAR. Sus opciones son True/False.
     - `aligner_salmon`: Realiza todo el proceso del pipeline del brazo de Salmon. Sus opciones son True/False.
     - `full_pipeline`: Realiza todo el pipeline. Sus opciones son True/False.
- `readlength`: Largo de los reads.
- `adapters`: Ruta a los adaptadores.

Para ver todos los parámetros disponibles se puede acceder a nextflow.config



## Procesos Detallados (PENDIENTE)

1. **FastQC**: Control de calidad de lecturas crudas y recortadas.
2. **Trimmomatic**: Recorte de adaptadores y filtrado de calidad.
3. **STAR**: Alineamiento de lecturas al genoma de referencia.
4. **Salmon**: Cuantificación de la expresión de transcritos.
5. **StringTie**: Ensamblaje de transcritos y estimación de abundancia.
6. **rMATS**: Detección y cuantificación de eventos de splicing diferencial.
7. **SUPPA2**: Análisis adicional de splicing diferencial.
8. **MultiQC**: Generación de informe de control de calidad integrado.

## Solución de Problemas
- **Pipeline se detiene inesperadamente**: Verifica los logs en el directorio `work/` para más detalles.
- **Errores de Docker**: Asegúrate de que los contenedores estén correctamente configurados y accesibles.

## Cita

Si utilizas este pipeline en tu investigación, por favor cítalo como:

```
[Reyes, Francisca]. (2024). Pipeline de Análisis de Splicing Diferencial. GitHub. https://github.com/fcastillo1/Splicing_Analysis
```

## Contacto

- [Reyes, Francisca] - fcastillor.19@gmail.com
- [Munita, Roberto] - robertomunita@gmail.com

- URL del Proyecto: [https://github.com/fcastillo1/differential-splicing-pipeline](https://github.com/fcastillo1/Splicing_Analysis)

## Agradecimientos

- [FastQC] (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [Trimmomatic] (http://www.usadellab.org/cms/?page=trimmomatic)
- [STAR] (https://github.com/alexdobin/STAR)
- [Salmon] (https://combine-lab.github.io/salmon/)
- [Samtools] (http://www.htslib.org/)
- [StringTie] (https://ccb.jhu.edu/software/stringtie/)
- [rMATS] (http://rnaseq-mats.sourceforge.net/)
- [Sashimi plots] (https://github.com/guigolab/rmats2sashimiplot)
- [Tximport] (https://bioconductor.org/packages/release/bioc/html/tximport.html)
- [SUPPA2] (https://github.com/comprna/SUPPA)
- [DESeq2] (https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
- [GSEA] (https://www.gsea-msigdb.org/gsea/index.jsp)
- [g:Profiler] (https://biit.cs.ut.ee/gprofiler/gost)
- [MultiQC] https://multiqc.info/)
- [Nextflow] (https://www.nextflow.io/)
- [nf-core] (https://nf-co.re/)

  
