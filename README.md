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


# Explícación de los parametros utilizados en cada herramienta

## FastQC
- Control de calidad 
   ```
  fastqc ${reads} --outdir . --threads ${task.cpus}
   ```
- Generación del archivo resumen (MultiQC)
   ```
  multiqc .
   ```
  Parámetros:
   - `reads`: Corresponde a los archivos de las lecturas en formato .fastq, en este caso se utiliza para los datos crudos y los que han pasado por el proceso de recorte (trimming).
   - `outdir`: Especifica el directorio donde se almacenarán los resultados. En este caso, el usuario define la ruta del directorio de salida. En base a esta ruta se creará una carpeta donde se guardarán los archivos generados por FastQC de los datos crudos y procesados.
   - `threads`: Indica el número de núcleos de CPU a utilizar para ejecutar el proceso de análisis. En este caso se utilizó el valor 4.
“.”: Indica que MultiQC y FastQC deben buscar los archivos de salida en el directorio actual para generar los reportes.

## Trimmomatic
El comando utilizado para realizar el trimming de las secuencias es:
- Para datos de secuenciación de lectura simple (SE)
  ```
   trimmomatic SE $reads ${sample_id}_trimmed.fastq.gz ILLUMINACLIP:${params.adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
  ```
- Para datos de secuenciación de lectura pareada (PE)
   ```
   trimmomatic PE ${reads[0]} ${reads[1]} ${sample_id}_1_paired_trimmed.fastq.gz ${sample_id}_1_unpaired_trimmed.fastq.gz ${sample_id}_2_paired_trimmed.fastq.gz ${sample_id}_2_unpaired_trimmed.fastq.gz ILLUMINACLIP:${params.adapters}:2:30:10      LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
  ```
- Parámetros:
   - `SE o PE`: Se ajusta el modo de ejecución de trimmomatic si los datos son de lectura simple (Single-End, SE) o lectura pareada (Paired-End, PE). En ese caso el parámetro utilizado fue Paired-End.
   - `reads`: En el caso de datos SE, corresponde al archivo .fastq de lectura simple. Para PE, se proporcionan dos archivos .fastq, uno para cada lectura (forward y reverse).
   - `sample_id`: Identificador de la muestra.
   -`ILLUMINACLIP`: Define los adaptadores que deben ser recortados. En este caso estaba definido en TruSeq3-PE. Los valores utilizados son 2 para el número de mismatches permitidos, 30 para el score de calidad mínima para cortar y penalización por desajustes de adaptadores de 10.
   - `LEADING`: Elimina las bases de baja calidad desde el extremo 5' de la lectura si su calidad es inferior a 3.
   - `TRAILING`: Elimina las bases de baja calidad desde el extremo 3' si su calidad es inferior a 3.
   - `SLIDINGWINDOW`: Aplica un corte si la calidad promedio de una ventana de 4 bases es menor a 15.
   - `MINLEN`: Este parámetro recorta de secuencias cuya longitud final sea menor a 36 nucleótidos.


## STAR
- Para realizar el alineamiento
   ```
   STAR --genomeDir $genomeDir --readFilesIn $read1 $read2 --readFilesCommand zcat --readMatesLengthsIn NotEqual --outFileNamePrefix ${prefix}_ --runThreadN ${task.cpus} --sjdbGTFfile $gtf --sjdbOverhang ${params.overhang} --                 alignSJDBoverhangMin ${params.sjdbOverhangMin} --outFilterScoreMinOverLread ${params.filterScore} --outFilterMatchNminOverLread ${params.filterScore} --outFilterMismatchNmax ${params.outFilterMismatchNmax} --outFilterMultimapNmax 20 --alignMatesGapMax 1000000 --outSAMattributes All --outSAMtype BAM Unsorted --outFilterType BySJout --twopassMode Basic --alignEndsType Local --alignIntronMax ${params.alignIntronMax} --quantMode GeneCounts
  ```

- Parámetros:

   - `readFilesIn`: Archivos de lecturas en formato .fastq.
   - `readFilesCommand`: Opción para descomprimir archivos fastq.gz.
   - `readMatesLengthsIn`: Indica si las longitudes de las lecturas pareadas pueden ser diferentes.
   - `outFileNamePrefix`: Prefijo para nombrar los archivos de salida.
   - `runThreadN`: Número de núcleos de CPU utilizados en la alineación.
   - `sjdbGTFfile`: Archivo GTF con anotaciones para guiar la alineación.
   - `sjdbOverhang`: Longitud del segmento de lectura. En este caso fueron 101 pares de bases.
   - `alignSJDBoverhangMin`: Longitud mínima de los fragmentos de lectura alineados en las junctions de splicing. En este caso, se establece en 8 el valor.
   - `outFilterScoreMinOverLread`: Umbral mínimo para la puntuación de alineación en relación con la longitud de la lectura. Se usa para filtrar alineaciones de baja calidad. En este caso, el valor es 0.66.
   - `outFilterMatchNminOverLread`: Umbral mínimo para la proporción de coincidencias en relación con la longitud de la lectura. También este valor se ha establecido en 0.66.
   - `outFilterMismatchNmax`: Número máximo permitido de desajustes en una lectura. En este caso, se permitió un máximo de 5 desajustes.
   - `outFilterMultimapNmax`: Número máximo permitido de mapeos múltiples para una lectura. En este caso serán 20 ubicaciones de mapeo.
   - `alignMatesGapMax`: Máxima distancia permitida entre lecturas emparejadas en la alineación. Este valor fue de 1.000.000 nucleótidos.
   - `outSAMattributes` All: Especifica que todos los atributos del alineamiento deben ser incluidos en los archivos de salida.
   - `outSAMtype BAM Unsorted`: Formato de los archivos de salida. En este caso, se generará un archivo BAM sin ordenar.
   - `outFilterType` BySJout: Tipo de filtrado aplicado a las lecturas basadas en la información de las junctions de splicing detectadas.
   - `twopassMode` Basic: Modo de alineación en dos fases
   - `alignEndsType` Local: Método de alineación de los extremos de las lecturas, permitiendo alineaciones locales en lugar de globales.
   - `alignIntronMax` 1000000: Longitud máxima permitida de los intrones para la alineación. El valor utilizado fue de 1,000,000 nucleótidos.
   - `quantMode` GeneCounts: Modo de cuantificación que cuenta las lecturas alineadas por gen.
 
## Samtools
- Ordenar el archivo BAM:
   ```
   samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam $bam
  ```
   - Parámetros:
      - `task.cpus`: Especifica el número de hilos de CPU a utilizar para en el proceso, en este caso se definió en 4 hilos.
      - `o`: Indica el nombre del archivo de salida que contendrá las lecturas ordenadas. El archivo final tiene el identificador de la muestra, y sorted.bam para indicar que el archivo está ordenado.
      bam: Especifica la ruta del archivo BAM a ordenar. 

- Indexar el archivo BAM
   ```
   samtools index ${sample_id}.sorted.bam
  ```
   - Parámetros:
      - `sample_id.sorted.bam: Especifica el archivo BAM en particular que es ordenado e indexado. 


## StringTie
- Merge de GTF con StringTie:
   ```
   stringtie --merge -G $reference_gtf -o stringtie_merged.gtf assembly_gtf_list.txt -p $task.cpus
  ```
   - Parámetros:
      - `assembly_gtf_list.txt`: Lista de todos los archivos .gtf en el directorio.
      - `merge`: Indica que StringTie debe realizar una fusión de múltiples archivos GTF en uno solo.
      - `G`: Especifica el archivo GTF proporciona anotaciones conocidas para guiar el proceso de fusión. 
      - `o`: Define el nombre del archivo GTF de salida, que contiene la anotación fusionada de las muestras.
      - `p`: Indica el número de hilos de CPU a utilizar para el proceso de fusión. En este caso se utilizó 4.

- Comparación con GffCompare:
   ```
   gffcompare -R -V -r $reference_gtf stringtie_merged.gtf
  ```
   - Parámetros:
      - `R`: Utiliza el archivo GTF de referencia para la comparación.
      - `V`: Produce un informe de las diferencias entre las anotaciones.
      - `r`: Especifica el archivo GTF de referencia.
      - `stringtie_merged.gtf`: El archivo GTF fusionado que se comparará contra el archivo de referencia.

- Ejecutar StringTie en archivos BAM:
   ```
   stringtie $bam -G $gtf -o ${sample_id}.gtf $strandedness -a 8 -p $task.cpus
  ```
  ```
   stringtie $bam -G $gtf -o ${sample_id}_for_DGE.gtf $strandedness -a 8 -e -p $task.cpus
  ```
   - Parámetros:
      - `$bam`: Especifica el archivo BAM de entrada con lecturas alineadas.
      - `G`: Utiliza el archivo GTF de referencia para guiar el ensamblaje.
      - `o`: Define el nombre del archivo GTF de salida.
      - `$strandedness`: Especifica la orientación de las lecturas. El valor puede ser fr-first stranded, fr-second strand y fr-unstranded, el cual fue utilizado para este proceso.
      - `a`: Establece el tamaño del fragmento para el ensamble. En este caso, fue de 8 bases.
      - `e`: Ejecuta el ensamblaje en modo de expresión.

- Convertir GFF a GTF con GffRead:
   ```
   gffread -E gffcmp.annotated.corrected.gff -T -o gffcmp.annotated.corrected.gtf
  ```
   - Parámetros:
      - `E`: Excluye las transcripciones extendidas (opcional), en este caso se usó un archivo que mejora las correcciones del archivo de anotación.
      - `T`: Convierte el archivo a formato GTF.
      - `o`: Archivo GTF de salida.

- Preparar matrices de conteo con PrepDE:
   ```
   prepDE.py -i sample_lst.txt -l ${params.readlength} -g ${run_prefix}_gene_count_matrix.csv -t ${run_prefix}_transcript_count_matrix.csv
  ```
   - Parámetros:
      - `i`: Archivo que contiene la lista de archivos GTF generados para cada muestra.
      - `l`: Longitud de las lecturas utilizadas en el análisis.
      - `g`: Archivo CSV de salida con la matriz de conteo de genes.
      - `t`: Archivo CSV de salida con la matriz de conteo de transcritos.


## BigWig
   ```
   bamCoverage -b $bam -o ${sample_id}.bw
  ```
- Parámetros:
   - `b`: Especifica el archivo BAM de entrada que contiene las lecturas alineadas.
   - `o`: Define el nombre del archivo de salida en formato BigWig (.bw). 

