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

   ```

2. Instala Nextflow:
   ```
   curl -s https://get.nextflow.io | bash
   ```

3. Si usas Docker, asegúrate de tener Docker instalado y en ejecución.
   - Verificar la instalación de Docker
   ```
   docker --version
   ```
   - Ver dockers activos
   ```
   docker ps
   ```

   - Ver todos los dockers (incluyendo detenidos) 
   ```
   docker ps -a 
   ```
   
   - Verificar espacio disponible
   ```
   df -h
   ```
   
   - Los montajes en Docker conectan directorios de tu computadora (host) con directorios dentro del contenedor, es decir para poder utilizar los archivos correctamente y que estén disponibles con los programas de docker se tienen que montar para esto se utiliza la opción -v, como se observa a continuación
   ```
   -v /ruta/en/tu/computadora:/ruta/en/el/contenedor
   ```

   - En el directorio del proyecto donde esta el main y el docker hay que construir la imagen para usar los programas
   ```
   docker build -t rnaseq_pipeline:latest .
   ```
   - Después de ejecuta el pipeline 
   ```
   sudo docker run \
     # Montajes de datos de entrada
     -v /home/francisca/Data_Intento_Giulia:/Data_Giulia \
     -v /mnt/disco_2/FReyes/Giulia_Data:/Giulia_Data \
     
     # Montajes de datos de referencia
     -v /mnt/disco_2/FReyes/Datos:/FReyes/Datos \
     
     # Montaje para resultados
     -v /mnt/disco_2/FReyes/results_control_treatment_star:/FReyes/results_control_treatment_star \
     
     # Montajes del pipeline
     -v /home/francisca/Proyecto/DifferentialSplicingAnalysis:/pipeline \
     -v $(pwd)/work:/pipeline/work \
     
     # Directorio de trabajo
     -w /pipeline \
     
     # Imagen y comando
     rnaseq_pipeline:latest \
     nextflow run main.nf \
     [resto de los parámetros...]
   ```

## Configuración

1. Archivo de configuración principal: `nextflow.config`
   - Aquí se definen los parámetros globales y las configuraciones de los procesos.

2. Archivos de configuración de Docker (Dockerfile)

4. Módulos: Directorio `modules/`
   - Contiene los scripts de los procesos individuales (STAR, Salmon, rMATS, etc.)

## Uso

### Preparación de Datos

1. Prepara tu archivo de muestras (`samplesheet.csv`) con el siguiente formato:
- Para descargas (SRR):
   ```
   # Paired-end SRR
   sample_id,fastq,fastq2,condition,trimmed_fastq,trimmed_fastq2
   SRR12463396,,,control,,,
   SRR12463397,,,treatment,,,

   # Single-end SRR
   sample_id,fastq,fastq2,condition,trimmed_fastq,trimmed_fastq2
   SRR12463396,,,control,,
   SRR12463397,,,treatment,,
   ```

- Para archivos crudos locales:
  ```
   # Paired-end crudos
   sample_id,fastq,fastq2,condition,trimmed_fastq,trimmed_fastq2
   sample1,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz,control,,,
   sample2,/path/to/sample2_R1.fastq.gz,/path/to/sample2_R2.fastq.gz,treatment,,,
   
   # Single-end crudos
   sample_id,fastq,fastq2,condition,trimmed_fastq,trimmed_fastq2
   sample1,/path/to/sample1.fastq.gz,,control,,
   sample2,/path/to/sample2.fastq.gz,,treatment,,
   ```

- Para archivos ya trimados:
  ```
   # Paired-end trimados
   sample_id,fastq,fastq2,condition,trimmed_fastq,trimmed_fastq2
   sample1,,,control,/path/to/sample1_1_trimmed.fastq.gz,/path/to/sample1_2_trimmed.fastq.gz
   sample2,,,treatment,/path/to/sample2_1_trimmed.fastq.gz,/path/to/sample2_2_trimmed.fastq.gz
   
   # Single-end trimados
   sample_id,fastq,fastq2,condition,trimmed_fastq,trimmed_fastq2
   sample1,,,control,/path/to/sample1_trimmed.fastq.gz,
   sample2,,,treatment,/path/to/sample2_trimmed.fastq.gz,
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

## Procesos Detallados

1. **FastQC**: Control de calidad de lecturas crudas y recortadas.
2. **Trimmomatic**: Recorte de adaptadores y filtrado de calidad.
3. **STAR**: Alineamiento de lecturas al genoma de referencia.
4. **Salmon**: Cuantificación de la expresión de transcritos.
6. **rMATS**: Detección y cuantificación de eventos de splicing diferencial.
7. **SUPPA2**: Análisis adicional de splicing diferencial.
8. **DESeq2**: Análisis adicional de expresión diferencial.
9. **GProfiler2**: Enriquecimiento funcional de genes y rutas biológicas.
10. **GSEA**: Enriquecimiento de genes
11. **MultiQC**: Generación de informe de control de calidad integrado.

## Solución de Problemas
- **Pipeline se detiene inesperadamente**: Verifica los logs en el directorio `work/` para más detalles.
- **Errores de Docker**: Asegúrate de que los contenedores estén correctamente configurados y accesibles.

## Cita

Si utilizas este pipeline en tu investigación, por favor cítalo como:

```
[Reyes, Francisca]. (2024). Pipeline de Análisis de Splicing Diferencial. GitHub. LINK
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
 
## rMATS
   ```
   rmats.py --b1 <bam1_file> --b2 <bam2_file> --gtf <reference_gtf> --od <output_directory> --tmp <temporary_directory> -t paired --readLength 101 --nthread 12 --libType fr-firststrand --novelSS --mil 50 --mel 500
  ```
- Parámetros:
   - `--b1`: Archivo BAM para el primer grupo de muestras.
   - `--b2`: Archivo BAM para el segundo grupo de muestras.
   - `--gtf`: Archivo GTF de referencia.
   - `--od`: Directorio donde se almacenarán los resultados de rMATS.
   - `--tmp`: Directorio para archivos temporales.
   - `-t`: Especifica que las lecturas son emparejadas, en este caso se indica que es paired-end.
   - `--readLength`: Longitud de las lecturas, el valor de este parámetro corresponde a 101.
   - `--nthread`: Número de hilos de CPU a utilizar para el procesamiento, se utilizó un valor de 12.
   - `--libType`: Tipo de biblioteca, en este caso, fr-firststrand para lecturas orientadas al primer extremo.
   - `--novelSS`: Permite la detección de sitios de empalme novedosos.
   - `--mil`: Número mínimo de lecturas para considerar un evento de splicing, en este caso se utilizó el valor de 50.
   - `--mel`: Número máximo de lecturas para considerar un evento de splicing, se definió en 500 para este proceso.


## rMATS (sashimiplot)
   ```
   rmats2sashimiplot --b1 $bam1 --b2 $bam2 --event-type $event_type -e "${gene}_${event_type}_events_with_header.txt" --exon_s 1 --intron_s 5 --l1 $label1 --l2 $label2 -o "${gene}_${event_type}_sashimi_plot" --group-info $grouping_file --   font_size 8 --fig-height 8 --fig-width 12
  ```

- Parámetros:
   - `--b1`: Archivo BAM para el primer grupo de muestras.
   - `--b2`: Archivo BAM para el segundo grupo de muestras.
   - `--event-type`: Tipo de evento de empalme a visualizar, las opciones a utilizar son SE, A5SS, A3SS, MXE, RI. En este caso para efectos del proyecto se realizan los sashimi plots solo con SE (Exon Skipping).
   - `--e`: Archivo con los eventos de empalme, generado por rMATS.
   - `--exon_s`: Tamaño del exón en el gráfico, el valor utilizado en este caso es 1.
   - `--intron_s`: Tamaño del intrón en el gráfico, el valor considerado en este caso es 5.
   - `--l1`: Etiqueta para el primer grupo de muestras.
   - `--l2`: Etiqueta para el segundo grupo de muestras.
   - `--o`: Nombre de los archivos de salida.
   - `--group-info`: Archivo para especificar las condiciones de agrupamiento en los gráficos. Se debe indicar el rango de las muestras para cada condición. Por ejemplo, si hay tres muestras de control y tres de tratamiento, el archivo será:
      control: 1-3
      tratamiento: 4-6
   - `--font-size`: Tamaño de fuente para el gráfico, en ese caso el valor definido fue de 8.
   - `--fig-height`: Altura del gráfico en pulgadas, el valor utilizado es de 8.
   - `--fig-width`: Ancho del gráfico en pulgadas, el valor en este caso fue de 12.

## Salmon 
- Salmon Index
   ```
   salmon index --threads $task.cpus -t gentrome.fa.gz -d decoys.txt -i salmon_index
  ```
   - Parámetros:
      - `--threads`: Número de hilos de CPU a utilizar, en este caso se definió con 12.
      - `-t`: Archivo de transcriptoma en formato FASTA o FASTQ, este se combina los transcritos y el genoma
      - `-d`: Archivo que contiene secuencias de decoy, las cuales son secuencias genómicas usadas para evitar alineaciones erróneas a regiones no transcriptómicas, mejorando la precisión de la cuantificación de transcritos. En este caso el archivo se generó por defecto.
      - `-i`: Directorio donde se almacenará el índice de Salmon.

- Salmon quant
   ```
   salmon quant --index $index --libType $strandedness $input_reads --geneMap $gtf --threads $task.cpus --validateMappings --numBootstraps 100 --seqBias --gcBias $args -o $prefix
  ```
   - Parámetros:
      - `--index`: Directorio que contiene el índice de Salmon.
      - `--libType`: Existen cuatro tipos de bibliotecas: ISR (first-strand, reverse), ISF (second-strand, forward), IU (unstranded), y A(unknown). En este caso, se usaron datos unstranded (IU) según el estudio de Ghandi et al. (2019). Para confirmar, se realizó un análisis con inferExperiment de BEDTools sobre 200,000 lecturas paired-end del archivo SRR8615452_Aligned.sortedByCoord.out, de las cuales un 3.86% no pudo determinarse su orientación. El análisis mostró que las fracciones de lecturas en las categorías para definir la hebra a la cual pertenecen son: 1++, 1--, 2+-, 2-+ (0.4808) y 1+-, 1-+, 2++, 2--(0.4806). Estos valores al ser muy similares confirmaron que el experimento es unstranded.
      - `--input_reads`: Archivos de lectura de entrada
      - `--geneMap`: Archivo GTF para mapear los genes
      - `--threads`: Número de hilos de CPU a utilizar 
      - `--validateMappings`: Activa la validación de mapeos.
      - `--numBootstraps`: Número de réplicas de bootstrap para estimación de errores, el valor utilizado fue 100.
      - `--seqBias`: Corrige sesgos específicos en los datos de entrada.
      - `--gcBias`: Corrige sesgos relacionados de contenido GC.
      - `-o`: Directorio de salida para los resultados.

## Suppa
- Generación de eventos.
   ```
   suppa.py generateEvents -i ${gtf_file} -o events -f ioe -e ${params.events}
  ```
   - Parámetros:
      - `-i`: Archivo GTF de entrada.
      - `-o`: Prefijo para los archivos de salida.
      - `-f`: Formato de eventos, en este caso se utilizó IOE (Inclusion of Events), donde contiene información sobre la inclusión o exclusión de un exón en eventos de splicing alternativo.
      - `-e`: Tipos de eventos a generar, entre las opciones son SE, A5/A3, MX, RI, AF/AL. En este caso se consideraron todos los eventos.

- Generación de isoformas.
   ```
   suppa.py generateEvents -i ${gtf_file} -o isoforms -f ioi
  ```
   - Parámetros:
      - `-i`: Archivo GTF de entrada.
      - `-o`: Prefijo para los archivos de salida.
      - `-f`: Formato de eventos, la opción IOI (Inclusion of Isoforms) indica la proporción de isoformas en los eventos de splicing y permite medir cómo influyen en las variantes de transcritos.

- Cálculo de PSI por evento.
   ```
   suppa.py psiPerEvent -i $ioe_file -e $tpm -o ${ioe_file.baseName}
  ```
   - Parámetros:
      - `-i`: Archivo IOE con eventos.
      - `-e`: Archivo TPM con valores obtenidos desde tximport.
      - `-o`: Prefijo para los archivos de salida.

- Cálculo de PSI por isoforma.
  ```
   suppa.py psiPerIsoform -g $gtf -e $tpm -o per_isoform
  ```
   - Parámetros:
      - `-g`: Archivo GTF de referencia.
      - `-e`: Archivo TPM con valores obtenidos desde tximport.
      - `-o`: Prefijo para los archivos de salida.


- Análisis de splicing diferencial.
  ```
  suppa.py diffSplice -m empirical -i ${ioe_file} -p ${psi_basalA} ${psi_basalB} -e ${tpm_basalA} ${tpm_basalB} --area ${params.area} --lower-bound ${params.lower_bound} --alpha ${params.alpha} --tpm-threshold ${params.tpm_threshold} --nan-   threshold ${params.nan_threshold} -gc -o diffsplice_${event_type}
  ```
   - Parámetros:
      - `-m`: Método para el cálculo diferencial, puede ser empirical/classical. En este caso se utilizó empirical.
      - `-i`: Archivo IOE con evento.
      - `-p`: Archivos con valores PSI para las diferentes condiciones.
      - `-e`: Archivos con valores TPM por cada una de las condiciones.
      - `--area`: Área alrededor del evento, el valor utilizado corresponde a 1000.
      - `--lower-bound`: Umbral inferior de significancia, el valor de este parámetro fue de 0.05.
      - `--alpha`: Nivel de significancia, el valor utilizado es 0.05
      - `--tpm-threshold`: Umbral de TPM, el valor definido es 0.
      - `--nan-threshold`: Umbral para valores NaN, el valor utilizado es 0.
      - `-gc`: Corrección de sesgo GC.
      - `-o`: Prefijo para los archivos de salida.

- Clustering de eventos de splicing.
  ```
  suppa.py clusterEvents --dpsi diffsplice_${event_type}.dpsi --psivec diffsplice_${event_type}.psivec --sig-threshold ${params.sig_threshold} --dpsi-threshold ${params.dpsi_threshold} --eps ${params.eps} --metric ${params.metric} --min-pts    ${params.min_pts} --groups ${params.groups} --clustering ${params.clustering_method} -o cluster${event_type}
  ```
   - Parámetros:
      - `--dpsi`: Archivo ΔPSI para clustering.
      - `--psivec`: Archivo vector PSI para clustering.
      - `--sig-threshold`: Umbral de significancia, el valor definido fue de 0.05.
      - `--dpsi-threshold`: Umbral de ΔPSI, en este caso el valor utilizado fue 0.05.
      - `--eps`: Parámetro épsilon para DBSCAN, el valor fue 0.05.
      - `--metric`: Métrica de las distancias utilizadas para el cálculo, las opciones son euclidean, manhattan, cosine. En este caso la métrica utilizada fue euclidean.
      - `--min-pts`: Número mínimo de puntos para clúster, en este caso el valor utilizado es 20.
      - `--groups`: Para el clustering de eventos, especifica los rangos de columnas para cada condición. Por ejemplo,  si la primera condición tiene tres muestras y la segunda condición también, el parámetro sería 1-3,4-6.
      - `--clustering`: Método de clustering, existen dos opciones las cuales son DBSCAN, OPTICS. En este caso el valor utilizado fue DBSCAN.
      - `-o`: Prefijo para los archivos de salida.


## DESeq2
  ```
  Rscript $projectDir/modules/deseq_differential.R --count_file '${counts}' --sample_file '${samplesheet}' --contrast '${contrast}' --output_prefix '${prefix}' --gtf_file '${gtf_file}' --lfc_threshold ${params.lfc_threshold ?: 0} --alpha ${params.alpha ?: 0.1} --p_adjust_method '${params.p_adjust_method ?: 'BH'}' --shrink_lfc ${params.shrink_lfc ?: 'true'} --cores ${task.cpus}
  ```
- Parámetros:
   - `--count`: Archivo con las cuentas de lectura, en este caso se utilizó el archivo salmon.merged.gene_counts_length_scaled.tsv proveniente desde tximport.
   - `--samplesheet`: Archivo con la información de las muestras.
   - `--contrast`: Contraste para la comparación, en este caso se definió como condiciones, con el parámetro de reference_level para la primera condición de las lecturas y target_level para la segunda condición a comparar.
   - `--output_prefix`: Prefijo para los archivos de salida.
   - `--gtf_file`: Archivo GTF de anotación.
   - `--lfc_threshold`: Umbral para el cambio logarítmico de la expresión (LFC), el valor definido fue de 0.58.
   - `--alpha`: Nivel de significancia, el valor utilizado fue de 0.05.
   - `--p_adjust_method`: Método de ajuste para p-valores. Las opciones corresponden a: Sin ajuste (none), utilizando el método de Benjamini-Hochberg (BH), método de Benjamini-Yekutieli para pruebas correlacionadas (BY), método de Holm-Bonferroni de ajuste secuencial (HOLM) y utilizando el método de Hochberg con un ajuste secuencial y conservador. En este caso se utilizó el método de BH.
   - `--shrink_lfc`: Indica si se debe realizar el ajuste los cambios en el logaritmo para reducir la varianza en el valor de Log FoldChange, en este caso el parámetro fue true.
   - `--cores`: Número de núcleos a utilizar, el valor utilizado fue 4.


## GSEA

  ```
  $gseaPath/gsea-cli.sh GSEAPreranked -rnk $rnk -gmx $gmt_file -out $prefix -set_min ${params.gsea_set_min} -set_max ${params.gsea_set_max} -collapse ${params.gsea_collapse} -mode ${params.gsea_mode} -create_svgs ${params.gsea_create_svgs}_include_only_symbols ${params.gsea_include_only_symbols} -make_sets ${params.gsea_make_sets} -plot_top_x ${params.gsea_plot_top_x} -rnd_seed ${params.gsea_rnd_seed} -zip_report ${params.gsea_zip_report} $args
  ```
- Parámetros:
   - `GSEAPreranked`: Permite ejecutar el análisis GSEA con pre-ranking.
   - `-rnk`: Archivo de ranking fue generado con los resultados de DESeq2.
   - `-gmx`: Archivo con los conjuntos de genes. En este caso se descargaron desde la base de datos MSigDB para Homo Sapiens. Los archivos utilizados fueron Hallmarks de cáncer (H), Pathways de categorías de genes curados (C2), Targets reguladores de genes (C3), Ontología de genes (C5), Pathways de oncogenes (C6) y Conjuntos de genes característicos de cada célula (C8).
   - `-out`: Prefijo para los archivos de salida.
   - `-set_min`: Tamaño mínimo de los conjuntos, en este caso fue de 10.
   - `-set_max`: Tamaño máximo de los conjuntos, el valor utilizado fue 1000.
   - `-collapse`: Indica si se colapsaron los conjuntos, en este caso fue false.
   - `-mode`: Modo de análisis, se utilizó con máxima probabilidad.
   - `-create_svgs`: Indica la creación de gráficos SVG, en este caso fue true.
   - `-include_only_symbols`: Incluye la simbología de genes, se utilizó true.
   - `-make_sets`: Indica la creación de conjuntos de genes, en este caso fueron creados.
   - `-plot_top_x`: Número de conjuntos para mostrar en los gráficos, el valor utilizado fue de 20.
   - `-rnd_seed`:  Generación de números aleatorios, se utilizó timestamp.
   - `-zip_report`: Indica si se debe comprimir el informe, el parámetro empleado fue false.

## GProfiler

  ```
  gost_results <- gost(query = de_genes, organism = "hsapiens", correction_method = "g_SCS", user_threshold = 0.05, domain_scope = "custom", custom_bg = background_genes, sources = c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC", "TF", "MIRNA", "HPA", "CORUM"))
  ```
- Parámetros:
   - `query`: Listado de genes diferenciales para analizar.
   - `organism`: Organismo de interés, en este caso fue Homo sapiens.
   - `correction_method`: Método de corrección para p-valores se decidió utilizar g_SCS, el cual es un método que ajusta los valores de p-value en función de las puntuaciones de enriquecimiento en una condición.
   - `padj_threshold`: Umbral para significancia, el valor utilizado fue de 0.05
   - `lfc_treshold`: Es el valor límite para filtrar eventos de splicing según su log Fold Change, en este caso fue definido con el valor 0.
   - `domain_scope`: Uso de un background definido de genes.
   - `sources`: Fuentes de información de los términos de enriquecimiento a incluir pueden ser biología molecular, procesos biológicos, etc. Las opciones disponibles son: GO:BP, GO:MF, GO:CC, KEGG, REAC, TF, MIRNA, HPA, CORUM todas estas fueron incluidas en este proceso.
