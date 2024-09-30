# Splicing_Analysis
# Pipeline de Análisis de Splicing Diferencial

## Descripción General

Este pipeline de Nextflow está diseñado para realizar un análisis completo de splicing diferencial en datos de RNA-seq. Integra varias herramientas bioinformáticas para proporcionar un análisis robusto y detallado del splicing alternativo.

## Características Principales

1. Control de calidad de lecturas antes y después del procesamiento (FastQC)
2. Recorte de adaptadores (Trimmomatic)
3. Alineamiento de lecturas (STAR)
4. Cuantificación de transcritos (Salmon)
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
- Mínimo 100 GB de espacio en disco
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
  (PENDIENTE AGREGAR ALGO MÁS)

3. Asegúrate de tener los siguientes archivos de referencia:
   - Genoma de referencia (formato FASTA)
   - Archivo de anotación (formato GTF)
   - Archivo del transcriptoma (formato FASTA)
   - Índice de STAR en caso de no saber realizarlo puedes hacer: (PENDIENTE)

### Ejecución del Pipeline

Ejecuta el pipeline con el siguiente comando: 

```
nextflow run main.nf -profile <docker/singularity/standard> \
    --samplesheet samplesheet.csv \
    --genome /ruta/al/genoma.fa \
    --gtf /ruta/a/anotacion.gtf \
    --outdir resultados
```

### Parámetros Principales

- `--samplesheet`: Ruta al archivo CSV de muestras (requerido)
- `--outdir`: Directorio de salida (por defecto: 'results')
- `--genome`: Ruta al archivo FASTA del genoma de referencia
- `--gtf`: Ruta al archivo GTF de anotación
- `--star_index`: Ruta al índice de STAR (opcional)
- `--salmon_index`: Ruta al índice de Salmon (opcional)
- `--single_end`: Usar si las lecturas son single-end (por defecto: false)
- `--stranded`: Especifica la orientación de la hebra ('first-strand', 'second-strand', o 'unstranded')

Para ver todos los parámetros disponibles:
```
nextflow run main.nf --help
```

## Estructura de Salida

```
resultados/
├── fastqc/
│   ├── raw/
│   └── trimmed/
├── trimmomatic/
├── star/
├── salmon/
├── stringtie/
├── rmats/
│   └── sashimiplot/
├── suppa/
└── multiqc/
    └── multiqc_report.html
```

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

- **Error de memoria en STAR**: Ajusta el parámetro `--limitBAMsortRAM` en la configuración de STAR.
- **Pipeline se detiene inesperadamente**: Verifica los logs en el directorio `work/` para más detalles.
- **Errores de Docker**: Asegúrate de que los contenedores estén correctamente configurados y accesibles.

## Contribuir

1. Fork el repositorio
2. Crea una nueva rama (`git checkout -b feature/AmazingFeature`)
3. Realiza tus cambios y haz commit (`git commit -m 'Add some AmazingFeature'`)
4. Push a la rama (`git push origin feature/AmazingFeature`)
5. Abre un Pull Request

## Cita

Si utilizas este pipeline en tu investigación, por favor cítalo como:

```
[Reyes, Francisca]. (2024). Pipeline de Análisis de Splicing Diferencial. GitHub. https://github.com/fcastillo1/Splicing_Analysis
```


## Contacto

[Reyes, Francisca] - fcastillor.19@gmail.com
[Munita, Roberto] - robertomunita@gmail.com

URL del Proyecto: [https://github.com/fcastillo1/differential-splicing-pipeline](https://github.com/fcastillo1/Splicing_Analysis)

## Agradecimientos

- [STAR](https://github.com/alexdobin/STAR)
- [Salmon](https://combine-lab.github.io/salmon/)
- [StringTie](https://ccb.jhu.edu/software/stringtie/)
- [rMATS](http://rnaseq-mats.sourceforge.net/)
- [SUPPA2](https://github.com/comprna/SUPPA)
- [Nextflow](https://www.nextflow.io/)
- [nf-core](https://nf-co.re/)
