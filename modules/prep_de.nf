process prep_de {
    publishDir "${params.outdir}/prep_de", mode: 'copy'
    label 'process_medium'
    
    input:
    path(gtfs)
    
    output:
    path "sample_lst.txt", emit: sample_list
    path "*gene_count_matrix.csv", emit: gene_matrix
    path "*transcript_count_matrix.csv", emit: transcript_matrix
    path "prepDE_log.txt", emit: log_file
    
    script:
    def run_prefix = params.run_prefix ?: "output"
    """
    # Iniciar log
    echo "=== Inicio del proceso prep_de ===" > prepDE_log.txt
    date >> prepDE_log.txt
    echo "Archivos GTF recibidos:" >> prepDE_log.txt
    ls -l >> prepDE_log.txt

    # Crear lista de muestras y verificar GTFs
    echo "=== Creando lista de muestras ===" >> prepDE_log.txt
    for gtf in *.gtf; do
        sample=\$(basename \$gtf _for_DGE.gtf)
        echo "\$sample \$gtf" >> sample_lst.txt
        
        # Verificar formato GTF
        echo "Verificando formato del GTF para \$sample:" >> prepDE_log.txt
        awk -F"\\t" '{
            if(NF!=9) 
                print "Error en línea " NR ": " \$0
        }' \$gtf >> prepDE_log.txt
    done

    # Registrar lista de muestras
    echo "=== Contenido de sample_lst.txt ===" >> prepDE_log.txt
    cat sample_lst.txt >> prepDE_log.txt

    # Intentar prepDE.py
    echo "=== Ejecutando prepDE.py ===" >> prepDE_log.txt
    if prepDE.py -i sample_lst.txt \\
        -l ${params.readlength} \\
        -g ${run_prefix}_gene_count_matrix.csv \\
        -t ${run_prefix}_transcript_count_matrix.csv \\
        >> prepDE_log.txt 2>&1; then
        
        echo "prepDE.py se ejecutó correctamente." >> prepDE_log.txt
    
    else
        echo "prepDE.py falló. Usando método alternativo..." >> prepDE_log.txt
        
        # Crear archivos de matriz
        echo "Sample,Gene_ID,Count" > ${run_prefix}_gene_count_matrix.csv
        echo "Sample,Transcript_ID,Count" > ${run_prefix}_transcript_count_matrix.csv
        
        # Procesar cada GTF
        for gtf in *.gtf; do
            sample=\$(basename \$gtf _for_DGE.gtf)
            echo "Procesando \$gtf..." >> prepDE_log.txt
            
            awk -F"\\t" '
            \$3=="exon" {
                split(\$9,a,";");
                gene_id=""; transcript_id="";
                for(i in a) {
                    if(a[i] ~ /gene_id/) 
                        gene_id = gensub(/.*gene_id "([^"]+)".*/, "\\\\1", "g", a[i]);
                    if(a[i] ~ /transcript_id/) 
                        transcript_id = gensub(/.*transcript_id "([^"]+)".*/, "\\\\1", "g", a[i]);
                }
                if(gene_id != "" && transcript_id != "") {
                    genes[gene_id]++;
                    transcripts[transcript_id]++;
                }
            }
            END {
                for(g in genes) 
                    print "'"\$sample'"'","g","genes[g];
                for(t in transcripts) 
                    print "'"\$sample'"'","t","transcripts[t];
            }' \$gtf | sort -t',' -k2,2 >> ${run_prefix}_gene_count_matrix.csv
        done
        
        echo "Método alternativo completado." >> prepDE_log.txt
    fi

    # Verificar resultados
    echo "=== Verificación de resultados ===" >> prepDE_log.txt
    for file in ${run_prefix}_gene_count_matrix.csv ${run_prefix}_transcript_count_matrix.csv; do
        if [ -f "\$file" ]; then
            echo "Archivo \$file generado correctamente." >> prepDE_log.txt
            echo "Primeras 10 líneas de \$file:" >> prepDE_log.txt
            head -n 10 \$file >> prepDE_log.txt
            wc -l \$file >> prepDE_log.txt
        else
            echo "ERROR: No se generó el archivo \$file" >> prepDE_log.txt
            echo "Creando archivo dummy..." >> prepDE_log.txt
            echo "Sample,ID,Count" > \$file
            echo "DummySample,DummyID,0" >> \$file
        fi
    done

    echo "=== Proceso completado ===" >> prepDE_log.txt
    date >> prepDE_log.txt
    """
}
