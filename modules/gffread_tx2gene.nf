process GFFREAD_TX2GENE {
    tag "$gtf"
    label 'process_low'
    publishDir "${params.outdir}/gffread", mode: 'copy'


    input:
    path gtf

    output:
    path "*.tx2gene.tsv", emit: tx2gene
    path "versions.yml" , emit: versions
    path "gffread_log.txt", emit: log

    script:
    """
    echo "Contenido del directorio de trabajo:" > gffread_log.txt
    ls -la >> gffread_log.txt

    echo "Verificando existencia del archivo GTF:" >> gffread_log.txt
    [ -f $gtf ] && echo "El archivo GTF existe" >> gffread_log.txt || echo "El archivo GTF no existe" >> gffread_log.txt

    echo "Primeras líneas del archivo GTF:" >> gffread_log.txt
    head -n 5 $gtf >> gffread_log.txt

    echo "Ejecutando comando para generar tx2gene..." >> gffread_log.txt

    awk -F'\\t' '\$3=="transcript" {
        split(\$9,a,";");
        tid=""; gid="";
        for(i in a) {
            if(a[i]~/transcript_id/) tid=gensub(/.*transcript_id "([^"]+)".*/, "\\\\1", "g", a[i]);
            if(a[i]~/gene_id/) gid=gensub(/.*gene_id "([^"]+)".*/, "\\\\1", "g", a[i]);
        }
        if(tid!="" && gid!="") print tid "\\t" gid;
    }' $gtf > ${gtf.baseName}.tx2gene.tsv

    echo "Comando completado." >> gffread_log.txt

    echo "Contenido del archivo tx2gene (primeras 10 líneas):" >> gffread_log.txt
    head -n 10 ${gtf.baseName}.tx2gene.tsv >> gffread_log.txt

    echo "Número total de líneas en tx2gene:" >> gffread_log.txt
    wc -l ${gtf.baseName}.tx2gene.tsv >> gffread_log.txt

    echo "Número de transcript_ids únicos:" >> gffread_log.txt
    cut -f1 ${gtf.baseName}.tx2gene.tsv | sort | uniq | wc -l >> gffread_log.txt

    echo "Número de gene_ids únicos:" >> gffread_log.txt
    cut -f2 ${gtf.baseName}.tx2gene.tsv | sort | uniq | wc -l >> gffread_log.txt

    if [ ! -s ${gtf.baseName}.tx2gene.tsv ]; then
        echo "Error: El archivo tx2gene está vacío" >> gffread_log.txt
        exit 1
    else
        echo "El archivo tx2gene se generó correctamente" >> gffread_log.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | head -n1)
    END_VERSIONS
    """
}
