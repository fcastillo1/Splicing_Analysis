process GFFREAD_TX2GENE {
    tag "$gtf"
    label 'process_low'
    publishDir "${params.outdir}/gffread", mode: 'copy'

    input:
    path gtf

    output:
    path "*.tx2gene.tsv", emit: tx2gene
    path "versions.yml", emit: versions

    script:
    """
    # Extraer transcript_id, gene_id y gene_name (si existe)
    awk -F'\\t' '\$3=="transcript" {
        split(\$9,a,";");
        tid=""; gid=""; name="";
        for(i in a) {
            gsub(/^ +/, "", a[i]); # Eliminar espacios iniciales
            if(a[i]~/transcript_id/) {
                tid=gensub(/.*transcript_id[ ="]+([^"]+)"?.*/, "\\\\1", "g", a[i]);
            }
            if(a[i]~/gene_id/) {
                gid=gensub(/.*gene_id[ ="]+([^"]+)"?.*/, "\\\\1", "g", a[i]);
            }
            if(a[i]~/gene_name|gene_symbol/) {
                name=gensub(/.*gene_name[ ="]+([^"]+)"?.*/, "\\\\1", "g", a[i]);
                if(name=="") {
                    name=gensub(/.*gene_symbol[ ="]+([^"]+)"?.*/, "\\\\1", "g", a[i]);
                }
            }
        }
        if(tid!="" && gid!="") print tid "\\t" gid "\\t" (name!="" ? name : gid);
    }' $gtf | sort -u > ${gtf.baseName}.tx2gene.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | head -n1)
    END_VERSIONS
    """
}
