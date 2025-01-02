process GPROFILER {
   tag "gProfiler on ${meta.id}"
   label 'process_medium'
   
   publishDir "${params.outdir}/gprofiler/${meta.id}", mode: 'copy'
   
   input:
   tuple val(meta), path(deseq_results)
   
   output:
   tuple val(meta), path("*_gprofiler_results.tsv"), emit: results
   tuple val(meta), path("*_enrichment_plots.pdf"), emit: plots
   tuple val(meta), path("*_summary.txt"), emit: summary
   tuple val(meta), path("*_enrichment_report.html"), emit: report
   path "versions.yml", emit: versions
   
   script:
   """
   #!/usr/bin/env Rscript

   # Cargar librerías con mensajes suprimidos
   suppressPackageStartupMessages({
       library(gprofiler2)
       library(tidyverse)
       library(ggplot2)
       library(igraph)
       library(RColorBrewer)
       library(ggridges)
       library(viridis)
   })

   # Configurar tema general
   theme_set(theme_bw() + theme(
       text = element_text(size = 12),
       plot.title = element_text(size = 14, face = "bold"),
       plot.subtitle = element_text(size = 10, face = "italic"),
       axis.text.x = element_text(angle = 45, hjust = 1)
   ))

   # Función mejorada para crear y dibujar redes
   create_and_plot_network <- function(terms_data, title) {
       n <- nrow(terms_data)
       if (n < 2) return(NULL)
       
       # Crear matriz de similitud
       similarity_matrix <- matrix(0, nrow = n, ncol = n)
       rownames(similarity_matrix) <- terms_data\$term_name
       colnames(similarity_matrix) <- terms_data\$term_name
       
       # Calcular similitudes de manera segura
       for(i in 1:n) {
           for(j in 1:n) {
               if (i != j) {
                   similarity_matrix[i,j] <- min(
                       terms_data\$intersection_size[i] / terms_data\$term_size[i],
                       terms_data\$intersection_size[j] / terms_data\$term_size[j]
                   )
               }
           }
       }
       
       # Crear grafo con nueva sintaxis
       graph <- graph_from_adjacency_matrix(
           similarity_matrix,
           weighted = TRUE,
           mode = "undirected",
           diag = FALSE
       )
       
       # Dibujar grafo con manejo de errores
       tryCatch({
           plot(graph,
               vertex.size = sqrt(pmax(1, terms_data\$intersection_size))*3,
               vertex.color = colorRampPalette(c("blue", "red"))(n),
               vertex.label.cex = 0.6,
               edge.width = E(graph)\$weight*2,
               layout = layout_with_fr,
               main = title)
       }, error = function(e) {
           message("No se pudo crear el gráfico de red para: ", title)
           message("Error: ", e\$message)
       })
   }

   # Leer y procesar datos
   deseq_data <- read.csv("${deseq_results}")
   deseq_data\$gene_id_clean <- gsub("\\\\.[0-9]+\$", "", deseq_data\$gene_id)

   # Filtrar genes significativos con manejo de NA
   sig_genes_up <- deseq_data %>%
       filter(!is.na(padj) & !is.na(log2FoldChange)) %>%
       filter(padj < ${params.gprofiler_user_threshold} & log2FoldChange >= 1) %>%
       pull(gene_id_clean)

   sig_genes_down <- deseq_data %>%
       filter(!is.na(padj) & !is.na(log2FoldChange)) %>%
       filter(padj < ${params.gprofiler_user_threshold} & log2FoldChange <= -1) %>%
       pull(gene_id_clean)

   sig_genes <- unique(c(sig_genes_up, sig_genes_down))

   # Análisis de enriquecimiento 
   gost_results <- gost(
       query = sig_genes,
       organism = "${params.gprofiler_organism}",
       sources = strsplit("${params.gprofiler_sources}", ",")[[1]],
       significant = TRUE,
       user_threshold = ${params.gprofiler_user_threshold},
       correction_method = "${params.gprofiler_correction_method}",
       domain_scope = "annotated",
       custom_bg = ${params.gprofiler_background == "" ? "NULL" : params.gprofiler_background},
       ordered_query = FALSE,
       exclude_iea = TRUE
   )

   # Procesar resultados con manejo seguro de datos
   result_df <- if(!is.null(gost_results)) {
       as.data.frame(gost_results\$result) %>%
           select(-any_of(c("parents", "intersection"))) %>%
           mutate(
               enrichment_ratio = as.numeric(intersection_size) / as.numeric(query_size),
               log10_pval = -log10(as.numeric(p_value))
           )
   } else {
       data.frame() # DataFrame vacío si no hay resultados
   }

   # Guardar resultados
   write.table(
       result_df,
       file = "${meta.id}_gprofiler_results.tsv",
       sep = "\\t",
       row.names = FALSE,
       quote = FALSE
   )

   # Crear resumen detallado
   sink("${meta.id}_summary.txt")
   cat("ANÁLISIS DE ENRIQUECIMIENTO FUNCIONAL\\n")
   cat("====================================\\n\\n")
   cat("ESTADÍSTICAS GENERALES\\n")
   cat("---------------------\\n")
   cat(sprintf("Total genes analizados: %d\\n", nrow(deseq_data)))
   cat(sprintf("Genes significativos (padj < ${params.gprofiler_user_threshold}): %d\\n", 
       sum(deseq_data\$padj < ${params.gprofiler_user_threshold}, na.rm = TRUE)))
   cat(sprintf("Genes up-regulados (log2FC >= 1): %d\\n", length(sig_genes_up)))
   cat(sprintf("Genes down-regulados (log2FC <= -1): %d\\n", length(sig_genes_down)))
   cat(sprintf("Total genes para enriquecimiento: %d\\n\\n", length(sig_genes)))
   cat("RESULTADOS DEL ANÁLISIS\\n")
   cat("----------------------\\n")
   cat(sprintf("Total términos enriquecidos: %d\\n\\n", nrow(result_df)))
   if(nrow(result_df) > 0) {
       cat("Distribución por categoría:\\n")
       print(table(result_df\$source))
   }
   sink()

   # Crear visualizaciones
   pdf("${meta.id}_enrichment_plots.pdf", width = 12, height = 10)

   # Volcano plot con manejo de NA
   p1 <- ggplot(deseq_data %>% filter(!is.na(padj) & !is.na(log2FoldChange)),
       aes(x = log2FoldChange, y = -log10(padj))) +
       geom_point(aes(color = padj < ${params.gprofiler_user_threshold} & 
                   (log2FoldChange >= 1 | log2FoldChange <= -1)),
               alpha = 0.6, size = 1) +
       scale_color_manual(values = c("grey70", "#E41A1C"),
                       labels = c("No significativo", "Significativo")) +
       geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue", alpha = 0.5) +
       geom_hline(yintercept = -log10(${params.gprofiler_user_threshold}), 
               linetype = "dashed", color = "blue", alpha = 0.5) +
       labs(
           title = "Volcano Plot",
           subtitle = sprintf("%d genes significativamente expresados (|log2FC| >= 1)", 
                           length(sig_genes)),
           x = "log2 Fold Change",
           y = "-log10(adjusted p-value)",
           color = "Estado"
       ) +
       theme(legend.position = "bottom")
   print(p1)

   if(nrow(result_df) > 0) {
       # GeneRatio plot
       result_df_with_direction <- result_df %>%
           mutate(
               GeneRatio = intersection_size / term_size,
               log2_enrichment = log2(enrichment_ratio),
               signed_ratio = case_when(
                   log2_enrichment >= 0 ~ GeneRatio,
                   log2_enrichment < 0 ~ -GeneRatio,
                   TRUE ~ 0
               ),
               Direction = case_when(
                   log2_enrichment >= 0 ~ "upregulated (FC >= 1)",
                   log2_enrichment < 0 ~ "downregulated (FC <= -1)",
                   TRUE ~ "unchanged"
               )
           ) %>%
           mutate(Direction = factor(Direction, 
                  levels = c("upregulated (FC >= 1)", "downregulated (FC <= -1)", "unchanged")))

       # Plots adicionales solo si hay resultados
       if(nrow(result_df_with_direction) > 0) {
           # Tomar los primeros 15 términos por categoría
           top_terms <- result_df_with_direction %>%
               group_by(source) %>%
               slice_max(order_by = abs(signed_ratio), n = 15) %>%
               ungroup()

           p2 <- ggplot(top_terms, aes(x = signed_ratio, y = reorder(term_name, abs(signed_ratio)))) +
               geom_point(aes(size = intersection_size, color = p_value)) +
               scale_color_viridis_c(
                   option = "plasma",
                   name = "p.adjust",
                   direction = -1
               ) +
               scale_size_continuous(name = "Count", range = c(2, 8)) +
               scale_x_continuous(
                   breaks = seq(-1, 1, by = 0.2),
                   limits = function(x) c(min(x) * 1.1, max(x) * 1.1)
               ) +
               geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", alpha = 0.5) +
               labs(
                   title = "Enriched Pathways",
                   x = "GeneRatio (negative values = down-regulation)",
                   y = NULL
               ) +
               theme(
                   axis.text.y = element_text(size = 8),
                   axis.text.x = element_text(angle = 45, hjust = 1),
                   plot.title = element_text(size = 12, face = "bold"),
                   legend.position = "right"
               )
           print(p2)

           # Manhattan plot
           print(gostplot(gost_results, 
                       interactive = FALSE, 
                       capped = TRUE) + 
               ggtitle("Manhattan Plot de Términos Enriquecidos"))

           # Ridge plot
           source_counts <- table(result_df\$source)
           sources_with_enough_data <- names(source_counts[source_counts >= 3])
           
           if(length(sources_with_enough_data) > 0) {
               filtered_df <- result_df %>%
                   filter(source %in% sources_with_enough_data)
               
               p3 <- ggplot(filtered_df, 
                         aes(x = log10_pval, y = source, fill = source)) +
                   geom_density_ridges(alpha = 0.6) +
                   scale_fill_brewer(palette = "Set3") +
                   labs(
                       title = "Enrichment Distribution by Source",
                       x = "-log10(p-value)",
                       y = "Source"
                   ) +
                   theme_ridges()
               print(p3)
           }

           # Dot plot por categoría
           top_10_by_source <- result_df %>%
               group_by(source) %>%
               slice_max(order_by = log10_pval, n = 10) %>%
               ungroup()

           p4 <- ggplot(top_10_by_source, 
                     aes(x = enrichment_ratio, y = reorder(term_name, log10_pval))) +
               geom_point(aes(size = as.numeric(intersection_size), 
                           color = log10_pval)) +
               facet_wrap(~source, scales = "free_y", ncol = 1) +
               scale_color_viridis_c(option = "plasma") +
               scale_size_continuous(range = c(2, 10)) +
               labs(
                   title = "Top 10 Términos Enriquecidos por Categoría",
                   x = "Ratio de Enriquecimiento",
                   y = NULL,
                   size = "Genes",
                   color = "-log10(p-value)"
               ) +
               theme(
                   axis.text.y = element_text(size = 8),
                   strip.background = element_rect(fill = "grey90"),
                   strip.text = element_text(face = "bold")
               )
           print(p4)

           # Networks para GO:BP, GO:CC y GO:MF
           for (go_type in c("GO:BP", "GO:CC", "GO:MF")) {
               if (go_type %in% result_df\$source) {
                   terms <- result_df %>%
                       filter(source == go_type) %>%
                       slice_max(n = 30, order_by = log10_pval)
                   
                   if (nrow(terms) >= 2) {
                       create_and_plot_network(terms, paste(go_type, "Term Network"))
                   }
               }
           }
       }
   }

   dev.off()

   # Crear reporte HTML
   sink("${meta.id}_enrichment_report.html")
   cat('<!DOCTYPE html>
   <html>
   <head>
       <title>Reporte de Enriquecimiento Funcional</title>
       <style>
           body { font-family: Arial, sans-serif; margin: 40px; }
           h1, h2 { color: #2c3e50; }
           .stats { background-color: #f7f9fc; padding: 20px; border-radius: 5px; }
           table { border-collapse: collapse; width: 100%; margin: 20px 0; }
           th, td { padding: 8px; text-align: left; border-bottom: 1px solid #ddd; }
           th { background-color: #f2f2f2; }
       </style>
   </head>
   <body>')
   
   cat("<h1>Reporte de Enriquecimiento Funcional</h1>")
   cat("<div class='stats'>")
   cat("<h2>Estadísticas Generales</h2>")
   cat(sprintf("<p>Total genes analizados: %d</p>", nrow(deseq_data)))
   cat(sprintf("<p>Genes up-regulados (log2FC >= 1): %d</p>", length(sig_genes_up)))
   cat(sprintf("<p>Genes down-regulados (log2FC <= -1): %d</p>", length(sig_genes_down)))
   cat(sprintf("<p>Términos enriquecidos: %d</p>", nrow(result_df)))
   cat("</div>")
   
   if(nrow(result_df) > 0) {
       cat("<h2>Top 20 Términos Enriquecidos</h2>")
       cat("<table>")
       cat("<tr><th>Término</th><th>Categoría</th><th>p-value</th><th>Genes</th><th>Ratio</th></tr>")
       
       top_20_terms <- result_df %>%
           arrange(as.numeric(p_value)) %>%
           head(20)
       
       for(i in 1:nrow(top_20_terms)) {
           cat(sprintf("<tr><td>%s</td><td>%s</td><td>%.2e</td><td>%s</td><td>%.2f</td></tr>\\n",
               top_20_terms\$term_name[i],
               top_20_terms\$source[i],
               as.numeric(top_20_terms\$p_value[i]),
               top_20_terms\$intersection_size[i],
               top_20_terms\$enrichment_ratio[i]))
       }
       
       cat("</table>")
   } else {
       cat("<p>No se encontraron términos enriquecidos significativos.</p>")
   }
   
   cat("</body></html>")
   sink()

   # Guardar versiones
   writeLines(
       c(
           "---",
           "gprofiler2:",
           paste0("    version: ", packageVersion("gprofiler2")),
           "ggplot2:",
           paste0("    version: ", packageVersion("ggplot2")),
           "tidyverse:",
           paste0("    version: ", packageVersion("tidyverse"))
       ),
       "versions.yml"
   )
   """
}
