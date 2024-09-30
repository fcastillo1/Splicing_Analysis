---
title: "Análisis de Expresión Diferencial en Splicing"
author: Francisca Reyes 
date: Agosto 2024
output: 
  flexdashboard::flex_dashboard:
    orientation: columns
    vertical_layout: fill
    theme: flatly
params:
  deseq2_results: "deseq2_results.csv"
  normalized_counts: "normalized_counts.csv"
  samplesheet: "samplesheet.csv"
  gtf_file: "annotation.gtf"
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)

library(flexdashboard)
library(ggplot2)
library(plotly)
library(DT)
library(dplyr)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(scales)
library(reshape2)
library(knitr)
library(dendextend)
library(htmlwidgets)
library(kableExtra)

# Definición de la función safe_read
safe_read <- function(file_path, ...) {
  tryCatch(
    {
      read.csv(file_path, ...)
    },
    error = function(e) {
      message(paste("Error al leer", file_path, ":", e$message))
      return(NULL)
    }
  )
}

# Cargar datos
deseq2_results <- safe_read(params$deseq2_results, row.names = 1)
if (is.null(deseq2_results)) {
  stop("No se pudo cargar el archivo de resultados de DESeq2")
}

normalized_counts <- safe_read(params$normalized_counts, row.names = 1)
if (is.null(normalized_counts)) {
  warning("No se pudo cargar el archivo de conteos normalizados")
} else {
  # Asegurarse de que los datos son numéricos
  normalized_counts <- apply(normalized_counts, 2, as.numeric)
}

samplesheet <- read.csv(params$samplesheet)

# Definir colores
color_palette <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")

# Definir la función format_pvalue
format_pvalue <- function(p) {
  ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
}
```

Resumen
===================================== 

Column {data-width=550}
-------------------------------------

### Descripción del Análisis

Este dashboard presenta los resultados del análisis de expresión diferencial realizado con DESeq2. El análisis compara los niveles de expresión génica entre dos condiciones, identificando genes que muestran cambios significativos en su expresión.

```{r}
sig_genes <- sum(deseq2_results$padj < 0.05, na.rm = TRUE)
up_reg <- sum(deseq2_results$padj < 0.05 & deseq2_results$log2FoldChange > 0.58, na.rm = TRUE)
down_reg <- sum(deseq2_results$padj < 0.05 & deseq2_results$log2FoldChange < 0.58, na.rm = TRUE)

# Generar el texto descriptivo
main_findings <- paste(
  "Se analizaron un total de", nrow(deseq2_results), "genes.",
  "Se identificaron", sig_genes, "genes diferencialmente expresados (padj < 0.05).",
  "De estos,", up_reg, "están sobre-expresados y", down_reg, "sub-expresados en la condición de tratamiento comparada con el control."
)

# Mostrar el texto
cat(main_findings)
```

### Tabla de datos analizados
```{r}
# Leer el archivo CSV (ajustar la ruta según sea necesario)
datos <- read.csv(params$samplesheet, sep = ",")

# Crear una tabla con solo las columnas sample_id y condition
tabla_resultado <- datos %>%
  select(sample_id, condition)

# Mostrar la tabla con formato kable
knitr::kable(tabla_resultado, caption = "Tabla de Identificadores de Muestra y Condiciones", 
             col.names = c("ID de Muestra", "Condición"), 
             format = "html", 
             table.attr = "class='table table-striped'")
``` 
-------------------------------------


Gráficos de Diagnóstico
=====================================   

Column {.tabset}
-------------------------------------

### MA Plot

```{r}
ma_plot <- ggplot(deseq2_results, aes(baseMean, log2FoldChange, 
                  color = factor(case_when(
                    is.na(padj) ~ "NA",
                    padj >= 0.05 ~ "No significativo",
                    padj < 0.05 & log2FoldChange > 0.58 ~ "Up-regulated",
                    padj < 0.05 & log2FoldChange < 0.58 ~ "Down-regulated"
                  )))) +
  geom_point(size = 1, alpha = 0.6) +
  scale_x_log10(
    labels = scales::trans_format("log10", scales::math_format(10^.x)),
    breaks = scales::trans_breaks("log10", function(x) 10^x)
  ) +
  geom_hline(yintercept = 0, color = "blue", linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = c(-1, 1), color = "blue", linetype = "dotted", size = 0.5) +
  labs(x = "Media de la expresión normalizada", 
       y = "Log2 Fold Change", 
       color = "Significancia",
       title = "MA Plot") +
  scale_color_manual(
    values = c("No significativo" = "grey50", 
               "Up-regulated" = "#00BFC4", 
               "Down-regulated" = "#F8766D", 
               "NA" = "black"),
    labels = c("No significativo", 
               "Up-regulated (padj < 0.05)", 
               "Down-regulated (padj < 0.05)", 
               "NA",
               "Umbral de Fold Change")
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
    axis.title = element_text(face = "bold", size = 16),
    legend.title = element_text(face = "bold", size = 14),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_line(color = "grey95")
  ) +
  scale_y_continuous(limits = c(-10, 10), breaks = seq(-10, 10, by = 5))

# Añadir etiquetas para los genes más significativos
top_genes <- deseq2_results %>%
  filter(padj < 0.05) %>%
  top_n(10, wt = abs(log2FoldChange))

ma_plot <- ma_plot +
  geom_text_repel(data = top_genes, 
                  aes(label = gene_name), 
                  size = 3, 
                  box.padding = 0.5, 
                  point.padding = 0.5,
                  force = 1,
                  max.overlaps = Inf)

ggplotly(ma_plot)
```

### Volcano Plot

```{r}
deseq2_results$category <- "Not Significant"
deseq2_results$category[deseq2_results$padj < 0.05 & deseq2_results$log2FoldChange > 0.58] <- "Upregulated"
deseq2_results$category[deseq2_results$padj < 0.05 & deseq2_results$log2FoldChange < -0.58] <- "Downregulated"

# Asumiendo que deseq2_results tiene las columnas necesarias
volcano_data <- deseq2_results %>%
  mutate(
    significance = case_when(
      padj < 0.05 & abs(log2FoldChange) > 0.58 ~ "Significant",
      TRUE ~ "Not Significant"
    ),
    label = ifelse(padj < 0.05 & abs(log2FoldChange) > 2, gene_name, "")  # Etiquetas para los genes más significativos
  )

volcano_plot <- ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj), color = category, 
                                         text = paste("Gene:", gene_name, 
                                                      "<br>log2FoldChange:", round(log2FoldChange, 2),
                                                      "<br>padj:", signif(padj, 3)))) +
  geom_point(alpha = 0.6, size = 1) +
  geom_vline(xintercept = c(-1, 1), color = "#F8766D", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), color = "#00BFC4", linetype = "dashed") +
  geom_text_repel(aes(label = label), box.padding = 0.5, max.overlaps = 20) +
  labs(x = "Log2 Fold Change", y = "-Log10 Adjusted p-value", color = "Category") +
  scale_color_manual(values = c( "Upregulated" = "#00BFC4", "Downregulated" = "#F8766D", "Not Significant" = "grey50"),
                     labels = c("Not Significant", "Upregulated", "Downregulated")) +
  theme_minimal() +
  theme(legend.position = "bottom")

# Convertir a plotly y personalizar la interactividad
volcano_plotly <- ggplotly(volcano_plot, tooltip = "text") %>%
  layout(legend = list(orientation = "h", y = -0.2)) %>%
  # Añadir selección de puntos
  highlight(on = "plotly_selected", off = "plotly_deselect", 
            selected = attrs_selected(marker = list(size = 10)))

# Mostrar el gráfico interactivo
volcano_plotly
```

### PCA Plot

```{r}
# Leer y preparar los datos
normalized_counts <- read.csv(params$normalized_counts, row.names = 1)
conditions <- read.csv(params$samplesheet)$condition

# Función para simular la transformación vst
vst_transform <- function(x) {
  sqrt(x + 3/8)
}

# Aplicar la transformación vst a los datos normalizados
vst_counts <- vst_transform(normalized_counts)

# Realizar PCA
pca_result <- prcomp(t(vst_counts), center = TRUE, scale. = FALSE)
pcaData <- as.data.frame(pca_result$x)
pcaData$condition <- factor(conditions)

# Calcular el porcentaje de varianza explicada
percentVar <- round(100 * (pca_result$sdev^2 / sum(pca_result$sdev^2)), 1)

# Determinar el número de condiciones únicas
unique_conditions <- levels(pcaData$condition)
n_conditions <- length(unique_conditions)

# Generar colores y formas para cada condición
condition_colors <- colorRampPalette(brewer.pal(9, "Set1"))(n_conditions)
condition_shapes <- rep(c(16, 17, 15, 18, 8, 13), length.out = n_conditions)

# Crear el gráfico
pca_plot <- ggplot(pcaData, aes(x = PC1, y = PC2, color = condition, shape = condition)) +
  geom_point(size = 3) +
  scale_color_manual(values = setNames(condition_colors, unique_conditions)) +
  scale_shape_manual(values = setNames(condition_shapes, unique_conditions)) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("Principal Component Analysis") +
  theme_bw() +
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_line(color = "grey95"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.title = element_text(face = "bold"),
    legend.position = "right",
    aspect.ratio = 1
  )

# Convertir a gráfico interactivo
pca_interactive <- ggplotly(pca_plot)
pca_interactive
```

### Heatmap

```{r}
if (!is.null(normalized_counts) && !is.null(deseq2_results)) {
  # Selecciona los 40 genes más significativos
  top_genes <- head(order(deseq2_results$padj), 40)
  mat <- normalized_counts[top_genes, ]
  
  # Log-transforma y escala los datos
  mat <- log2(mat + 1)
  mat <- t(scale(t(mat)))
  
  # Usar gene_name en lugar de gene_id
  rownames(mat) <- deseq2_results$gene_name[top_genes]
  
  # Realizar clustering jerárquico para genes y muestras
  hc_genes <- hclust(dist(mat))
  hc_samples <- hclust(dist(t(mat)))
  
  # Reordenar la matriz basada en el clustering
  mat_ordered <- mat[hc_genes$order, hc_samples$order]
  
  # Preparar los datos para plotly
  melted_mat <- reshape2::melt(mat_ordered)
  colnames(melted_mat) <- c("Gene", "Sample", "Expression")
  
  # Crear el heatmap interactivo con plotly
  heatmap_plotly <- plot_ly(
    data = melted_mat,
    x = ~Sample,
    y = ~Gene,
    z = ~Expression,
    type = "heatmap",
    colors = colorRamp(c("navy", "white", "firebrick3")),
    hoverinfo = "text",
    text = ~paste("Gene:", Gene, "<br>Sample:", Sample, "<br>Expression:", round(Expression, 2))
  ) %>%
  layout(
    title = "Top 40 Differentially Expressed Genes",
    xaxis = list(title = "Samples", tickangle = 40),
    yaxis = list(title = "Genes"),
    margin = list(l = 100, r = 20, b = 100, t = 50)
  )
  
  # Mostrar el heatmap interactivo
  heatmap_plotly
} else {
  cat("No se pudo generar el heatmap debido a la falta de datos necesarios.")
}
```

### Dispersion Plot

```{r}
# Crear gene_stats si no existe
if (!exists("gene_stats")) {
  
  # Asegurarse de que deseq2_results y normalized_counts existen
  if (exists("deseq2_results") && exists("normalized_counts")) {
    common_genes <- intersect(rownames(deseq2_results), rownames(normalized_counts))
    
    gene_stats <- data.frame(
      baseMean = rowMeans(normalized_counts[common_genes, ]),
      padj = deseq2_results[common_genes, "padj"],
      row.names = common_genes
    )
    
    # Calcular la dispersión
    gene_stats$dispersion <- apply(normalized_counts[common_genes, ], 1, var) / gene_stats$baseMean
    
    # Eliminar filas con NA o valores infinitos
    gene_stats <- gene_stats[is.finite(gene_stats$baseMean) & 
                             is.finite(gene_stats$dispersion) &
                             !is.na(gene_stats$padj), ]
  } else {
    stop("No se pueden crear gene_stats porque faltan deseq2_results o normalized_counts")
  }
}

if (exists("gene_stats") && nrow(gene_stats) > 0) {
  
  # Calcular la línea de tendencia
  fit <- loess(log10(dispersion) ~ log10(baseMean), data = gene_stats, span = 0.3)
  gene_stats$fitted <- 10^predict(fit)
  
  # Identificar outliers (simulando dispersión final)
  gene_stats$outlier <- gene_stats$dispersion > quantile(gene_stats$dispersion, 0.99)
  
# Asumiendo que gene_stats ya está creado y contiene los datos necesarios
  dispersion_plot <- ggplot(gene_stats, aes(x = baseMean, y = dispersion)) +
    geom_point(aes(color = ifelse(padj < 0.05, "Dispersión final", "Dispersión genética")), alpha = 0.4, size = 1) +
    geom_smooth(aes(color = "Tendencia Ajustada (fitted)"), method = "loess", se = FALSE) +
    scale_x_log10(labels = trans_format("log10", math_format(10^.x)),
                  breaks = 10^seq(-1, 5, by = 2)) +
    scale_y_log10(labels = trans_format("log10", math_format(10^.x)),
                  breaks = 10^seq(-8, 1, by = 2)) +
    annotation_logticks() +
    scale_color_manual(values = c("Dispersión final" = "blue", "Dispersión genética" = "black", "Tendencia Ajustada (fitted)" = "red"),
                      name = "Type") +
    labs(x = "Media de la expresión normalizada", 
        y = "Dispersión estimada",
        title = "Estimaciones de Dispersión",
        color = "Tipo") +
    theme_bw() +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 12, face = "bold"),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10))

  # Convertir a plotly para hacerlo interactivo
  dispersion_plotly <- ggplotly(dispersion_plot, tooltip = c("x", "y", "color")) %>%
    layout(legend = list(title = list(text = "<b>Tipo</b>")))

  # Mostrar el gráfico interactivo
  dispersion_plotly
} else {
  print("No se pudo generar el gráfico de dispersión. Verifica que gene_stats existe y tiene datos.")
}
```

Resultados de expresión
=====================================   

Column {.tabset}
-------------------------------------

### Genes Diferencialmente Expresados
```{r}
datatable(deseq2_results, 
          options = list(pageLength = 15, scrollX = TRUE),
          filter = 'top')
```

### Top 50 Genes Sobre-expresados
```{r}
top_up <- deseq2_results %>%
  filter(padj < 0.05) %>%
  arrange(desc(log2FoldChange)) %>%
  head(50)

datatable(top_up, options = list(scrollX = TRUE))
```

### Top 50 Genes Sub-expresados
```{r}
top_down <- deseq2_results %>%
  filter(padj < 0.05) %>%
  arrange(log2FoldChange) %>%
  head(50)

datatable(top_down, options = list(scrollX = TRUE))
```

Análisis de Enriquecimiento Funcional con gProfiler
===================================================

Column {.tabset}
-------------------------------------

### Resumen de Resultados

```{r gprofiler_summary, fig.width=12, fig.height=8}
gprofiler_file <- "deseq2.gprofiler2.all_enriched_pathways.tsv"
if (file.exists(gprofiler_file)) {
  gprofiler_results <- read.delim(gprofiler_file, stringsAsFactors = FALSE)
  
  valueBox(nrow(gprofiler_results), 
           caption = "Total de términos enriquecidos", 
           icon = "fa-list", 
           color = "primary")
  
  valueBox(sum(gprofiler_results$p_value < 0.05), 
           caption = "Términos significativos (p < 0.05)", 
           icon = "fa-check", 
           color = "success")
  
  top_1000 <- gprofiler_results %>%
    arrange(p_value) %>%
    head(1000) %>%
    mutate(p_value_fmt = format_pvalue(p_value))
  
  DT::datatable(
    top_1000 %>% select(source, term_name, p_value_fmt, term_size, query_size, intersection_size),
    caption = "Top términos más significativamente enriquecidos",
    options = list(pageLength = 10, scrollX = TRUE),
    rownames = FALSE
  )
} else {
  cat("No se encontraron resultados de gProfiler.")
}
```

### Gráfico de Barras
```{r}
if (exists("gprofiler_results")) {
  top_terms <- head(gprofiler_results[order(gprofiler_results$p_value), ], 50)

  p <- ggplot(top_terms, aes(x = reorder(term_name, -log10(p_value)), y = -log10(p_value), fill = source)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_brewer(palette = "Set3") +
    labs(title = "Top 50 Términos Más Enriquecidos",
         x = "Término",
         y = "-log10(valor p)",
         fill = "Fuente") +
    theme_minimal(base_size = 12) +
    theme(axis.text.y = element_text(size = 8),
          legend.position = "bottom")

  ggplotly(p)
}
```

### Gráfico de Dispersión
```{r}
if (exists("gprofiler_results")) {
  p <- ggplot(gprofiler_results, aes(x = term_size, y = -log10(p_value), color = source, size = intersection_size)) +
    geom_point(alpha = 0.7) +
    scale_x_log10(labels = comma) +
    scale_size(range = c(1, 10), name = "Tamaño de intersección") +
    labs(title = "Distribución de Términos Enriquecidos",
         x = "Tamaño del término (escala logarítmica)",
         y = "-log10(valor p)",
         color = "Fuente") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "right")

  ggplotly(p)
}
```

### Gráfico de Enriquecimiento
```{r}
png_file <- "deseq2.gprofiler2.gostplot.png"
if (file.exists(png_file)) {
  knitr::include_graphics(png_file)
} else {
  cat("La imagen del gráfico de gProfiler no está disponible.")
}
```


Métodos y Parámetros
=====================================

### Descripción del Análisis

Este análisis de expresión diferencial se realizó utilizando el paquete DESeq2 en R. Los pasos principales del análisis fueron:

1. Importación de datos de conteo de genes.
2. Creación del objeto DESeqDataSet.
3. Estimación de factores de tamaño.
4. Estimación de dispersión.
5. Test de Wald para expresión diferencial.
6. Corrección de valores p por pruebas múltiples.

### Parámetros Utilizados

- Umbral de fold change: 0.58 (en escala log2)
- Umbral de valor p ajustado: 0.05
- Método de ajuste de valores p: Benjamini-Hochberg (FDR)

### Software y Versiones

```{r}
sessionInfo()
```
