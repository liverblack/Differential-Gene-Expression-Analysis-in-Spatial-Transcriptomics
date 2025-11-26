library(DESeq2)
library(ggplot2)
library(ggrepel)
library(dplyr)
raw <- read.csv("D:/bioinfo/GSM5591748_df.csv", header = TRUE)
rownames(raw) <- paste0("S",1:nrow(raw))
cts <- raw[, -ncol(raw)]
coldata <- data.frame(layer = raw[, ncol(raw)])
rownames(coldata) <- rownames(raw)
cts <- t(cts) 
coldata$layer <- factor(coldata$layer)
dds <- DESeqDataSetFromMatrix(
  countData = cts,
  colData = coldata,
  design = ~ layer
)
colData(dds)$layer <- factor(
  colData(dds)$layer,
  levels = c("0", "1", "2", "3")
)
dds <- DESeq(dds)

contrast_pairs <- list(
  c("layer", "1", "0"),
  c("layer", "2", "0"),
  c("layer", "3", "0"),
  c("layer", "1", "2"),
  c("layer", "2", "3"),
  c("layer", "1", "3")
)

results_list <- lapply(contrast_pairs, function(pair) {
  res <- results(dds, contrast = pair)
  comparison_name <- paste(pair[2], "vs", pair[3])
  list(res = res, name = comparison_name)
})

res_combined <- bind_rows(lapply(results_list, function(x) {
  res <- as.data.frame(x$res)
  res %>%
    mutate(
      comparison = x$name,  
      gene = rownames(res),  
      significant = ifelse(
        padj < 0.05 & abs(log2FoldChange) > 1,
        "Significant",
        "Not Significant"
      ),
      direction = case_when(
        significant == "Significant" & log2FoldChange > 0 ~ "Up",
        significant == "Significant" & log2FoldChange < 0 ~ "Down",
        TRUE ~ "Not Significant"
      )
    ) %>%
    na.omit()  
}))

ggplot(res_combined, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(
    aes(color = direction, alpha = direction),
    size = 1.5
  ) +
  scale_color_manual(
    name = "Expression",
    values = c(
      "Up" = "#E64B35",     # 显著上调（红色）
      "Down" = "#3C5488",  # 显著下调（蓝色）
      "Not Significant" = "#999999"  # 非显著（灰色）
    )
  ) +
  scale_alpha_manual(
    values = c(
      "Up" = 0.8,
      "Down" = 0.8,
      "Not Significant" = 0.3
    ),
    guide = "none"  
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_text_repel(
    data = res_combined %>%
      group_by(comparison) %>%
      filter(direction != "Not Significant") %>%
      top_n(-5, padj),  # 每个对比标注前5个显著基因
    aes(label = gene),
    size = 3,
    max.overlaps = 20
  ) +
  facet_wrap(
    ~ comparison,
    ncol = 3,       
    scales = "free"     
  ) +
  theme_bw() +
  labs(
    x = "log2(Fold Change)",
    y = "-log10(Adjusted p-value)",
    title = "Combined Volcano Plot (All Comparisons)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom"
  )