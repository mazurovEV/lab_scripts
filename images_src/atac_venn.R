
source('lib.R')

library(VennDiagram)

###

fantom_df  <- read.delim(paste0(DATA_DIR, 'ENSG00000229852.fantom_two_aso.DE_Summary.txt'), as.is=TRUE, header=TRUE)
himorna_df <- read.delim(paste0(DATA_DIR, 'ENSG00000229852_H3K27ac.himorna_peaks_anno.txt'), as.is=TRUE, header=TRUE)
atac_df    <- read.delim(paste0(DATA_DIR, 'ENSG00000229852.atac_diff_peaks_anno.txt'), as.is=TRUE, header=TRUE)

head(himorna_df)

fantom_v <- fantom_df %>%
  filter(log2FC < 0) %>%
  pull(geneID) %>%
  unique()

himorna_v <- himorna_df %>%
  filter(!is.na(feature)) %>%
  # ENSG00000189337.17  =>  ENSG00000189337
  mutate(geneID = gsub('\\.\\d+$', '', feature, perl=TRUE)) %>%
  pull(geneID) %>%
  unique()

atac_v <- atac_df %>%
  filter(!is.na(feature)) %>%
  # ENSG00000189337.17  =>  ENSG00000189337
  mutate(geneID = gsub('\\.\\d+$', '', feature, perl=TRUE)) %>%
  pull(geneID) %>%
  unique()

common3_v <- Reduce(intersect, list(fantom_v, himorna_v, atac_v))

# Do not create log file: https://stackoverflow.com/a/36812214/310453
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

venn.diagram(
  x = list(fantom_v, himorna_v, atac_v),
  category.names = c(
    paste("FANTOM DEGs,\nn = ", length(fantom_v)),
    paste("HiMoRNA genes,\nn = ", length(himorna_v)),
    paste("ATAC-seq genes, n = ", length(atac_v))),
  cat.fontface = "bold",
  cat.pos = c(-15, 15, 180),
  cat.dist = c(0.05, 0.05, 0.02),
  filename = paste0(OUT_DIR, 'atac_venn.png'),
  output=TRUE
)


# Save common genes to file
# aso_07_df %>%
#   filter(geneSymbol %in% common3_v) %>%
#   select(geneSymbol) %>%
#   unique() %>%
#   write.table(file = paste0(DATA_DIR, "chaserr_common3.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

# n_07 = length(aso_07_v)
# n_10 = length(aso_10_v)
# both_aso = length(intersect(aso_07_v, aso_10_v))
# all_genes = 20000
# phyper(
#   both_aso, # number of white balls drawn without replacement from an urn
#   n_10, # the number of white balls in the urn
#   all_genes - n_10, # the number of black balls in the urn
#   n_07, # the number of balls drawn from the urn.
#   lower.tail = FALSE
# )

