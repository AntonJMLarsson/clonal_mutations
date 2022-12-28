library(ggplot2)
library(circlize)
prefix <- "PBMC_donor7"
mgatk_se <- readRDS("mgatk_out/final/mgatk.rds")
mgatk_se_filter <- mgatk_se[, mgatk_se$sample[mgatk_se$depth >= 10]]

source("/home/antonl/programs/mgatk_scripts/variant_calling.R") 

# mgatk_se is the Summarized Experiment .rds file
# That is automatically produced from running
# The mgatk CLI python package 

mut_se <- call_mutations_mgatk(mgatk_se_filter)

misc_df <- data.frame(rowData(mut_se))


#circos.genomicLines(cov_per_pos_table[,c('loc','locplus1')], cov_per_pos_table$cov_per_position)
cov <- assays(mgatk_se)[["coverage"]]
cov_per_position = Matrix::rowMeans(cov)
cov_per_pos_table <- data.table(cov_per_position)
cov_per_pos_table$chrom <- "MT"
cov_per_pos_table$loc <- rownames(cov_per_pos_table)
cov_per_pos_table$loc <- as.integer(cov_per_pos_table$loc)
cov_per_pos_table$locplus1 <- as.integer(cov_per_pos_table$loc) + 1
cov_per_pos_table$cov_per_position <- log2(cov_per_pos_table$cov_per_position)
cov_per_pos_table <- cov_per_pos_table[,c('chrom','loc','locplus1','cov_per_position')]
cov_per_pos_table <- cov_per_pos_table[,c('chrom','loc','cov_per_position')]

write.csv(cov_per_pos_table, file = paste(prefix, 'mito_coverage.csv', sep='_'))

pdf(paste(prefix, "MT_coverage.pdf", sep='_'))
circos.initializeCircularGenome('MT', 16569, major.by = 1000)
circos.genomicTrack(cov_per_pos_table, panel.fun = function(region, value, ...) {
  i = getI(...)
  circos.genomicLines(region, value, col = i, ...)
})
text(0, 0, "MT", cex = 1)
dev.off()

print(mean(2**cov_per_pos_table$cov_per_position))

count_df <- data.frame(as.matrix(t(assays(mut_se)[["count"]])))
coverage_df <- data.frame(as.matrix(t(assays(mut_se)[["coverage"]])))
mito_df <- data.frame(as.matrix(t(assays(mut_se)[["allele_frequency"]])))

write.csv(mito_df, file=paste(prefix, "mitochondrial_allele_frequency_all_variants.csv", sep='_'))
write.csv(coverage_df, file=paste(prefix, "mitochondrial_coverage_all_variants.csv", sep='_'))
write.csv(count_df, file=paste(prefix, "mitochondrial_count_all_variants.csv", sep='_'))
write.csv(misc_df, file=paste(prefix, "mutation_summary_statistics.csv", sep='_'))