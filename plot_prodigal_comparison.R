library(gggenes)
library(dplyr)

setwd('C:/Users/Rick/course_output')

# function to parse prodigal annotation file
parse.prodigal <- function(f) {
  print(f)
  read.table(f) %>%
    mutate(V1 = gsub('>', '', V1)) %>%
    tidyr::separate(V1, sep = '_', into = c('id' ,'start', 'end', 'strand')) %>%
    mutate(start = as.numeric(start), end = as.numeric(end))
}

# Load the default annotation
default.genes <- parse.prodigal('genes_default.txt')
default.genes$group <- 'default'

# Load the meta annotations
meta.genes <- parse.prodigal('genes_meta.txt')
meta.genes$group <- 'meta'

# Combine the data
genes <- rbind(default.genes, meta.genes)

# Check when they predict the same
combined <- genes %>% 
  rowwise() %>%
  mutate(tag = paste(start, end, strand, sep = '_')) %>%
  group_by(tag) %>%
  summarise(start = unique(start), end = unique(end), strand = unique(strand), groups = paste0(sort(unique(group)), collapse = ','))

# Plot this
ggplot(combined, aes(xmin = start , xmax = end, y = groups, fill = groups)) +
  geom_gene_arrow() +
  facet_wrap(~ groups, scales = "free", ncol = 1) +
  scale_fill_brewer(palette = "Set3") +
  theme_genes()








