library(ggplot2)
library(dplyr)
library(cowplot)

# Read the cluster file 
cluster.anno <- read.table(file.choose(), comment.char = '', sep = '\t', header = T, na.strings = "")

# Specify output name
out.file <- 'cluster_sizes.pdf'

# Count the number of (our) proteins per cluster
our.proteins <- cluster.anno %>% 
  group_by(reference) %>%
  summarise(
    n.proteins = n_distinct(member), 
    n.prodigal.proteins = sum(grepl('bin_', member)),
    anno = unique(annot),
    prodigal = any(grepl('bin_', member))
  ) %>% 
  unique() %>% 
  mutate(row = row_number()) %>%
  rowwise() %>%
  mutate(anno = paste(row, anno, collapse = " "))

# plot i
with.refseq <- ggplot(our.proteins %>% filter(n.proteins > 100), 
                      aes(
                        x= reorder(anno, n.proteins), 
                        y = n.proteins)
) +
  geom_bar(stat='identity', fill = 'blue', alpha = 0.5)+
  coord_flip() + 
  theme_bw() +
  xlab('#proteins') +
  ylab('PHROG annotation') +
  ggtitle('Number of proteins in cluster (with RefSeq)') +
  theme(plot.title = element_text(hjust = 0.5))


without.refseq <- ggplot(our.proteins %>%  filter(prodigal == T & n.prodigal.proteins > 5), 
                         aes(x=reorder(anno, n.prodigal.proteins), 
                             y=n.prodigal.proteins)
) +
  geom_bar(fill='red', stat='identity', alpha = 0.5) +
  coord_flip() +
  theme_bw() +
  xlab('#proteins') +
  ylab('PHROG annotation') +
  ggtitle('Number of proteins in cluster (without RefSeq)') +
  theme(plot.title = element_text(hjust = 0.5))

cowplot::plot_grid(with.refseq, without.refseq, rel_widths = c(0.6, 0.4))

pdf(out.file, width = 12, height = 8) 
cowplot::plot_grid(with.refseq, without.refseq, rel_widths = c(0.6, 0.4))
dev.off()
print(paste0("Plot saved in ", getwd()))
