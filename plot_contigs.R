# Load the packages
library(dplyr)
library(gggenes)
library(ggplot2)
library(gridExtra)


###################################################################
# Simple function to escape || in a string such that filtering 
# (e.g. grepl) does not interpret it as a RegEx
###################################################################
fix.name <- function(name){return(gsub('\\|\\|', '\\\\|\\\\|', name))}

###################################################################
# Subsets dataframe for a specific contig identifier
###################################################################
get.contig <- function(contig.id) {
  contig <- df %>% filter(grepl(fix.name(contig.id), df$query)) %>%
    mutate(orientation= ifelse(strand == '+', TRUE, FALSE)) %>%
    mutate(molecule = contig.id)
}

###################################################################
# Plot genes on one or more molecules
###################################################################
plot.genes <- function(contigs) {
  ggplot(contigs, aes(xmin = start, xmax = end, y = molecule, fill = annot,  forward = orientation)) +
    geom_gene_arrow() + 
    facet_wrap(molecule ~ 'genome', scales = "free", ncol = 1) +
    theme_genes() +
    guides(fill="none")
}

###################################################################
# Calculates length of the predicted proteins in each group
# and plots the density distribution of those
###################################################################
plot.lengths <- function(contigs) {
  lens <- contigs %>% mutate(len = end - start)
  names <- unique(lens$molecule)
  mean1 <- round(lens %>% filter(molecule == names[1]) %>% pull(len) %>% mean(),0)
  mean2 <- round(lens %>% filter(molecule == names[2]) %>% pull(len) %>% mean(),0)
  label.pos <- 0.0018
   ggplot(lens, aes(x=len, fill = molecule)) + geom_density(alpha = 0.4) +
    theme_bw() +
    geom_vline(xintercept = mean1, size = 0.3, col = 'red', alpha = 0.8) +
    geom_vline(xintercept = mean2, size = 0.3, col = 'blue', alpha = 0.8) +
    geom_text(aes(x=mean1, y = label.pos, label = mean1), hjust= -0.5, col = 'red') +
    geom_text(aes(x=mean2, y = label.pos, label = mean2), hjust= -0.5, col = 'blue') +
    scale_fill_manual(values = c('red', 'blue')) +
    xlab('Predicted protein length') +
    ylab('Density')
}

###################################################################
# Searches "annotations" that are split in multiple predicted 
# proteins (in one or both contigs). Then it plots the number of 
# predicted proteins per contig
###################################################################
plot.anno <- function(contigs) {
  grouped.anno <- contigs %>% 
    filter(annot != '') %>%
    group_by(molecule, annot) %>% 
    summarise(n.proteins = n_distinct(start))
  multi.anno <- grouped.anno %>% filter(n.proteins > 1) %>% head(10) # max 10
  grouped.anno <- grouped.anno %>% filter(annot %in% multi.anno$annot)
  
  ggplot(grouped.anno, aes(x=annot, y = n.proteins, fill = molecule)) +
    geom_bar(position="dodge", stat="identity", alpha = 0.4) +
    theme_bw() +
    coord_flip() +
    ylab('Number of predicted ORFs') +
    xlab('Annotation') +
    scale_fill_manual(values = c('red', 'blue')) 
}


# Ask user for the hhsearch results file
df <- read.table(file.choose(), sep = '\t', header = T, fill = T, comment.char = '')


# Infinite loop to ask user to provide contig identifiers
while (TRUE) {
  # Ask for user input 
  contig.id.1 <- readline("First contig id: ")  
  contig.id.2 <- readline("Second contig id: ") 
  out.name    <- readline("Output path (end with.pdf): ")
  
  # Run contigs from dataframe
  contig1 <- get.contig(contig.id.1)
  contig2 <- get.contig(contig.id.2)
  contigs <- rbind(contig1, contig2)
  
  # Produce plots
  genes   <- plot.genes(contigs)
  lengths <- plot.lengths(contigs)
  annots  <- plot.anno(contigs) 
  
  # Actually plot it
  grid.arrange(genes, lengths, annots, layout_matrix = rbind(c(1,1), c(2,3)))
  
  # Also save it just in case of weird symbols again
  pdf(out.name, width = 12, height = 8) 
  grid.arrange(genes, lengths, annots, layout_matrix = rbind(c(1,1), c(2,3)))
  dev.off() 
  
  print('Done!')
  continue = readline("Continue (Y/N): ")
  if (continue == 'N') {
    break
  }
}
















