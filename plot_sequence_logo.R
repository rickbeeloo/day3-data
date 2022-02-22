library(bioseq)
library(insect)
library(ggplot2)
library(cowplot) # Just to plot two logos in the same plot

# Load the contig FASTA sequences
# we can then use substr() to get parts of the contig
contig.seq <- bioseq::read_fasta('bin_01.fasta')

# Function to plot sequence logo for prodigal predictions
plot.logo <- function(df) {
  # Get the sequences for normal mode
  seqs = c()
  for (i in 1:nrow(df)) {
    row = df[i,]
    
    # Get the upstream seq
    if (row$strand == '+') {
      upstream_seq = substr(contig.seq, row$start-30, row$start+2)
    
    # Get the upstream sequence for the - strand, in this case
    # we also have to take the reverse complement (insect:rc)
    } else {
      end = row$end
      upstream_seq = insect::rc(substr(contig.seq,end-2, end + 30))
    }
    
    # Check if we have the full sequence, for example it could be just
    # on the end of a contig and hence not 33 nucleotides, we will exclude these 
    if (nchar(upstream_seq) == 33) {
      seqs = append(seqs, upstream_seq)
    }
  }
  
  ggseqlogo::ggseqlogo(seqs)
}

p1 <- plot.logo(default.genes) + ggtitle('Normal')
p2 <- plot.logo(meta.genes) + ggtitle('Meta')
cowplot::plot_grid(p1, p2, labels = 'auto')