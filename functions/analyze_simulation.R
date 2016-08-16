library(dplyr)
library(biogram)
library(ggplot2)


sim_quipt <- read.csv("./results/sim1.csv") %>%
  mutate(motif_len = nchar(gsub(".", "", ngrams2df(as.character(motif))[["ngram"]], fixed = TRUE))) 

# p-values of true motifs
tm_pvals <- select(sim_quipt, l_seq, n_seq, criterion, p.value, motif_len)
ggplot(tm_pvals, aes(x = criterion, y = p.value, color = as.factor(motif_len))) +
  geom_boxplot() +
  facet_grid(l_seq ~ n_seq, labeller = label_both)
