source("./functions/simulate_sequence.R")

library(Boruta)

# define alphabet
alph <- as.character(1L:4)

n_seq = 1000
l_seq = 8

motifs <- generate_motif(alph, 5)
dat <- generate_seqs(n_seq, l_seq, motifs)

Boruta(dat, c(rep(1, n_seq), rep(0, n_seq)))
